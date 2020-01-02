# The NumNomsS method for multiple testing problems.
# Based on Gladwin, Derks et al. (2012) and Gladwin & Vink, 2018.
# This is a "no smoke without fire" test that trades off precision for statistical power.
# The NumNomS set-wise test considers the number of nominally significant tests.
# A null hypothesis distribution of this number is determined by permutation testing, 
# such that correlations between variables being tested can be exploited to improve power.
# The software also provides a test-wise permutation-based familywise-corrected
# critical p-value, pPBFW, that should be at least as powerful as Bonferroni correction.
# Additionally, while this is not principled at this point, the software reports which tests
# are significant when using Bonferroni correction over the expected number of nominally significant
# tests. This is intended to be used only after protection by NumNomS significance, to reduce
# the number of false test-wise positives after a false NumNomS positive.
# Note that the number of iterations can be adjusted to address any concerns with the use
# of permutation tests to a trivial level, relative to the inherent uncertainity of sampling.

nSigFunc <- function(ToTest) {
  # Adjust this function to perform the desired tests
  # It has to return the Output list as specified
  pvec <- apply(ToTest, 2, function(v){htest<-t.test(v); htest$p.value})
  pvec <- pvec
  nSig <- length(which(pvec < .05))
  Output <- list(nSig, pvec)
}

Permed0 <- function(DF) {
  # Adjust this function to randomly permute the data as needed
  flipper0 <- matrix(runif(nObservations), ncol = 1)
  flipper0[flipper0 < 0.5] = -1
  flipper0[flipper0 >= 0.5] = 1
  flipper0 <- flipper0 %*% matrix(1, ncol = nVariables)
  ToTest <- DF * flipper0
}

NumNomS <- function(DF) {
  # Generate null hypothesis distribution of number-of-nominally-significant tests
  nIterations = 1000
  nSigvec <- c()
  smallest_p_vec <- c()
  print("Generating null hypothesis distribution...")
  for (n in 1:nIterations) {
    if (n %% floor(nIterations/10) == 0) {
      print(paste(floor(100*n/nIterations), " %"))
    }
    # Randomly flip sign of all variables of observations (rows)
    ToTest <- Permed0(DF)
    O <- nSigFunc(ToTest)
    nSigvec <- c(nSigvec, O[[1]])
    smallest_p_vec <- c(smallest_p_vec, min(O[[2]]))
  }
  # Compare observed number of nomially significant tests with H0 distribution
  OObs <- nSigFunc(DF)
  p <- length(which(nSigvec >= OObs[[1]])) / length(nSigvec)
  # Get pPBFWC_, the minimal p-value over the group of tests with 5% familywise error rate 
  pPBFWC_ <- quantile(smallest_p_vec, 0.05)
  # Return:
  #   p-value of number of nominally significant results
  #   vector of nominally significant of tests
  #   vector of significant tests with Bonferroni-correction for expected number of nominally significant tests
  #   Permutation-based FWE-controlling p-value
  #   vector of significant tests using pPBFWC_
  nNomDiv <- mean(nSigvec)
  if (nNomDiv == 0) {
    nNomDiv <- 1
  }
  Output <- list(p, ((p < 0.05) + 0) * (OObs[[2]] < 0.05) + 0, ((p < 0.05) + 0) * (OObs[[2]] < (0.05/nNomDiv)) + 0, pPBFWC_, (OObs[[2]] < pPBFWC_) + 0)
}

# Tests with simulated data
nIterations <- 1000
nVariables <- 200
nObservations <- 100
common_signal_strength <- 0.5
pvec <- c()
sigM <- matrix(NA, nrow = 0, ncol = nVariables)
corrSigM <- matrix(NA, nrow = 0, ncol = nVariables)
pPBFWC_vec <- c()
pPBFWC_SigM <- matrix(NA, nrow = 0, ncol = nVariables)
for (n in 1:nIterations) {
  print(paste("Iteration ", n, " of ", nIterations))
  # Generate a data frame
  M <- matrix(NA, nrow = nObservations, ncol = 0)
  varnames <- c()
  for (iColumn in 1:nVariables) {
    v <- rnorm(nObservations)
    M <- cbind(M, v)
    varnames <- c(varnames, paste("var", iColumn))
  }
  common_signal <- matrix(rnorm(prod(dim(M))), nrow <- nObservations)
  M <- (1 - common_signal_strength) * M + common_signal_strength * common_signal
  DF <- data.frame(M)
  names(DF) <- varnames
  # Determine which variables are labelled significant
  O <- NumNomS(DF)
  # Store results
  pvec <- c(pvec, O[[1]])
  sigM <- rbind(sigM, O[[2]])
  corrSigM <- rbind(corrSigM, O[[3]])
  pPBFWC_vec <- c(pPBFWC_vec, O[[4]])
  pPBFWC_SigRow <- O[[5]]
  pPBFWC_SigRow[FALSE] <- 0
  pPBFWC_SigM <- rbind(pPBFWC_SigM, pPBFWC_SigRow)
}
report0 <- function() {
  # Report simulation results
  # Proportion of samples with significant NumNomS test
  print(paste("Proportion of samples with a significant NumNomS test: ", length(which(pvec < .05)) / length(pvec)))
  # Proportion of samples with any significant pPBFWC_-level tests
  y <- apply(pPBFWC_SigM, 1, function(r){max(r)})
  print(paste("Proportion of samples with any significant pPBFWC_ tests: ", length(which(y > 0)) / length(y)))
  # Expected number of significant tests, uncorrected
  y <- apply(sigM, 1, function(r){sum(r)})
  print(paste("Expected number of significant tests, uncorrected: ", mean(y)))
  # Expected number of significant tests, corrected
  y <- apply(corrSigM, 1, function(r){sum(r)})
  print(paste("Expected number of significant tests, corrected: ", mean(y)))
  # Critical p-value improvement with pPBFWC_ versus Bonferroni
  print(paste("Critical p-value improvement with pPBFWC_ versus Bonferroni: ", mean(pPBFWC_vec/(0.05/nVariables))))
}
report0()
