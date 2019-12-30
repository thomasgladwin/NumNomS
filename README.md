# NumNomS

Based on Gladwin, Derks, et al. (2012) and Gladwin & Vink (2018).

This is a "no smoke without fire" multiple testing correction that trades off precision for statistical power.

The NumNomS set-wise test considers the Number of Nominally Significant tests.

A null hypothesis distribution of this number is determined by permutation testing, such that correlations between variables being tested can be exploited to improve power.

The software also provides a test-wise permutation-based familywise-corrected critical p-value, pPBFW, that should be at least as powerful as Bonferroni correction.

Additionally, while this is not principled at this point, the software reports which tests are significant when using Bonferroni correction over the expected number of nominally significant tests. This is intended to be used only after protection by NumNomS significance, to reduce the number of false test-wise positives after a false NumNomS positive to a similar level as for a single test.

Note that the number of iterations can be adjusted to address any concerns with the use of permutation tests to a trivial level, relative to the inherent uncertainity of sampling.

NumNomS.R contains the core NumNomS function, adjustable functions to change the performed tests and permutation, and simulation code.
