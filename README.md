
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pwrml

<!-- badges: start -->
<!-- badges: end -->

This package offers implementation of analytical and sampling-based
power analyses for the Wald, LR, score, and gradient tests in the
frameworks of linear hypotheses and marginal maximum likelihood
estimation.

The vignette “demo” gives an introduction of basic features,
“adding\_hypotheses” contains a tutorial on setting up custom
hypotheses. Further details can be found in our paper “Power analysis
for the Wald, LR, score, and gradient tests in a marginal maximum
likelihood framework: Applications in IRT”.

## Installation

``` r
# library(devtools)
# install_github('anonymized')
```

## To-Dos:

-   speed up calculation of score statistic using internals from the
    mirt functions, e.g. attributes(fitted*r**e**s*)Internals
