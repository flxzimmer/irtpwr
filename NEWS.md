# irtpwr 1.1.0

This release fixes several bugs, some of which change results, makes the
analytical method much faster, and drops three dependencies. Results that
are not listed under "Changes that affect results" are unchanged.

## New features

* `calc.pvalues()` is now exported: it runs the Wald, LR, Score and
  Gradient tests on observed data for a hypothesis set up with
  `setup.hypothesis()`. It was announced in 1.0.3 but could not be called.

## Changes that affect results

* `type = "basic"`, analytical method: the restricted maximum likelihood
  under H0 (d1 = 0) now optimises all free item parameters, not only the
  first slope. The LR, score, and gradient noncentrality parameters were
  too large before, so **the required sample sizes were too small: expect
  about 3% to 14% larger N** (measured for 18 parameter sets with 3 to 7
  items; larger with more items). For example, the LR test needs 1014
  instead of 951 persons for power .8 in the 5-item example of the
  "adding_hypotheses" vignette, which now sets a seed. Wald results are
  unchanged. The template in the "adding_hypotheses" vignette was fixed in
  the same way; custom hypotheses copied from the old template have the
  same problem and should be updated.
* `type = "PCMvsGPCM"`, analytical method: the Fisher information now fixes
  the GPCM thresholds at their given values instead of estimating them from
  a random sample. Wald and score results therefore no longer depend on the
  random seed. They change by a few percent in either direction, depending
  on the item parameters and on the seed the previous result used: in 220
  combinations of parameter sets (3 to 6 items) and seeds, N changed by
  -6% to +9% for Wald and by -7% to +7% for score; 90% of the changes were
  between -3% and +5%, and half were smaller than 1.5%. Examples with 4
  items and seed 1: N 842 -> 888 for Wald (+5.5%) and 469 -> 476 for
  score (+1.5%) for one parameter set, 437 -> 438 for Wald and 311 -> 311
  for score for another. LR and gradient results are unchanged.
* Sampling-based method: the gradient statistic computed without the score
  statistic (`stat = "Gradient"` without `"Score"`) was wrong for 3PL
  models and failed for multidimensional models. It now equals the gradient
  statistic computed together with the score statistic.
* `calc.pvalues()` with DIF hypotheses: the score and gradient statistics
  now use the group variable. Before, they assumed that the data were
  sorted by group into two equal halves, and were wrong for unsorted or
  unbalanced data.
* Sampling-based method with DIF hypotheses and an odd `sampling.npers`:
  the generated data now have exactly `sampling.npers` rows (one row was
  missing before).
* `calc.time()` now uses the correct number of items for DIF and PCMvsGPCM
  hypotheses, and the number of categories for PCMvsGPCM. Estimates for
  DIF hypotheses were far too high before (32 times too high for 7 items).
* `calc.time()` is recalibrated for the faster analytical method. It now
  assumes a fixed overhead plus a part proportional to the number of
  response patterns, and times short runs twice (the first call in a
  session includes start-up costs). With the faster method, the previous
  extrapolation (run time times n.cat^(additional items)) gave estimates
  that were typically about 3 times too high (up to 8 times); now they are
  typically about 1.3 times too high. In our test runs, the actual time
  was usually between about half of the estimate and 25% more; for 14 or
  more items the estimate can be too low (actual time up to about 1.5
  times the estimate).

## Bug fixes

* When `stat` was not given in the order Wald, LR, Score, Gradient, the
  result names did not match the values (e.g. `stat = c("LR", "Wald")`
  swapped the labels). Unknown statistics now give an informative error.
* The analytical method is not available for 3PL models. It now stops with
  an informative error (use `method = "sampling"`). Before, it failed with
  an obscure error (LR, score, gradient) or returned a wrong Wald result.
* The analytical method is not available for multidimensional models
  (`type = "multi_basic"`). It now stops with an informative error (use
  `method = "sampling"`), as for 3PL models. Before, it failed with an
  obscure error (e.g. "non-conformable arguments").
* `setup.hypothesis()` gives an informative error for an unknown `type`.
* `setup.hypothesis("DIF2PL")` now works when the groups differ in only
  one parameter.
* `PCMvsGPCM`: the leading zero column of `d` is detected more reliably.
  `setup.hypothesis()` now stops with an informative error if `d` does
  not describe items with exactly 3 response categories (e.g. a vector
  `d` or more than two threshold columns), which is all `PCMvsGPCM`
  supports. Before, `irtpwr()` failed later with an obscure error
  ("subscript out of bounds").
* `irtpwr()` with an unknown `method` now gives an informative error
  (options are "analytical" and "sampling"); before, it failed with
  "object 'ncps' not found".
* `irtpwr()` with an empty `stat` now gives an informative error instead
  of silently returning an empty result.
* `irtpwr()` and `summary()` now check `power`, `alpha` and `N`: `power`
  and `alpha` must be single numbers strictly between 0 and 1 (proportions,
  e.g. `power = 0.8`), and `N` must be a single positive number. Before,
  `power = 80` or `alpha = 5` were accepted and `summary()` failed later
  with an obscure error ("f() values at end points not of opposite sign"),
  and a vector `N` (e.g. `N = c(300, 600)`) gave a summary table with
  wrong row labels.
* When the analytical method runs out of memory (too many items), the hint
  to use the sampling-based method is now followed directly by the original
  error, instead of repeating the whole computation first.
* `calc.pvalues()`: corrected example and documentation; the p values are
  now returned visibly.
* Vignettes: the printed `head(dat)` output in the "power_analysis"
  vignette now matches the data the code generates, and the example in
  the "adding_hypotheses" vignette sets a seed.
* Documentation: the README, the "demo" and "power_analysis" vignettes
  and `?setup.hypothesis` now list all six hypothesis types (`"basic"`,
  `"3PL_basic"` and `"multi_basic"` were missing; the "power_analysis"
  vignette said there were only three types).
  `?setup.hypothesis` now refers to the "adding_hypotheses" vignette for
  custom hypotheses instead of a section that did not exist, and no longer
  mentions the internal function `coef_short()`. `?irtpwr` and
  `?summary.irtpwrresult` state that `power` and `alpha` are proportions
  and that `power` and `N` must be single numbers.

## Performance

The speedups do not change the random number stream. Apart from the new
quadrature (see the next point), they give identical results, checked
against the previous version with `identical()`.

* Gauss-Hermite quadrature now uses its own implementation of the
  Golub-Welsch algorithm instead of `spatstat.random::gauss.hermite()`. At the
  orders irtpwr uses (20, 30 and 50), it is more accurate than the previous
  polynomial-root method. Results change only far beyond printed precision
  (relative differences below 1e-8); sample sizes and power values are unchanged.

* The analytical method is 11 to 105 times faster, with larger gains for
  more items and for GPCM models. The Gauss-Hermite nodes are computed once
  per order instead of on every call, the model derivatives are no longer
  derived symbolically at run time, and the integration over response
  patterns is vectorised with a small compiled kernel.
* The sampling-based method is about 2.6 times faster with the default
  settings, mainly because response patterns are counted without hashing
  every row, and only once per call.

## Dependencies

* Removed the dependencies on spatstat.random, Deriv, and digest.
* irtpwr now requires R >= 4.1.0.
* irtpwr now contains compiled C code, so installing it from source needs
  a C compiler (binary packages from CRAN do not).

# irtpwr 1.0.3

* Added a function to calculate statistics and p values for hypothesis tests (`calc.pvalues`)

# irtpwr 1.0.2

* Fixed missing words in a vignette

# irtpwr 1.0.1

* Added a new vignette "Power Analysis for the Wald, LR, Score, and Gradient Tests using irtpwr" authored by Jan Radek. It gives a more thorough power analysis walkthrough using the package. 

# irtpwr 1.0.0

* We thank Nebil Bel Hedi for extensive testing of the package, especially with respect to the accuracy of the power estimates.


