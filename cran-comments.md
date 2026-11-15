## Submission

This is an update of irtpwr from 1.0.3 to 1.1.0. See NEWS.md. In short:

* Bug fixes, some of which change results. Most notably, the analytical
  LR, score and gradient sample sizes for `type = "basic"` were too small.
* The analytical method is 11 to 105 times faster, the sampling-based
  method about 2.6 times, with identical results.
* A small C routine (registered with `R_registerRoutines()`).
* The dependencies on spatstat.random, Deriv and digest are dropped.
* `calc.pvalues()` is exported.
* `Depends: R (>= 4.1.0)`, which fixes the NOTE on all CRAN check flavours
  for 1.0.3 (use of the pipe `|>`).

## Test environments

* local: Ubuntu 26.04 LTS, R 4.5.2; also with `_R_CHECK_DEPENDS_ONLY_=true`,
  and with AddressSanitizer and UndefinedBehaviorSanitizer for the C code
* GitHub Actions: macOS (R release), Windows (R release), Ubuntu (R devel,
  release, oldrel-1)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies on CRAN.
