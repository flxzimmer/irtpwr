# B6: analytical 3PL is not supported -> informative error
test_that("analytical 3PL gives an informative error", {
  altpars <- list(a = c(1.1, seq(.5, 1.4, length.out = 4)),
                  d = c(.1, seq(1.3, -1.3, length.out = 4)),
                  g = c(.22, rep(.2, 4)))
  hyp <- setup.hypothesis("3PL_basic", altpars)
  expect_error(irtpwr(hyp), "not supported for 3PL")
  expect_error(irtpwr(hyp, stat = "Wald"), "not supported for 3PL")
})

# B7: unknown hypothesis type -> clear error
test_that("unknown hypothesis type gives an informative error", {
  altpars <- list(a = c(1, 1.2, .8), d = c(0, .3, -.2))
  expect_error(setup.hypothesis("nonsense", altpars),
               "Unknown hypothesis type 'nonsense'")
})

# unknown method / empty stat -> clear errors (before: "object 'ncps' not
# found" and a silent empty result)
test_that("irtpwr checks 'method' and 'stat'", {
  hyp <- hyp_1pl3()
  expect_error(irtpwr(hyp, method = "foo"),
               "should be one of .*analytical.*sampling")
  expect_error(irtpwr(hyp, stat = character(0)), "'stat' is empty")
  expect_error(irtpwr(hyp, stat = NULL), "'stat' is empty")
  # the resolved method is stored, partial matching as in match.arg
  res <- irtpwr(hyp, stat = "Wald", method = "anal")
  expect_identical(res$method, "analytical")
  expect_identical(names(res$ncps), "Wald")
})

# PCMvsGPCM supports exactly 3 categories: errors come from
# setup.hypothesis(), not later from irtpwr()
test_that("PCMvsGPCM with other than 3 categories stops in setup.hypothesis", {
  a <- c(.8, 1, 1.3)
  expect_error(setup.hypothesis("PCMvsGPCM", list(a = a, d = c(.5, .2, -.3))),
               "exactly 3 response categories")
  d4 <- cbind(0, c(.5, .2, -.3), c(.1, -.4, -.8), c(-.2, -.6, -1))
  expect_error(setup.hypothesis("PCMvsGPCM", list(a = a, d = d4)),
               "4 response categories")
})

# power and alpha must be single proportions in (0, 1); before, power = 80
# or alpha = 5 were accepted and summary() failed later with "f() values at
# end points not of opposite sign"
test_that("irtpwr and summary check 'power' and 'alpha'", {
  hyp <- hyp_1pl3()
  for (p in list(80, 0, 1, -0.2, c(.8, .9), NA_real_, "0.8")) {
    expect_error(irtpwr(hyp, stat = "Wald", power = p, alpha = .05),
                 "'power' must be a single number strictly between 0 and 1")
  }
  for (a in list(5, 0, 1, c(.05, .01), NA_real_)) {
    expect_error(irtpwr(hyp, stat = "Wald", power = .8, alpha = a),
                 "'alpha' must be a single number strictly between 0 and 1")
  }
  res <- irtpwr(hyp, stat = "Wald", power = .8, alpha = .05)
  expect_error(summary(res, power = 80), "'power' must be a single number")
  expect_error(summary(res, alpha = 5), "'alpha' must be a single number")
  # valid inputs unchanged
  expect_identical(summary(res)$result,
                   summary(irtpwr(hyp, stat = "Wald"), power = .8,
                           alpha = .05)$result)
})

# N must be a single positive number; before, N = c(300, 600) gave a summary
# table with wrongly recycled row labels
test_that("irtpwr and summary check 'N'", {
  hyp <- hyp_1pl3()
  for (n in list(c(300, 600), 0, -100, NA_real_, Inf, "600")) {
    expect_error(irtpwr(hyp, stat = "Wald", N = n, alpha = .05),
                 "'N' must be a single positive number")
  }
  res <- irtpwr(hyp, stat = "Wald", N = 600, alpha = .05)
  expect_error(summary(res, N = c(300, 600)),
               "'N' must be a single positive number")
  expect_error(summary(res, N = 0), "'N' must be a single positive number")
  expect_length(summary(res, N = 300)$result, 1)
})

# analytical multidimensional (multi_basic) is not supported -> informative
# error (before: "non-conformable arguments" / "tab does not match idx")
test_that("analytical multi_basic gives an informative error", {
  altpars <- list(a = cbind(a1 = c(1, 1, 1), a2 = c(.5, .5, 0)),
                  d = c(0, 0, 0))
  hyp <- setup.hypothesis("multi_basic", altpars)
  expect_error(irtpwr(hyp), "not supported for multidimensional models")
  expect_error(irtpwr(hyp, stat = "Wald"),
               "not supported for multidimensional models")
})
