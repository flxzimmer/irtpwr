# B8: calc.pvalues returns its p values visibly; the documented call works
test_that("calc.pvalues example call works and returns p values visibly", {
  dat <- mirt::expand.table(mirt::LSAT7)
  mirtfit <- mirt::mirt(dat, 1, verbose = FALSE)
  hyp <- setup.hypothesis(type = "1PLvs2PL", altpars = mirtfit)
  set.seed(1)
  res <- withVisible(ns("calc.pvalues")(dat = dat, hyp = hyp))
  expect_true(res$visible)
  expect_named(res$value, c("Wald", "LR", "Score", "Gradient"))
  expect_true(all(res$value > 0 & res$value < 1))
})
