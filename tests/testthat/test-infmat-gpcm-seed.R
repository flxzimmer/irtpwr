# B3: the analytical (Fisher) information for gpcm must not depend on the
# RNG seed, i.e. the thresholds must be fixed at their given values.
test_that("Fisher infmat for gpcm does not depend on the seed", {
  gp <- list(a = c(.8, 1, 1.3, 1.6),
             d = cbind(0, c(.5, .2, -.3, .9), c(.1, -.4, -.8, .6)))
  set.seed(1)
  i1 <- ns("infmat")(gp, method = "Fisher", itemtype = "gpcm")
  set.seed(2)
  i2 <- ns("infmat")(gp, method = "Fisher", itemtype = "gpcm")
  expect_equal(i1, i2, tolerance = 1e-8)
  hg <- setup.hypothesis("PCMvsGPCM", gp)
  set.seed(1); w1 <- ns("ncp.wald")(hg)
  set.seed(2); w2 <- ns("ncp.wald")(hg)
  expect_equal(w1, w2, tolerance = 1e-8)
})
