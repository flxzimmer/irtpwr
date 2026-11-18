# B5: the restricted ML for the 'basic' hypothesis (d1 = 0) must optimise
# all free parameters, so the expected score is zero for every parameter
# except the restricted d1.
test_that("basic maximizeL gives a stationary point of the restricted model", {
  altpars <- list(a = c(.8, 1.3, 1.6), d = c(.4, .2, -.6))
  hyp <- setup.hypothesis("basic", altpars)
  parsr <- hyp$type$maximizeL(hyp)
  expect_identical(parsr$d[1], 0)
  lx <- ns("ncp.score")(hyp, parsr = parsr)$lx
  # lx is ordered a1, d1, a2, d2, a3, d3; d1 is restricted
  expect_lt(max(abs(lx[-2])), 1e-4)
  expect_gt(abs(lx[2]), 1e-2)
})
