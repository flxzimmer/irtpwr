# B4: in calc.pvalues (and score_obs/grad_obs in general), the DIF Score and
# Gradient statistics must use the group variable, not "first half / second
# half" of the rows.

test_that("DIF Score/Gradient p values do not depend on row order or group sizes", {
  set.seed(1)
  g1 <- g2 <- list(a = rlnorm(5, sdlog = .2), d = rnorm(5))
  g2$a[1] <- g2$a[1]^2
  g2$d[1] <- g2$d[1] + .5
  hd <- setup.hypothesis("DIF2PL", list(g1, g2))
  set.seed(2)
  d <- ns("setup.data")(hd, n = 2000)
  # drop some rows of group A so that the groups are unbalanced
  keep <- c(301:2000)
  dat <- as.data.frame(d$data)[keep, ]
  grp <- d$group[keep]
  set.seed(3)
  p_sorted <- ns("calc.pvalues")(dat, hd, group = grp)
  o <- sample(nrow(dat))
  set.seed(3)
  p_shuffled <- ns("calc.pvalues")(dat[o, ], hd, group = grp[o])
  expect_equal(p_shuffled, p_sorted, tolerance = 1e-4)
  expect_true(all(p_sorted[c("Score", "Gradient")] < 0.01))
})
