# B9: calc.time must extrapolate from the actual number of items (and
# categories) of the example hypothesis, also for DIF and GPCM hypotheses.
# Expected values follow the run time model in calc.time():
# 1.2 * time * (64 + n.cat^n.items) / (64 + n.cat^n.items.ex).

test_that("calc.time uses the right item and category counts", {
  # a fake clock: every (mocked) irtpwr() run takes exactly 0.2 seconds
  calls <- 0
  clock <- 0
  local_mocked_bindings(
    elapsed = function() clock,
    irtpwr = function(hyp, ...) {
      calls <<- calls + 1
      clock <<- clock + 0.2
      list(ncps = NA)
    })
  expected <- function(n.cat, n, n.ex) 1.2 * 0.2 * (64 + n.cat^n) / (64 + n.cat^n.ex)
  g1 <- g2 <- list(a = c(1, 1.2, .8, 1.1, .9), d = c(0, .3, -.2, .5, -.4))
  g2$a[1] <- 1.5
  g2$d[1] <- .5
  hd <- setup.hypothesis("DIF2PL", list(g1, g2))
  t_dif <- calc.time(hd, n.items = 10)
  e <- expected(2, 10, 5)                         # 2.72 (with n.items.ex = 0: 4.0)
  expect_equal(t_dif, e)
  expect_equal(calls, 2)  # short runs are timed twice
  hg <- setup.hypothesis("PCMvsGPCM", list(
    a = c(.8, 1, 1.3, 1.6),
    d = cbind(0, c(.5, .2, -.3, .9), c(.1, -.4, -.8, .6))))
  t_gpcm <- calc.time(hg, n.items = 7)
  e <- expected(3, 7, 4)                          # 3.73 (with base 2: 0.58)
  expect_equal(t_gpcm, e)
  h1 <- hyp_1pl3()
  t_1pl <- calc.time(h1, n.items = 3)
  expect_lt(t_1pl, 0.6)
})
