# DIF2PL with exactly one differing parameter (one-row Amat) must work
test_that("DIF2PL with a single differing parameter works", {
  g1 <- g2 <- list(a = c(1, 1.2, .8), d = c(0, .3, -.2))
  g2$d[1] <- .5
  hd <- setup.hypothesis("DIF2PL", list(g1, g2))
  expect_identical(dim(hd$resmod$Amat), c(1L, 7L))
  r <- irtpwr(hd)$ncps
  expect_named(r, c("Wald", "LR", "Score", "Gradient"))
  expect_true(all(r > 0))
  # one parameter: all four ncps are close to each other
  expect_lt(diff(range(r)) / mean(r), 0.1)
})
