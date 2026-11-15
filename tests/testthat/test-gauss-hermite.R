test_that("Gauss-Hermite nodes integrate N(0, 1) moments exactly up to degree 2n - 1", {
  dfact <- function(k) if (k <= 0) 1 else prod(seq(2 * k - 1, 1, by = -2))
  for (n in c(1, 5, 20, 30)) {
    nd <- irtpwr:::gh_nodes(n)
    expect_equal(sum(nd$ww), 1, tolerance = 1e-13)
    expect_true(all(diff(nd$xx) > 0))
    for (k in 0:(n - 1)) {
      expect_equal(sum(nd$ww * nd$xx^(2 * k)), dfact(k), tolerance = 1e-12)
      expect_equal(sum(nd$ww * nd$xx^(2 * k + 1)), 0, tolerance = 1e-10 * dfact(k + 1))
    }
  }
  expect_equal(irtpwr:::gauss_hermite(function(x) exp(x), order = 30), exp(0.5), tolerance = 1e-13)
})
