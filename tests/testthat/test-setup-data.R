# B11: multigroup setup.data generates exactly n rows, also for odd n
test_that("multigroup setup.data returns n rows and a matching group vector", {
  g1 <- g2 <- list(a = c(1, 1.2, .8), d = c(0, .3, -.2))
  g2$a[1] <- 1.5
  g2$d[1] <- .5
  hd <- setup.hypothesis("DIF2PL", list(g1, g2))
  for (n in c(5, 6, 101)) {
    d <- ns("setup.data")(hd, n = n)
    expect_identical(nrow(d$data), as.integer(n))
    expect_length(d$group, n)
    expect_identical(as.vector(table(d$group)),
                     as.integer(c(ceiling(n / 2), floor(n / 2))))
  }
})

test_that("even n gives the same data as before (equal halves)", {
  g1 <- g2 <- list(a = c(1, 1.2, .8), d = c(0, .3, -.2))
  g2$a[1] <- 1.5
  g2$d[1] <- .5
  hd <- setup.hypothesis("DIF2PL", list(g1, g2))
  set.seed(10)
  d <- ns("setup.data")(hd, n = 40)
  set.seed(10)
  ref <- rbind(mirt::simdata(a = g1$a, d = g1$d, Theta = matrix(rnorm(20)),
                             itemtype = "2PL"),
               mirt::simdata(a = g2$a, d = g2$d, Theta = matrix(rnorm(20)),
                             itemtype = "2PL"))
  expect_identical(d$data, ref)
  expect_identical(d$group, rep(c("A", "B"), each = 20))
})
