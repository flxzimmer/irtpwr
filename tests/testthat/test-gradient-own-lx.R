# B2: in the sampling approach, the Gradient statistic computed on its own
# must equal the one computed together with the Score statistic (which
# hands its lx over).

grad_alone_vs_with_score <- function(hyp, n) {
  data <- ns("setup.data")(hyp, n = n)
  fitted <- ns("mml.fit")(hyp, data = data, infmat.unres = "ApproxFisher",
                          infmat.res = "ApproxFisher", approx.npers = n,
                          SE.type = "Fisher")
  c(alone = ns("stat_obs")(fitted, "Gradient")[["Gradient"]],
    with_score = ns("stat_obs")(fitted, c("Score", "Gradient"))[["Gradient"]])
}

test_that("3PL: Gradient alone equals Gradient with Score", {
  skip_on_cran()
  altpars <- list(a = c(1.1, seq(.5, 1.4, length.out = 4)),
                  d = c(.1, seq(1.3, -1.3, length.out = 4)),
                  g = c(.22, rep(.2, 4)))
  hyp <- setup.hypothesis("3PL_basic", altpars)
  set.seed(3)
  r <- grad_alone_vs_with_score(hyp, 2000)
  expect_identical(r[["alone"]], r[["with_score"]])
})

test_that("multidimensional: Gradient alone works and equals Gradient with Score", {
  skip_on_cran()
  altpars <- list(a = cbind(a1 = c(1.2, .8, 1, 1.5, .9),
                            a2 = c(.6, .9, .4, .7, 0)),
                  d = c(.3, .3, -.4, .8, -.2))
  altpars$d[2] <- .8
  hyp <- setup.hypothesis("multi_basic", altpars)
  set.seed(4)
  r <- grad_alone_vs_with_score(hyp, 2000)
  expect_identical(r[["alone"]], r[["with_score"]])
})
