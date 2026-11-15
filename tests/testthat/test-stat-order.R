# B1: result labels must follow the values, whatever the order of `stat`

test_that("analytical ncps are labelled correctly for any order of stat", {
  hyp <- hyp_1pl3()
  ref <- irtpwr(hyp, stat = c("Wald", "LR"))$ncps
  rev <- irtpwr(hyp, stat = c("LR", "Wald"))$ncps
  expect_named(ref, c("Wald", "LR"))
  expect_identical(rev[c("Wald", "LR")], ref)
  expect_equal(unname(ref["Wald"]), ns("ncp.wald")(hyp))
  sg <- irtpwr(hyp, stat = c("Gradient", "Score"))$ncps
  expect_identical(sg[["Score"]], ns("ncp.score")(hyp)$ncp)
})

test_that("sampling ncps are labelled correctly for any order of stat", {
  hyp <- setup.hypothesis("1PLvs2PL", list(a = c(0.5, 1.5, 2.5),
                                           d = c(-0.5, 0.2, 0.9)))
  set.seed(1)
  ref <- irtpwr(hyp, stat = c("Wald", "LR"), method = "sampling",
                sampling.npers = 1000, approx.npers = 1000)$ncps
  set.seed(1)
  rev <- irtpwr(hyp, stat = c("LR", "Wald"), method = "sampling",
                sampling.npers = 1000, approx.npers = 1000)$ncps
  expect_named(ref, c("Wald", "LR"))
  expect_true(ref[["Wald"]] != ref[["LR"]])
  expect_identical(rev[c("Wald", "LR")], ref)
})

test_that("unknown statistics give an informative error", {
  expect_error(irtpwr(hyp_1pl3(), stat = "wald"), "Unknown statistic")
})
