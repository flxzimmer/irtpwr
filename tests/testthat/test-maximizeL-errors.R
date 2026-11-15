# B10: an allocation error in maximizeL gives the hint once and is then
# rethrown; maximizeL is not called a second time.

test_that("allocation errors in maximizeL are rethrown after the hint", {
  hyp <- hyp_1pl3()
  calls <- 0
  hyp$type$maximizeL <- function(hyp) {
    calls <<- calls + 1
    stop("cannot allocate vector of size 12.3 Gb")
  }
  expect_message(
    expect_error(irtpwr(hyp, stat = "LR"), "cannot allocate vector"),
    "too many items")
  expect_identical(calls, 1)
})

test_that("other errors in maximizeL propagate unchanged", {
  hyp <- hyp_1pl3()
  hyp$type$maximizeL <- function(hyp) stop("some other problem")
  expect_no_message(expect_error(irtpwr(hyp, stat = "LR"), "some other problem"))
})
