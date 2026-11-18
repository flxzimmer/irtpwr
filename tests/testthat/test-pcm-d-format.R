# PCMvsGPCM: d with or without the leading zero column is detected robustly

setup_gpcm <- function(d) {
  setup.hypothesis("PCMvsGPCM", list(a = c(.8, 1, 1.3), d = d))
}

test_that("d with and without the zero column give the same hypothesis", {
  thr <- cbind(c(.5, .2, -.3), c(.1, -.4, -.8))
  ref <- setup_gpcm(cbind(0, thr))
  expect_identical(ref$unresmod$parsets$d, cbind(0, thr))
  expect_identical(setup_gpcm(thr)$unresmod, ref$unresmod)
  expect_identical(setup_gpcm(thr)$resmod$Amat, ref$resmod$Amat)
})

test_that("an all-zero first threshold column is not mistaken for the zero column", {
  thr <- cbind(c(0, 0, 0), c(.1, -.4, -.8))
  h <- setup_gpcm(thr)
  expect_identical(h$unresmod$parsets$d, cbind(0, thr))
  expect_identical(ncol(h$resmod$Amat), 3L * 3L)
  # with the zero column it stays as given
  expect_identical(setup_gpcm(cbind(0, thr))$unresmod$parsets$d, cbind(0, thr))
})

test_that("mirt-style column names are handled", {
  thr <- cbind(d1 = c(0, 0, 0), d2 = c(.1, -.4, -.8))
  expect_identical(unname(setup_gpcm(thr)$unresmod$parsets$d),
                   unname(cbind(0, thr)))
  named0 <- cbind(d0 = 0, d1 = c(.5, .2, -.3), d2 = c(.1, -.4, -.8))
  expect_identical(setup_gpcm(named0)$unresmod$parsets$d, named0)
})

test_that("d for other than 3 categories stops early with a clear error", {
  # vector d: one threshold per item, 2 categories
  expect_error(setup_gpcm(c(.5, .2, -.3)),
               "'d' is a vector.*exactly 3 response categories")
  # one threshold column, with or without the zero column: 2 categories
  expect_error(setup_gpcm(cbind(d0 = 0, d1 = c(.5, .2, -.3))),
               "2 response categories.*exactly 3")
  expect_error(setup_gpcm(cbind(c(.5, .2, -.3))),
               "2 response categories.*exactly 3")
  # three threshold columns: 4 categories
  thr3 <- cbind(c(.5, .2, -.3), c(.1, -.4, -.8), c(-.2, -.6, -1))
  expect_error(setup_gpcm(thr3), "4 response categories.*exactly 3")
  expect_error(setup_gpcm(cbind(0, thr3)), "4 response categories.*exactly 3")
})
