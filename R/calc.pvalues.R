#' Carry out Hypothesis Tests on Observed Data
#'
#' Fits the unrestricted and the restricted model of a hypothesis to observed
#' data and carries out the Wald, likelihood ratio (LR), score, and gradient
#' tests of the null hypothesis. This is the step after the power analysis:
#' once the data are collected, the same hypothesis object can be used to test
#' the hypothesis.
#'
#' The models are fitted with marginal maximum likelihood (mirt). The
#' statistics are asymptotically chi-square distributed under the null
#' hypothesis, with degrees of freedom equal to the number of restrictions
#' (rows of \code{hyp$resmod$Amat}). For hypotheses with two groups (such as
#' \code{DIF2PL}), the group variable must be given. The Wald and score statistics
#' use an information matrix computed for two groups of equal size, so they
#' are best suited to roughly balanced groups.
#'
#' @param dat data frame containing the item responses (and, optionally, the
#'   group variable).
#' @param hyp Hypothesis object created by the \code{\link{setup.hypothesis}} function.
#' @param group Optional argument for hypotheses with groups (e.g. \code{DIF2PL}).
#'   Either a vector with the group of each row of \code{dat}, or the name of the
#'   group variable in the data frame \code{dat}.
#'
#' @return Named numeric vector with the p values of the Wald, LR, Score, and
#'   Gradient tests.
#' @export
#'
#' @examples
#' library(mirt)
#'
#' ## Rasch (1PL) vs 2PL model on the LSAT7 data
#' dat <- expand.table(LSAT7)
#' mirtfit <- mirt(dat, 1, verbose = FALSE)
#' hyp <- setup.hypothesis(type = "1PLvs2PL", altpars = mirtfit)
#' calc.pvalues(dat = dat, hyp = hyp)
#'
#' ## DIF in the first item of a 2PL model, with a group variable
#' set.seed(1)
#' a <- rlnorm(5, sdlog = 0.2)
#' d <- rnorm(5)
#' d2 <- d
#' d2[1] <- d2[1] + 0.5
#' difdat <- rbind(simdata(a, d, 500, itemtype = "2PL"),
#'                 simdata(a, d2, 500, itemtype = "2PL"))
#' difdat <- data.frame(difdat, group = rep(c("A", "B"), each = 500))
#' hyp <- setup.hypothesis(type = "DIF2PL",
#'   altpars = list(list(a = a, d = d), list(a = a, d = d2)))
#' calc.pvalues(dat = difdat, hyp = hyp, group = "group")
calc.pvalues <- function(dat, hyp, group = NULL){

  datx = list()
  if(length(group)==1) { # preparing data and group variable if given as character
    datx$data <- dat[, !names(dat) %in% group]
    datx$group = dat[group] |> unlist()
  } else {
    datx$data = dat
    datx$group = group
  }

  df <-  nrow(hyp$resmod$Amat)
  fitted <- mml.fit(data = datx,hyp = hyp)
  stats <- stat_obs(fitted)
  lpvals <- stats::pchisq(stats,df,lower.tail=F)
  lpvals
}


