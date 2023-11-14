#' Carry out Hypotheses Tests
#'
#' @param dat data frame containing the data.
#' @param hyp Hypothesis Object created by the setup.hypothesis function.
#' @param group Optional argument for when there are groups in the hypothesis. Numeric vector of group variable or, alternatively, name of the group variable in the data frame dat.
#'
#' @return
#' @export
#'
#' @examples
#'
#' library(mirt)
#'
#' ## 1PL vs 2PL
#' dat <- expand.table(LSAT7)
#' mirtfit <- mirt(dat,1,verbose = FALSE)
#
#' # Setting up the alternative hypothesis
#' hyp <- setup.hypothesis(type = '1PLvs2PL', altpars = mirtfit)
#'
#' # Testing the hypothesis
#' test <- calc.pvalues(data = dat, hyp = hyp, group = group)
#' test
#'
#'
#' ## DIF in a 2PL model
#'
#' dat <- expand.table(LSAT7)
#' group <- as.factor(c(rep(0, times = 250), rep(1, times = 750)))
#'
#' # Setting up the alternative hypothesis
#' group1 <- group2 <- list(a = rlnorm(5, sdlog = 0.2), d = rnorm(5))
#' group2$a[1] <- (group2$a[1])^2
#' group2$d[1] <- group2$d[1] + 0.5
#' altpars <- list(group1, group2)
#' hyp <- setup.hypothesis(type = "DIF2PL", altpars = altpars)
#'
#' # Testing the hypothesis
#' test <- calc.pvalues(dat = dat, hyp = hyp, group = group)
#' test
#'
calc.pvalues <- function(dat = dat, hyp = hyp, group = NULL){

  datx = list()
  if(length(group)==1) { # preparing data and group variable if given as character
    datx$data <- dat[, !names(dat) %in% group]
    datx$group = dat[group] |> unlist()
  } else {
    datx$data = dat
    datx$group = group
  }

  df <-  nrow(hyp$resmod$Amat)
  fitted <- irtpwr:::mml.fit(data = datx,hyp = hyp)
  stats <- irtpwr:::stat_obs(fitted)
  lpvals <- pchisq(stats,df,lower.tail=F)
}


