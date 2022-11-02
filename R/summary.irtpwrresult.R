#' Summary of the Power Analysis
#'
#' Output the resulting power or sample size for each statistic
#'
#' @param object Object of class irtpwrresult as created by the irtpwr function
#' @param ... additional arguments to be passed
#'
#' @return An object of class summary.irtpwrresult
#' @export
#'
#' @examples
#'
#' library(mirt)
#' dat <- expand.table(LSAT7)
#' mirtfit <- mirt(dat,1,verbose = FALSE)
#' hyp <- setup.hypothesis(type = "1PLvs2PL", altpars = mirtfit)
#' res <- irtpwr(hyp=hyp,alpha=.05)
#' summary(res)
#'
summary.irtpwrresult <- function(object, ..., power = NULL, N = NULL, alpha = NULL) {

  ds <- object

  if (!is.null(power)) {ds$power = power; ds$N = NULL}
  if (!is.null(N)) {ds$N = N; ds$power = NULL}
  if (!is.null(alpha)) ds$alpha = alpha


  # both power and N specified
  if(!is.null(ds$power) & !is.null(ds$N)) stop("only one out of N and power can be specified as something other than NULL")


  if (is.null(ds$N)) {
    if (is.null(ds$alpha)|is.null(ds$power)) stop("At least one argument is missing. Check if you specified two out of N, power, alpha")
    ds$result = sapply(ds$ncps, function(ncp) calc.N(ds$hyp,ncp,power = ds$power,alpha=ds$alpha))
  }

  if (is.null(ds$power)) {
    if (is.null(ds$alpha)|is.null(ds$N)) stop("At least one argument is missing. Check if you specified two out of N, power, alpha")
    ds$result = sapply(ds$ncps, function(ncp) calc.power(ds$hyp,ncp,ssize = ds$N,alpha=ds$alpha))
  }

  names(ds$result) = names(ds$ncps)
  class(ds) <- "summary.irtpwrresult"

  return(ds)
}




#' Print Summary of the search result
#'
#' @param x Object of class irtpwrresult as created by the irtpwr function
#' @param ... additional arguments to be passed
#'
#' @return An object of class summary.irtpwrresult
#' @export
#'
#' @examples
#'
#' library(mirt)
#' dat <- expand.table(LSAT7)
#' mirtfit <- mirt(dat,1,verbose = FALSE)
#' hyp <- setup.hypothesis(type = "1PLvs2PL", altpars = mirtfit)
#' res <- irtpwr(hyp=hyp,alpha=.05)
#' summary(res)
#'
print.summary.irtpwrresult <- function(x, ...) {

  method = x$method

  method <- switch(method, analytical = "Analytical",
                      sampling = "Sampling-Based")

  if (is.null(x$power)) sentence = paste0("Power for N = ",x$N," (alpha = ",x$alpha,"):")
  if (is.null(x$N)) sentence = paste0("Sample sizes for power = ",x$power," (alpha = ",x$alpha,"):")

  # Introductory Sentence: Power for N = X
  cat("\n",sentence,"\n")

  # Power or Sample Size as a table
  if (is.null(x$N)) mat = data.frame(Statistic = names(x$ncps),N=x$result)
  if (is.null(x$power)) mat = data.frame(Statistic = names(x$ncps),Power=round(x$result,4),row.names=NULL)


  cat("\n")

  print(mat,row.names=F)

  cat("\nMethod: ", paste(method, sep = "\n",
                             collapse = "\n"), sep = "")

  cat("\n")
  cat("\n")

  invisible(x)
}
