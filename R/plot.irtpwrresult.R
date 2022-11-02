#' Plot Power Curves
#'
#' Plot the power curves using the resulting object of the irtpwr function.
#'
#' @param x Object of class irtpwrresult as created by the irtpwr function.
#' @param bounds integer vector. the first entry is the lower bound of the x-axis in the plot (sample size). The second entry is the upper bound. By default, these values are chosen to cover a power range of .5 to .95.
#' @param ... additional arguments to be passed.
#'
#' @return A ggplot object
#' @export
#'
#' @examples #Load a simulation function
#'
#' library(mirt)
#' dat <- expand.table(LSAT7)
#' mirtfit <- mirt(dat,1,verbose = FALSE)
#' hyp <- setup.hypothesis(type = "1PLvs2PL", altpars = mirtfit)
#' res <- irtpwr(hyp=hyp,alpha=.05,power =.8)
#' plot(res)
#'
plot.irtpwrresult <- function(x, bounds  = NULL, ...) {

  if(is.null(bounds)) {
    #specify some default bounds
    bounds = c(
      min(calc.N(x$hyp,x$ncps,power = .5,alpha=x$alpha)),
      max(calc.N(x$hyp,x$ncps,power = .95,alpha=x$alpha))
    )
  }

  ns = seq(bounds[1],bounds[2])
  a =  lapply(x$ncps, function(ncp) calc.power(x$hyp,ncp,ssize = ns,alpha=x$alpha))
  a = do.call(rbind,a)
  dat = data.frame(power = as.numeric(a), stat = rep(names(x$ncps),times = length(ns)),N=rep(ns,each=length(x$ncps)))

  pl = ggplot2::ggplot() + ggplot2::geom_line(ggplot2::aes(x=dat$N,y=dat$power,color=dat$stat)) + ggplot2::xlab("Sample Size") + ggplot2::ylab("Power") + ggplot2::labs(color = "Statistic")



  print(pl)
  invisible(pl)

}
