# simtools ----------------------------------------------------------------

#' Create artificial dataset from alternative hypothesis
#'
#' @param hyp Hypothesis Object created from the setup_hypothesis function
#' @param n integer, number of persons
#' @param dist.fun function to generate the person parameters
#'
#' @return
#' @export
#'
#' @examples
#'
#' altpars <- list(
#' a = rlnorm(5,sdlog = .4),
#' d = rnorm(5))
#' hyp <- setup_hypothesis(type = "1PLvs2PL", altpars = altpars)
#' data <- setup.data(hyp=hyp,n=500)
#'
setup.data = function(hyp, n,dist.fun=rnorm) {

  distfun = function(x) {dist.fun(x) %>% matrix(.,ncol=1)}

  if (isTRUE(hyp$unresmod$multigroup)) { # Multigroup Model

    df = lapply(hyp$unresmod$parsets,function(pars) mirt::simdata(a = pars$a,d = pars$d,Theta = distfun(n/2),itemtype = hyp$unresmod$itemtype)) %>% do.call(rbind,.)
    group=rep(c("A","B"),each=n/2)
    re=list(data = df,group=group)

  } else { # SingleGroup Model

    pars = hyp$unresmod$parsets
    is.multi = "a2" %in% colnames(pars$a)
    if(is.multi) {
      distfun = function(x) {dist.fun(2*x) %>% matrix(.,ncol=2)}
    }

    if (is.null(pars$g)) {pars$g = 0}
# browser()
    df = mirt::simdata(a = pars$a,
                       d = pars$d,
                       guess=pars$g,
                       Theta = distfun(n),
                       itemtype = hyp$unresmod$itemtype)
    re=list(data = df)
  }

    return(re)
}

