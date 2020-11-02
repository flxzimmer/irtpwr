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
#'
setup.data = function(hyp, n,dist.fun=rnorm) {
# for dist.type: replace with function input for a cleaner function (outsourcing the function creation)

    distfun = function(x) {dist.fun(x) %>% matrix(.,ncol=1)}

  if (isTRUE(hyp$unresmod$multigroup)) {
    df = lapply(hyp$unresmod$parsets,function(pars) simdata(a = pars$a,d = pars$d,Theta = distfun(n/2),itemtype = hyp$unresmod$itemtype)) %>% do.call(rbind,.)
    group=rep(c("A","B"),each=n/2)
    re=list(data = df,group=group)
  } else {
    if (is.null(hyp$unresmod$g)) {hyp$unresmod$g = 0}
    df = mirt::simdata(a = hyp$unresmod$a,d = hyp$unresmod$d,guess=hyp$unresmod$g,Theta = distfun(n),itemtype = hyp$unresmod$itemtype)
    re=list(data = df)
  }

    return(re)
}

