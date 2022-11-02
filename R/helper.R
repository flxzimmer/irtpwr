

logit <- function(p) { log(p/(1-p)) }
logitinv<- function(q) { 1/(1+exp(-q)) }


#' Calculate the computation time needed for the analytical method
#'
#' @param hyp Hypothesis object as created by the setup.hypothesis function
#' @param n.items Number of items
#'
#' @return Numeric
#' @export
#'
#' @examples
#'
#' library(mirt)
#' dat <- expand.table(LSAT7)
#' mirtfit <- mirt(dat,1,verbose = FALSE)
#' hyp <- setup.hypothesis(type = "1PLvs2PL", altpars = mirtfit)
#' calc.time(hyp,n.items=7)
#'
calc.time = function(hyp,n.items) {

  n.items.ex <- length(hyp$unresmod$parsets$d)

  ptm <- proc.time() #start timer
  ncps <- irtpwr(hyp=hyp)$ncps
  time <- proc.time() - ptm #stop timer

  time <- as.numeric(time[3])

  # time is t * 2^3
  re <- time *2^(n.items-n.items.ex)

  return(re)
}



findlambda = function(lx,A) {

  fn = function(k) {sum((lx + t(A)%*%k)^2)}

  k = optim(rep(0,nrow(A)),fn,method="BFGS")$par

  return(k)
}


coef_short = function(mirtfit,itemtype=NULL) {
  # extract coefs from mirtfit

  if (is.null(itemtype)) {
    itemtype <- mirt::extract.mirt(mirtfit, what = "itemtype")[1]
    }

  coefs <- mirt::coef(mirtfit,simplify=TRUE)

  if (mirt::extract.mirt(mirtfit, "ngroups") == 1) {

    a = as.data.frame(coefs$items)
    ismulti = "a2" %in% names(a)

    if (itemtype=="gpcm") { # Other Form for GPCM
      re =  list(a = a$a1,
                 d = as.matrix(a[grep("d",names(a))]),
                 g = rep(0,length(a$a1)),
                 itemtype=itemtype)
    } else if (ismulti) { # Other Form for mulidimensional
      re =  list(a = as.matrix(a[grep("a",names(a))]),
                 d = as.matrix(a[grep("d",names(a))]),
                 g = rep(0,length(a$a1)),
                 itemtype=itemtype)
    } else { # default form
      re =  list(a = a$a1,
                 d  = a$d,
                 g = a$g,
                 itemtype=itemtype)
    }


  } else { # multiple groups

    re = list()
    for (i in 1:length(coefs)) { # For all groups

      a = as.data.frame(coefs[[i]]$items)

      if (itemtype=="gpcm") { # Other Form for GPCM
        re[[i]] =  list(a = a$a1,
                        d=cbind(rep(0,length(a$a1)),a$d1,a$d2),
                        g = rep(0,length(a$a1)),
                        itemtype=itemtype)

      } else { # default form
        re[[i]] =  list(a = a$a1,
                        d  = a$d,
                        g = a$g,
                        itemtype=itemtype)
      }

    }
  }

  return(re)
}


#' Transform parameters to a longer format
#'
#' This is a helper function used to generate custom hypotheses. See the "adding_hypotheses" vignette.
#'
#'
#'
#' @param pars list of parameters. Can also be coefficients from a model fitted by mirt. In this case, the from.mirt argument has to be set to TRUE
#'
#' @param from.mirt logical, treat as coefficients from a model fitted by mirt if TRUE
#' @param itemtype character, type of the item as string, e.g. "2PL"
#'
#' @return numeric vector
#' @export
#'
#' @examples
#'
#' pars = list(a= c(1,1,1),d=c(0,0,0))
#' pars.long(pars,itemtype="2PL")
#'
pars.long = function(pars, itemtype, from.mirt=FALSE) {

  if(!from.mirt) {
    is.multi = "a2" %in% colnames(pars$a)

    if (itemtype=="2PL"&!is.multi) {
      re = rbind(pars$a,pars$d) %>% as.numeric()
    } else if (itemtype=="2PL"&is.multi) {
      re = cbind(pars$a,pars$d) %>% t() %>% as.numeric()
    } else if (itemtype=="3PL") {
      re = rbind(pars$a,pars$d,pars$g) %>% as.numeric()
    } else if (itemtype=="gpcm") {
      re = cbind(pars$a,pars$d[,2:ncol(pars$d)]) %>% t() %>% as.numeric()
    }
  }

  if(from.mirt) {

      if (itemtype=="2PL"& class(pars)=="MultipleGroupClass") {
        pars = lapply(mirt::coef(pars,simplify=TRUE),function(x) x$items[,1:2] %>% t() %>% as.numeric())
        re = c(pars$A,pars$B)

    } else if (itemtype=="2PL") {

        ps = mirt::coef(pars,simplify=TRUE)$items
        is.multi = "a2" %in% colnames(ps)

        if (!is.multi) {
        re = mirt::coef(pars,simplify=TRUE)$items[,1:2] %>% t() %>% as.numeric()
        }
        if (is.multi) {
          re = mirt::coef(pars,simplify=TRUE)$items[,1:3] %>% t() %>% as.numeric()
          re = re[-(length(re)-1)] # since it is not estimated
        }

    } else if (itemtype=="3PL") {

      re = mirt::coef(pars,simplify=TRUE)$items[,1:3] %>% t() %>% as.numeric()

    } else if (itemtype=="gpcm") {

      nkatx = max(mirt::extract.mirt(pars, "K") - 1)
      a1 = mirt::coef(pars,simplify=TRUE)$items
      re =  a1[,c(1,(ncol(a1)-nkatx+1):ncol(a1))] %>% t() %>% as.numeric()
    }
  }

  return(re)
}


get_ncp = function(chii,df) {

  #Extract ncp from vector of chi² distributed values
  #
  # Used in the sampling based ncps

  ncp_x = mean(chii)-df # mean for compatibility with chii as vector of multiple values
  chii[chii<=0] = .00000001
  if (ncp_x<0) {return(list(ncp=0,sd = NA))}
  return(list(ncp=ncp_x))
}



calc.N = function(hyp,ncp,power=.8,alpha=.05) {
  # Calculate sample size
  #
  # @param hyp Hypothesis Object created from the setup.hypothesis function
  # @param ncp numeric, Noncentrality parameter for n=1
  # @param power numeric, desired power, e.g. .8 (default)
  # @param alpha numeric, alpha niveau, e.g. .05 (default)
  #
  # @return
  # @export
  #
  # @examples
  #
  # dat <- expand.table(LSAT7)
  # mirtfit <- mirt(dat,1,verbose = FALSE)
  # hyp <- setup.hypothesis(type = "1PLvs2PL", altpars = mirtfit)
  # ncps <- calculate_ncps(hyp=hyp)
  # ssize(hyp=hyp,ncp=ncps,alpha=.05,power=.80)
  #

  df = nrow(hyp$resmod$Amat)
  qcentral <- stats::qchisq(p = 1 - alpha, df = df)
  func <-
    function (x) {
      1-  power - stats::pchisq(q = qcentral, df = df, ncp = x)
    }
  w0 <-
    stats::uniroot(
      f = func,
      interval = c(0, 1000),
      tol = .Machine$double.eps ^ 0.5
    )$root

  re <- as.numeric(ceiling(w0 / ncp))
  names(re) = names(ncp)
  return(re)
}



calc.power = function(hyp,ncp,ssize,alpha=.05) {
  # Calculate power
  #
  # @param hyp Hypothesis Object created from the setup.hypothesis function
  # @param ncp numeric, Noncentrality parameter for n=1
  # @param ssize interger, sample size
  # @param alpha numeric, alpha niveau, e.g. .05 (default)
  #
  # @return
  # @export
  #
  # @examples
  #
  # dat <- expand.table(LSAT7)
  # mirtfit <- mirt(dat,1,verbose = FALSE)
  # hyp <- setup.hypothesis(type = "1PLvs2PL", altpars = mirtfit)
  # ncps <- calculate_ncps(hyp=hyp)
  # power(hyp=hyp,ncp=ncps,alpha=.05,ssize=500)
  #


  df = nrow(hyp$resmod$Amat)
  crit = stats::qchisq(1-alpha,df = df, ncp = 0) %>% as.numeric()
  re = 1 - stats::pchisq(q = crit, df = df, ncp = ncp*ssize)
  return(re)
}
