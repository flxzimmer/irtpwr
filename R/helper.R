


calctime = function(hyp,n.items) {

  n.items.ex = length(hyp$unresmod$parsets$d)

  ptm <- proc.time() #start timer
  ncps <- calculate_ncps(hyp=hyp)
  time = proc.time() - ptm #stop timer

  time = as.numeric(time[3])

  # time is t * 2^3
  re = time *2^(n.items-n.items.ex)

  return(re)
}



findlambda = function(lx,A,v=2) {

  if (v==1) {
  fn = function(lambda) {sum(abs(lx + t(A)%*%lambda))}
  }
  if (v>=2) {
    fn = function(lambda) {sum((lx + t(A)%*%lambda)^v)}
  }
  lambda = optim(rep(0,nrow(A)),fn,method="BFGS")$par
  # lambda = optim(rep(0,nrow(A)),fn,method="Nelder-Mead")$par
  # t(lambda) %*% (A %*% sigma %*% t(A)) %*% lambda %>% as.numeric()

  return(lambda)
}

#' extract coefs from mirtfit
#'
#' @param mirtfit object created from mirt()
#' @param itemtype optional, itemtype as string
#'
#' @return
#' @export
#'
#' @examples
#'
#' dat <- expand.table(LSAT7)
#' mirtfit <- mirt(dat,1)
#' pars = coef_short(mirtfit)
#'
coef_short = function(mirtfit,itemtype=NULL) {
  # extracts coefs from mml coef output

  if (is.null(itemtype)) {
    itemtype=mirtfit@Model$itemtype[1]
  }

  coefs = mirt::coef(mirtfit,simplify=TRUE)

  if (mirtfit@Data$ngroups==1) {

    a = as.data.frame(coefs$items)

    if (itemtype=="gpcm") { # Other Form for GPCM
      re =  list(a = a$a1, d=cbind(rep(0,length(a$a1)),a$d1,a$d2),g = rep(0,length(a$a1)),itemtype=itemtype)

    } else { # default form
      re =  list(a = a$a1, d  = a$d,g = a$g,itemtype=itemtype)
    }


  } else { # multiple groups

    re = list()
    for (i in 1:length(coefs)) { # For all groups

      a = as.data.frame(coefs[[i]]$items)

      if (itemtype=="gpcm") { # Other Form for GPCM
        re[[i]] =  list(a = a$a1, d=cbind(rep(0,length(a$a1)),a$d1,a$d2),g = rep(0,length(a$a1)),itemtype=itemtype)

      } else { # default form
        re[[i]] =  list(a = a$a1, d  = a$d,g = a$g,itemtype=itemtype)
      }

    }
  }

  # if (length(re)==1) {   # shorten if possible
  #   re = re[[1]]
  # }


  return(re)

  # a = as.data.frame(mirt::coef(mirtfit,simplify=TRUE)$items)
  # if (is.null(itemtype)) {
  #   itemtype=mirtfit@Model$itemtype[1]
  # }
  #
  # re =  list(a = a$a1, d  = a$d,g = a$g,itemtype=itemtype)
  # if (itemtype=="gpcm") {
  #   re =  list(a = a$a1, d=cbind(rep(0,length(a$a1)),a$d1,a$d2),g = rep(0,length(a$a1)),itemtype=itemtype)
  # }
  #
  # return(re)
}


#' extract coefs from mirtfit or element created by shortcoef
#'
#' @param mirtfit object created from mirt
#' @param itemtype optional, itemtype as string
#' @param shortcoef alternative, object created from shortcoef
#'
#' @return
#' @export
#'
#' @examples
pars.long = function(pars,itemtype,from.mirt=FALSE) {

  if(!from.mirt) {
    if (itemtype=="2PL") {
      re = rbind(pars$a,pars$d) %>% as.numeric()
    } else if (itemtype=="3PL") {
      re = rbind(pars$a,pars$d,pars$g) %>% as.numeric()
    } else if (itemtype=="gpcm") {
      re = cbind(pars$a,pars$d[,2:ncol(pars$d)]) %>% t() %>% as.numeric()
    }
  }

  if(from.mirt) {
    if (itemtype=="2PL"& pars@Call[[1]]=="mirt::multipleGroup") {
      pars = lapply(mirt::coef(pars,simplify=TRUE),function(x) x$items[,1:2] %>% t() %>% as.numeric())
      re = c(pars$A,pars$B[1:2])
    } else if (itemtype=="2PL") {
      re = mirt::coef(pars,simplify=TRUE)$items[,1:2] %>% t() %>% as.numeric()
    } else if (itemtype=="3PL") {
      re = mirt::coef(pars,simplify=TRUE)$items[,1:3] %>% t() %>% as.numeric()
    } else if (itemtype=="gpcm") {
      nkatx = max(pars@Data$data)
      a1 = mirt::coef(pars,simplify=TRUE)$items
      re =  a1[,c(1,(ncol(a1)-nkatx+1):ncol(a1))] %>% t() %>% as.numeric()
    }
  }

  return(re)
}



#' Extract ncp from vector of chi² distributed values
#'
#' Used in the sampling based ncps
#'
#' @param chii numeric vector
#' @param df integer, degrees of freedom
#'
#' @return
#' @export
#'
#' @examples
get_ncp = function(chii,df) {
  ncp_x = mean(chii)-df
  chii[chii<=0] = .00000001
  if (ncp_x<0) {return(list(ncp=0,sd = NA))}
  # chi_ncp <- MASS::fitdistr(chii,"chi-squared",start=list(ncp=0),
  #                     method="Brent",df=df,lower=0,upper=10000000)
  # return(list(ncp=chi_ncp$estimate,sd = chi_ncp$sd))
  return(list(ncp=ncp_x))
}



#' Calculate sample size
#'
#' @param hyp Hypothesis Object created from the setup_hypothesis function
#' @param ncp numeric, Noncentrality parameter for n=1
#' @param power numeric, desired power, e.g. .8 (default)
#' @param alpha numeric, alpha niveau, e.g. .05 (default)
#'
#' @return
#' @export
#'
#' @examples
#'
#' dat <- expand.table(LSAT7)
#' mirtfit <- mirt(dat,1)
#' pars = coef_short(mirtfit)
#' hyp <- setup_hypothesis(type = "1PLvs2PL", altpars = pars)
#' ncps <- calculate_ncps(hyp=hyp)
#' ssize(hyp=hyp,ncp=ncps,alpha=.05,power=.80)
#'
ssize = function(hyp,ncp,power=.8,alpha=.05) {
  df = nrow(hyp$resmod$Amat)
  qcentral <- qchisq(p = 1 - alpha, df = df)
  func <-
    function (x) {
      1-  power - pchisq(q = qcentral, df = df, ncp = x)
    }
  w0 <-
    uniroot(
      f = func,
      interval = c(0, 1000),
      tol = .Machine$double.eps ^ 0.5
    )$root

  re <- as.numeric(ceiling(w0 / ncp))
  names(re) = names(ncp)
  return(re)
}



#' Calculate power
#'
#' @param hyp Hypothesis Object created from the setup_hypothesis function
#' @param ncp numeric, Noncentrality parameter for n=1
#' @param ssize interger, sample size
#' @param alpha numeric, alpha niveau, e.g. .05 (default)
#'
#' @return
#' @export
#'
#' @examples
#'
#' dat <- expand.table(LSAT7)
#' mirtfit <- mirt(dat,1)
#' pars = coef_short(mirtfit)
#' hyp <- setup_hypothesis(type = "1PLvs2PL", altpars = pars)
#' ncps <- calculate_ncps(hyp=hyp)
#' power(hyp=hyp,ncp=ncps,alpha=.05,ssize=500)
#'
power = function(hyp,ncp,ssize,alpha=.05) {
  df = nrow(hyp$resmod$Amat)
  crit = qchisq(1-alpha,df = df, ncp = 0) %>% as.numeric()
  re = 1 - pchisq(q = crit, df = df, ncp = ncp*ssize)
  return(re)
}



