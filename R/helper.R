
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

  #don't output g parameter if they are all zero?

  a = as.data.frame(mirt::coef(mirtfit,simplify=TRUE)$items)
  if (is.null(itemtype)) {
    itemtype=mirtfit@Model$itemtype[1]
  }

  re =  list(a = a$a1, d  = a$d,g = a$g,itemtype=itemtype)
  if (itemtype=="gpcm") {
    re =  list(a = a$a1, d=cbind(rep(0,length(a$a1)),a$d1,a$d2),g = rep(0,length(a$a1)),itemtype=itemtype)
  }

  return(re)
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
coef.long = function(mirtfit=NULL,itemtype=NULL,shortcoef=NULL) {
  # extracts coefs from mml coef output
  if (is.null(shortcoef)) {
    if (itemtype=="2PL"& mirtfit@Call[[1]]=="mirt::multipleGroup") {
      pars = lapply(mirt::coef(mirtfit,simplify=TRUE),function(x) x$items[,1:2] %>% t() %>% as.numeric())
      re = c(pars$A,pars$B[1:2])
    } else if (itemtype=="2PL") {
      re = mirt::coef(mirtfit,simplify=TRUE)$items[,1:2] %>% t() %>% as.numeric()
    } else if (itemtype=="3PL") {
      re = mirt::coef(mirtfit,simplify=TRUE)$items[,1:3] %>% t() %>% as.numeric()
    } else if (itemtype=="gpcm") {
      nkatx = max(mirtfit@Data$data)
      a1 = mirt::coef(mirtfit,simplify=TRUE)$items
      re =  a1[,c(1,(ncol(a1)-nkatx+1):ncol(a1))] %>% t() %>% as.numeric()
    }
  } else {
    if (itemtype=="2PL") {
      re = rbind(shortcoef$a,shortcoef$d) %>% as.numeric()
    } else if (itemtype=="3PL") {
      re = rbind(shortcoef$a,shortcoef$d,shortcoef$g) %>% as.numeric()
    } else if (itemtype=="gpcm") {
      re = cbind(shortcoef$a,shortcoef$d[,2:ncol(shortcoef$d)]) %>% t() %>% as.numeric()
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



