

logit <- function(p) { log(p/(1-p)) }
logitinv<- function(q) { 1/(1+exp(-q)) }



#' Calculate the time needed
#'
#' @param hyp
#' @param n.items
#'
#' @return
#' @export
#'
#' @examples
calctime = function(hyp,n.items) {

  n.items.ex <- length(hyp$unresmod$parsets$d)

  ptm <- proc.time() #start timer
  ncps <- calculate_ncps(hyp=hyp)
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
#' pars <- coef_short(mirtfit)
#'
coef_short = function(mirtfit,itemtype=NULL) {
  # extracts coefs from mml coef output

  # vllt nützlich? mirtfit@Internals[["shortpars"]]

  if (is.null(itemtype)) {
    itemtype <- extract.mirt(mirtfit, what = "itemtype")[1]
    }

  coefs <- mirt::coef(mirtfit,simplify=TRUE)

  if (extract.mirt(mirtfit, "ngroups") == 1) {

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


#' extract coefs from mirtfit or element created by shortcoef
#'
#' @param pars Parameter Set
#' @param from.mirt Boolean, is the parameter set a mirt model?
#' @param itemtype optional, itemtype as string
#'
#' @return
#' @export
#'
#' @examples
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

    # if (itemtype=="2PL"& class(pars)=="MultipleGroupClass" & pars@Model$nest== pars@Model$nestpars - pars@Data$nitems -5) {
    #
    #     pars = lapply(mirt::coef(pars,simplify=TRUE),function(x) x$items[,1:2] %>% t() %>% as.numeric())
    #     re = c(pars$A,pars$B[seq(12,length(pars$B),2)])
    #
    # } else

      if (itemtype=="2PL"& class(pars)=="MultipleGroupClass") {
        pars = lapply(mirt::coef(pars,simplify=TRUE),function(x) x$items[,1:2] %>% t() %>% as.numeric())
        re = c(pars$A,pars$B[1:2])

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

      nkatx = max(extract.mirt(pars, "K") - 1)
      # nkatx = max(pars@Data$data)
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
#' mirtfit <- mirt(dat,1,verbose = FALSE)
#' hyp <- setup_hypothesis(type = "1PLvs2PL", altpars = mirtfit)
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
#' mirtfit <- mirt(dat,1,verbose = FALSE)
#' hyp <- setup_hypothesis(type = "1PLvs2PL", altpars = mirtfit)
#' ncps <- calculate_ncps(hyp=hyp)
#' power(hyp=hyp,ncp=ncps,alpha=.05,ssize=500)
#'
power = function(hyp,ncp,ssize,alpha=.05) {
  df = nrow(hyp$resmod$Amat)
  crit = qchisq(1-alpha,df = df, ncp = 0) %>% as.numeric()
  re = 1 - pchisq(q = crit, df = df, ncp = ncp*ssize)
  return(re)
}


#' Title
#'
#' @param e
#' @param x
#'
#' @return
#' @export
#'
#' @examples
save.error = function(e,x) {
  num = sample(1:1000,1)
  file_string =paste0("errors/err_",num,".Rdata")
  save(e,x,file=file_string)
}

