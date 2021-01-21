# Statistics --------------------------------------------------------------

#' Fit restricted and unrestriced model by ML using mirt
#'
#' @param hyp Hypothesis Object created from the setup_hypothesis function
#' @param data dataframe with data to be fitted
#' @param infmat.method Method used for calculating the infmat used for the variance-covariance matrix, e.g. "Oakes" (default)
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
#' fitted <- mml.fit(data = data,hyp = hyp)
#'
mml.fit = function(hyp,data,infmat.method = "Oakes",free_mean = FALSE) {

      if (is.data.frame(data)) {

  } else {
    group = data$group
    data = data$data
  }


  if (isTRUE(free_mean)) {invariance = "free_mean"} else {invariance = ""}

  if (isTRUE(hyp$unresmod$multigroup)) { # Use a different function for multigroup models
    unres = mirt::multipleGroup(data,model = hyp$unresmod$model,itemtype = hyp$unresmod$itemtype,SE = TRUE,SE.type = infmat.method,technical = list(NCYCLES = 5000),verbose = FALSE,group=group,invariance=invariance)
  } else {
    unres =  mirt::mirt(data,model = hyp$unresmod$model,itemtype = hyp$unresmod$itemtype,SE = TRUE,SE.type = infmat.method,technical = list(NCYCLES = 5000),verbose = FALSE)
  }

  if (isTRUE(hyp$resmod$multigroup)) {
    res = mirt::multipleGroup(data,model = hyp$resmod$model,itemtype = hyp$resmod$itemtype,technical = list(NCYCLES = 5000),verbose = FALSE,group=group,invariance=invariance)
  } else {
    res =  mirt::mirt(data,model = hyp$resmod$model,itemtype = hyp$resmod$itemtype,technical = list(NCYCLES = 5000),verbose = FALSE)
  }

    re = list(unres=unres,res=res,hyp=hyp)
  return(re)
}


#' Calculate Wald statistic
#'
#' @param fitted Object created from mml.fit
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
#' fitted <- mml.fit(data = data,hyp = hyp)
#' wald_obs(fitted)
#'
wald_obs = function(fitted) {
  resmod = fitted$hyp$resmod
  pars=coef.long(mirtfit = fitted$unres,itemtype = resmod$itemtype)
  A = resmod$Amat

  dif = A%*%pars-resmod$cvec
  sigma=mirt::vcov(fitted$unres)
  browser()
  re = t(dif) %*% solve(A%*% sigma %*% t(A)) %*% dif %>% as.numeric()
  return(re)
}

#' Calculate LR statistic
#'
#' @param fitted Object created from mml.fit
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
#' fitted <- mml.fit(data = data,hyp = hyp)
#' lr_obs(fitted)
#'
lr_obs = function(fitted) {
  re = 2*(mirt::logLik(fitted$unres)-mirt::logLik(fitted$res))
  return(re)
}


#' Calculate Score statistic
#'
#' @param fitted Object created from mml.fit
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
#' fitted <- mml.fit(data = data,hyp = hyp)
#' score_obs(fitted)
#'
score_obs = function(fitted) {

  resmod = fitted$hyp$resmod
  unresmod = fitted$hyp$unresmod
  parsr = coef_short(fitted$res,itemtype = resmod$itemtype)
  load.functions(resmod$itemtype)
  patterns = fitted$res@Data$data
  hashs = apply(patterns,1,digest::digest) %>% factor(.,levels=unique(.))
  rownames(patterns)=hashs
  up = unique(patterns)
  freq = table(hashs)

  if (isTRUE(unresmod$multigroup))  {

    freq1 = table(hashs[1:(nrow(patterns)/2)])
    freq2 = table(hashs[(nrow(patterns)/2+1):nrow(patterns)])

    ly = lapply(rownames(up),function(i) {l = ldot(up[i,],parsr); list(freq[i] * l,freq1[i] * l,freq2[i] * l)} )
    lx = lapply(ly,function(x) x[[1]]) %>% do.call(rbind,.) %>% colSums() %>% array(.,dim=c(length(.),1))
    lx1 = lapply(ly,function(x) x[[2]]) %>% do.call(rbind,.) %>% colSums() %>% array(.,dim=c(length(.),1))
    lx2 = lapply(ly,function(x) x[[3]]) %>% do.call(rbind,.) %>% colSums() %>% array(.,dim=c(length(.),1))

    lx = c(lx1[1:2,],lx[3:length(lx),],lx2[1:2,]) %>% array(.,dim=c(length(.),1))

    method = "mirtOakesmultigroup"

  } else {

    ly = lapply(rownames(up),function(i) ldot(up[i,],parsr)*freq[i] )
    lx = do.call(rbind,ly) %>% colSums() %>% array(.,dim=c(length(.),1))

    if (resmod$itemtype == "3PL") {method="customFisher"} else {method = "mirtOakes"}
  }

  sigma = (infmat(parsr,method=method,resmod = resmod,data=patterns)*nrow(patterns)) %>% solve()


  re = t(lx) %*% sigma %*% lx %>% as.numeric()
  return(re)
}


