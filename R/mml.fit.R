

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
mml.fit = function(hyp,data,infmat.unres = "Fisher",infmat.res="Fisher",free_mean = FALSE,approx.npers=10^6) {
  # For the unrestricted model, you have the choice between observed and expected infomat
  # Currently, only expected infmats are implemented for the unrestricted models, these are more accurate and can also be calculated for large numbers of items (ApproxFisher)
  # For the restricted model, usually only the expected infmat makes sense since it has a special form depending on the unrestricted model and is therefore calculated seperately

  if (!is.data.frame(data)) {
    group = data$group
    data = data$data
  }

  multigroup.unres = isTRUE(hyp$unresmod$multigroup)
  multigroup.res = isTRUE(hyp$resmod$multigroup)

  if (isTRUE(free_mean)) {invariance = "free_mean"} else {invariance = ""}

  # if (infmat.unres == "ApproxFisher") { # Then Infmat for the unrestricted model is NOT calculated by mirt
  #   SE.calc = FALSE
  #   SE.type = "Fisher" # placeholder to avoid error
  # } else { #Infmat for the unrestricted model is calculated by mirt (Oakes or Fisher)
  #   SE.calc = TRUE
  #   SE.type = infmat.unres
  # }

  if (multigroup.unres) {
    unres = mirt::multipleGroup(data,model = hyp$unresmod$model,itemtype = hyp$unresmod$itemtype,technical = list(NCYCLES = 5000),verbose = FALSE,group=group,invariance=invariance)
  }

  if (!multigroup.unres) {
    unres =  mirt::mirt(data,model = hyp$unresmod$model,itemtype = hyp$unresmod$itemtype,technical = list(NCYCLES = 5000),verbose = FALSE)
  }

  if (multigroup.res) {
    res = mirt::multipleGroup(data,model = hyp$resmod$model,itemtype = hyp$resmod$itemtype,technical = list(NCYCLES = 5000),verbose = FALSE,group=group,invariance=invariance)
  }

  if (!multigroup.res) {
    res =  mirt::mirt(data,model = hyp$resmod$model,itemtype = hyp$resmod$itemtype,technical = list(NCYCLES = 5000),verbose = FALSE)
    }

  # Caluclating Infmats

  parsr = coef_short(res,itemtype = hyp$resmod$itemtype)
  pars = coef_short(unres,itemtype = hyp$unresmod$itemtype)

    unres@vcov = (infmat(pars,method=infmat.unres,approx.npers = approx.npers,multigroup = multigroup.unres,model = hyp$unresmod$model) * nrow(data) )%>% solve()

    res@vcov = (infmat(parsr,method=infmat.unres,approx.npers = approx.npers,multigroup = multigroup.unres,model = hyp$unresmod$model) * nrow(data) )%>% solve()


    re = list(unres=unres,res=res,hyp=hyp)
  return(re)
}


#' Title
#'
#' @param fitted
#'
#' @return
#' @export
#'
#' @examples
stat_obs = function(fitted,stat=c("Wald","LR","Score","Gradient")) {

  a = rep(NA,4) # Prepare Result Vector

  lx = NULL # Prepare empty lx Vector for Gradient stat

  if ("Wald" %in% stat) {a[1] = wald_obs(fitted)}
  if ("LR" %in% stat) {a[2] = lr_obs(fitted)}
  if ("Score" %in% stat) {
    sctemp = score_obs(fitted,export.lx = TRUE)
    a[3] =sctemp$stat
    lx = sctemp$lx
  }
  if ("Gradient" %in% stat) {a[4] = grad_obs(fitted,lx=lx)}

  re = a[!is.na(a)]
  names(re) = stat

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
  pars=pars.long(pars = fitted$unres,itemtype = resmod$itemtype,from.mirt=TRUE)
  A = resmod$Amat

  dif = A%*%pars-resmod$cvec
  sigma=mirt::vcov(fitted$unres)
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
score_obs = function(fitted,export.lx=FALSE) {

  resmod = fitted$hyp$resmod
  unresmod = fitted$hyp$unresmod
  parsr = coef_short(fitted$res,itemtype = resmod$itemtype)
  load.functions(resmod$itemtype)
  patterns = fitted$res@Data$data
  hashs = apply(patterns,1,digest::digest) %>% factor(.,levels=unique(.))
  rownames(patterns)=hashs
  up = unique(patterns)
  freq = table(hashs)

  if (isTRUE(resmod$multigroup))  { # multigroup model

    freq1 = table(hashs[1:(nrow(patterns)/2)])
    freq2 = table(hashs[(nrow(patterns)/2+1):nrow(patterns)])

    parsr = parsr[[1]]

    ly = lapply(rownames(up),function(i) {l = ldot(up[i,],parsr); list(freq[i] * l,freq1[i] * l,freq2[i] * l)} )
    lx = lapply(ly,function(x) x[[1]]) %>% do.call(rbind,.) %>% colSums() %>% array(.,dim=c(length(.),1))
    lx1 = lapply(ly,function(x) x[[2]]) %>% do.call(rbind,.) %>% colSums() %>% array(.,dim=c(length(.),1))
    lx2 = lapply(ly,function(x) x[[3]]) %>% do.call(rbind,.) %>% colSums() %>% array(.,dim=c(length(.),1))

    lx = c(lx1[1:2,],lx[3:length(lx),],lx2[1:2,]) %>% array(.,dim=c(length(.),1))


  } else { # single group model

    ly = lapply(rownames(up),function(i) ldot(up[i,],parsr)*freq[i] )
    lx = do.call(rbind,ly) %>% colSums() %>% array(.,dim=c(length(.),1))
  }

  sigma = mirt::vcov(fitted$res)

  re = t(lx) %*% sigma %*% lx %>% as.numeric()

  if (isTRUE(export.lx)) {
    re = list(stat=re,lx=lx)
  }

#   #alternative expression of score statistic
#   A = resmod$Amat
#   lambda1= findlambda(lx,A,v=1)
#   lambda2 = findlambda(lx,A,v=2)
#   re1 =  t(lambda1) %*% (A %*% sigma %*% t(A)) %*% lambda1 %>% as.numeric()
#   re2 =  t(lambda2) %*% (A %*% sigma %*% t(A)) %*% lambda2 %>% as.numeric()
#   re;re1;re2
# browser()

  return(re)
}



#' Calculate Gradient statistic
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
#' grad_obs(fitted)
#'
grad_obs = function(fitted,lx=NULL) {

  resmod = fitted$hyp$resmod
  unresmod = fitted$hyp$unresmod
  parsr = coef_short(fitted$res,itemtype = resmod$itemtype)
  load.functions(resmod$itemtype)
  patterns = fitted$res@Data$data
  hashs = apply(patterns,1,digest::digest) %>% factor(.,levels=unique(.))
  rownames(patterns)=hashs
  up = unique(patterns)
  freq = table(hashs)

  multigroup = isTRUE(unresmod$multigroup)

  if (multigroup & is.null(lx)) { # Multigroup Model

    freq1 = table(hashs[1:(nrow(patterns)/2)])
    freq2 = table(hashs[(nrow(patterns)/2+1):nrow(patterns)])

    parsr = parsr[[1]]

    ly = lapply(rownames(up),function(i) {l = ldot(up[i,],parsr); list(freq[i] * l,freq1[i] * l,freq2[i] * l)} )
    lx = lapply(ly,function(x) x[[1]]) %>% do.call(rbind,.) %>% colSums() %>% array(.,dim=c(length(.),1))
    lx1 = lapply(ly,function(x) x[[2]]) %>% do.call(rbind,.) %>% colSums() %>% array(.,dim=c(length(.),1))
    lx2 = lapply(ly,function(x) x[[3]]) %>% do.call(rbind,.) %>% colSums() %>% array(.,dim=c(length(.),1))

    lx = c(lx1[1:2,],lx[3:length(lx),],lx2[1:2,]) %>% array(.,dim=c(length(.),1))

  }

  if (!multigroup & is.null(lx)) { # Single Group Model
    ly = lapply(rownames(up),function(i) ldot(up[i,],parsr)*freq[i] )
    lx = do.call(rbind,ly) %>% colSums() %>% array(.,dim=c(length(.),1))
  }

  pars=pars.long(pars = fitted$unres,itemtype = resmod$itemtype,from.mirt=TRUE)
  A = resmod$Amat
  dif = (A%*%pars-resmod$cvec)

  lambda = findlambda(lx,A)

  # fn = function(lambda) {sum((lx + t(A)%*%lambda)^2)}
  # lambda = optim(rep(0,nrow(A)),fn)$par

  re = t(lambda) %*% dif %>% as.numeric() %>% abs()

  # fn = function(lambda) {sum(abs(lx + t(A)%*%lambda))}
  # lambda = optim(rep(0,nrow(A)),fn)$par
  #
  # re = t(lambda) %*% dif %>% as.numeric() %>% abs()
  # re
  return(re)
}


