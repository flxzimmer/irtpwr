

#' Calculate ncps (analytical/sampling based, Wald/LR/Score)
#'
#' Wrapper function for the different ncp functions
#'
#' @param hyp Hypothesis Object created from the setup_hypothesis function
#' @param stat character vector containing the statistics to be calculated. Options are "Wald","LR","Score". Also, "all" (default) for all three.
#' @param n integer, number of persons that the ncp refers to (default 1)
#' @param sampling boolean, TRUE gives the simulation based statistics instead of the analytical
#' @param n.pers passed to the ncp.sim function
#'
#' @return numerical vector containing the ncps
#' @export
#'
#' @examples
#' dat <- expand.table(LSAT7)
#' mirtfit <- mirt(dat,1)
#' pars = coef_short(mirtfit)
#' hyp <- setup_hypothesis(type = "1PLvs2PL", altpars = pars)
#' ncps <- calculate_ncps(hyp=hyp)
#'
calculate_ncps = function(hyp,stat=c("Wald","LR","Score","Gradient"),n=1,sampling=FALSE,sampling.npers = 10^4,approx.npers=10^4) {
  # implement ncp.sim for other than all three stats

  # if (stat[1]=="all") {stat= c("Wald","LR","Score","Gradient")}

  if (sampling) {
    data = setup.data(hyp,n = sampling.npers)
    fitted = mml.fit(hyp, data = data,infmat.unres = "ApproxFisher",infmat.res="ApproxFisher",approx.npers=approx.npers)
    obs =stat_obs(fitted,stat)
    ncps = sapply(obs,function(x) get_ncp(x,df=nrow(hyp$resmod$Amat))$ncp)
    re = ncps/sampling.npers
    }

  if (!sampling) {

    a = rep(NA,4) # Prepare Result Vector

    lx = NULL # Prepare empty lx Vector for Gradient stat

    # Calculate ML Restricted parameters centrally for use in Score and LR and Gradient
    if ( "LR" %in% stat | "Score" %in% stat | "Gradient" %in% stat ) {parsr = hyp$type$maximizeL(hyp)}

    if ("Wald" %in% stat) {a[1] = ncp.wald(hyp = hyp)}
    if ("LR" %in% stat) {a[2] = ncp.lr(hyp = hyp,parsr=parsr)}
    if ("Score" %in% stat) {
      sctemp = ncp.score(hyp = hyp,parsr=parsr)
      a[3] =sctemp$ncp
      lx = sctemp$lx
    }
    if ("Gradient" %in% stat) {a[4] = ncp.grad(hyp = hyp,parsr=parsr,lx=lx)}

    re = a[!is.na(a)]
    names(re) = stat
  }
  return(re)
}



#' Analytical ncp for the Wald statistic
#'
#' @param hyp Hypothesis Object created from the setup_hypothesis function
#' @param n integer, number of persons that the ncp refers to (default 1)
#'
#' @return
#' @export
#'
#' @examples
ncp.wald = function(hyp,n=1) {
  resmod = hyp$resmod
  unresmod = hyp$unresmod

  multigroup = isTRUE(unresmod$multigroup)

  # Fisher implementation is currently flawed in mirt:
  # if (multigroup) {method = "ApproxFisher"} else {method = "Fisher"}
  method = "Fisher"

  A = resmod$Amat
  dif = A%*%unresmod$longpars-resmod$cvec
  sigma = infmat(unresmod$parsets,method=method,model = hyp$unresmod$model, multigroup = multigroup,itemtype = hyp$unresmod$itemtype) %>% solve()
  re = t(dif) %*% solve(A%*% sigma %*% t(A)) %*% dif %>% c()
  return(re*n)
}


#' Analytical ncp for the LR statistic
#'
#' @param hyp Hypothesis Object created from the setup_hypothesis function
#' @param n integer, number of persons that the ncp refers to (default 1)
#' @param parsr Restricted parameters
#'
#' @return
#' @export
#'
#' @examples
ncp.lr = function(hyp,n=1,parsr= NULL) {

  resmod = hyp$resmod
  unresmod = hyp$unresmod

  if (is.null(parsr)) {
    parsr = hyp$type$maximizeL(hyp)
  }
  res = c()
  unres = c()
  load.functions(resmod$itemtype)

  if (isTRUE(unresmod$multigroup)) {

    pars1 = unresmod$parsets[[1]]
    pars2 = unresmod$parsets[[2]]
    n.kat = max(ncol(pars1$d),2)
    patterns = as.matrix(expand.grid(lapply(1:resmod$n.items,function(x) 0:(n.kat-1))))

    for (i in 1:nrow(patterns)) {
      gr = g(patterns[i,],parsr)
      g1 = g(patterns[i,],pars1)
      g2 = g(patterns[i,],pars2)
      res[i] = log(gr) * (.5 * g1 + .5 * g2)
      unres[i] = .5 * log(g1) * g1 + .5 * log(g2) * g2
    }

  } else {

    pars = unresmod$parsets
    n.kat = max(ncol(pars$d),2)
    patterns = as.matrix(expand.grid(lapply(1:resmod$n.items,function(x) 0:(n.kat-1))))

    for (i in 1:nrow(patterns)) {
      gr = g(patterns[i,],parsr)
      gu = g(patterns[i,],pars)
      res[i] = log(gr) * (gu)
      unres[i] = log(gu) * (gu)
    }

  }

  res = sum(res)
  unres = sum(unres)

  re = 2*(unres-res) %>% c()

  return(re*n)

}


#' Analytical ncp for the Score statistic
#'
#' @param hyp Hypothesis Object created from the setup_hypothesis function
#' @param n integer, number of persons that the ncp refers to (default 1)
#' @param parsr Restricted parameters
#'
#' @return
#' @export
#'
#' @examples
ncp.score = function(hyp,n=1,parsr=NULL) {
  resmod = hyp$resmod
  unresmod = hyp$unresmod

  if (is.null(parsr)) {
    parsr = hyp$type$maximizeL(hyp)
  }

  load.functions(resmod$itemtype)
  multigroup = isTRUE(unresmod$multigroup)

  # Fisher implementation is currently flawed in mirt:
  # if (multigroup) {method = "ApproxFisher"} else {method = "Fisher"}
  method = "Fisher"

  if (multigroup) { # Multigroup Model

    pars1 = unresmod$parsets[[1]]
    pars2 = unresmod$parsets[[2]]
    n.kat = max(ncol(pars1$d),2)
    patterns = lapply(1:resmod$n.items,function(x) 0:(n.kat-1)) %>% expand.grid() %>% as.matrix() %>% split(.,seq(nrow(.)))

    freq12 = lapply(patterns,function(x) c(g(x,pars1)/2,g(x,pars2)/2)) %>% do.call(rbind,.)
    freq = cbind(rowSums(freq12),freq12)
    ly = lapply(1:length(patterns),function(i) {
      # paste(100*i/length(patterns)%>%round(.,1),"% finished") %>% print()
      l = ldot(patterns[[i]],parsr); list(freq[i,1] * l,freq[i,2] * l,freq[i,3] * l)} )
    lx = lapply(ly,function(x) x[[1]]) %>% do.call(rbind,.) %>% colSums() %>% array(.,dim=c(length(.),1))
    lx1 = lapply(ly,function(x) x[[2]]) %>% do.call(rbind,.) %>% colSums() %>% array(.,dim=c(length(.),1))
    lx2 = lapply(ly,function(x) x[[3]]) %>% do.call(rbind,.) %>% colSums() %>% array(.,dim=c(length(.),1))

    lx = c(lx1[1:2,],lx[3:length(lx),],lx2[1:2,]) %>% array(.,dim=c(length(.),1))

    # transform restricted pars for infmat calculation
    parsr = list(parsr,parsr)

  } else { # Single Group Model

    # if(resmod$itemtype == "3PL") {method="customFisher"}

    pars = unresmod$parsets
    n.kat = max(ncol(pars$d),2)
    patterns = lapply(1:resmod$n.items,function(x) 0:(n.kat-1)) %>% expand.grid() %>% as.matrix() %>% split(.,seq(nrow(.)))
    i=0
    ly = lapply(patterns,function(x) {ldot(x,parsr) * g(x,pars)})
    lx = do.call(rbind,ly) %>% colSums() %>% array(.,dim=c(length(.),1))
  }

  sigma = infmat(parsr,method=method,model = hyp$unresmod$model, multigroup = multigroup,itemtype = hyp$unresmod$itemtype) %>% solve()

  ncp = t(lx) %*% sigma %*% lx %>% c()

  re = list(ncp = ncp*n,lx=lx)

  return(re)
}

#' Analytical ncp for the Score statistic
#'
#' @param hyp Hypothesis Object created from the setup_hypothesis function
#' @param n integer, number of persons that the ncp refers to (default 1)
#' @param parsr Restricted parameters
#'
#' @return
#' @export
#'
#' @examples
ncp.grad = function(hyp,n=1,parsr=NULL,lx=NULL) {
  resmod = hyp$resmod
  unresmod = hyp$unresmod

  if (is.null(parsr)) {
    parsr = hyp$type$maximizeL(hyp)
  }

  load.functions(resmod$itemtype)
  multigroup = isTRUE(unresmod$multigroup)

  if (multigroup & is.null(lx)) { # Multigroup Model

    pars1 = unresmod$parsets[[1]]
    pars2 = unresmod$parsets[[2]]
    n.kat = max(ncol(pars1$d),2)
    patterns = lapply(1:resmod$n.items,function(x) 0:(n.kat-1)) %>% expand.grid() %>% as.matrix() %>% split(.,seq(nrow(.)))

    freq12 = lapply(patterns,function(x) c(g(x,pars1)/2,g(x,pars2)/2)) %>% do.call(rbind,.)
    freq = cbind(rowSums(freq12),freq12)
    ly = lapply(1:length(patterns),function(i) {
      # paste(100*i/length(patterns)%>%round(.,1),"% finished") %>% print()
      l = ldot(patterns[[i]],parsr); list(freq[i,1] * l,freq[i,2] * l,freq[i,3] * l)} )
    lx = lapply(ly,function(x) x[[1]]) %>% do.call(rbind,.) %>% colSums() %>% array(.,dim=c(length(.),1))
    lx1 = lapply(ly,function(x) x[[2]]) %>% do.call(rbind,.) %>% colSums() %>% array(.,dim=c(length(.),1))
    lx2 = lapply(ly,function(x) x[[3]]) %>% do.call(rbind,.) %>% colSums() %>% array(.,dim=c(length(.),1))

    lx = c(lx1[1:2,],lx[3:length(lx),],lx2[1:2,]) %>% array(.,dim=c(length(.),1))

    # transform restricted pars for infmat calculation
    parsr = list(parsr,parsr)

  }

  if (!multigroup & is.null(lx)) { # Single Group Model

    pars = unresmod$parsets
    n.kat = max(ncol(pars$d),2)
    patterns = lapply(1:resmod$n.items,function(x) 0:(n.kat-1)) %>% expand.grid() %>% as.matrix() %>% split(.,seq(nrow(.)))
    i=0
    ly = lapply(patterns,function(x) {ldot(x,parsr) * g(x,pars)})
    lx = do.call(rbind,ly) %>% colSums() %>% array(.,dim=c(length(.),1))
  }

  A = resmod$Amat
  dif = A%*%unresmod$longpars-resmod$cvec

  lambda = findlambda(lx,A)

  # fn = function(lambda) {sum(abs(lx + t(A)%*%lambda))}
  # lambda = optim(rep(0,nrow(A)),fn)$par

  re = t(lambda) %*% dif %>% as.numeric() %>% abs()

  return(re*n)
}



