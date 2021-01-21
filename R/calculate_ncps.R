

#' Calculate ncps (analytical/sampling based, Wald/LR/Score)
#'
#' Wrapper function for the different ncp functions
#'
#' @param hyp Hypothesis Object created from the setup_hypothesis function
#' @param stat character vector containing the statistics to be calculated. Options are "Wald","LR","Score". Also, "all" (default) for all three.
#' @param n integer, number of persons that the ncp refers to (default 1)
#' @param simbased boolean, TRUE gives the simulation based statistics instead of the analytical
#' @param n.pers passed to the ncp.sim function
#' @param runs passed to the ncp.sim function
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
calculate_ncps = function(hyp,stat="all",n=1,simbased=FALSE,n.pers = 10^4,runs =1) {
  # implement ncp.sim for other than all three stats

  if (stat=="all") {stat= c("Wald","LR","Score")}

  if (simbased) {
    re = ncp.sim(hyp,stat=stat,n=n,n.pers=n.pers,runs=runs)
  }

  if (!simbased) {
    # Calculate ML Restricted parameters centrally for use in Score and LR
    if ( "LR" %in% stat | "Score" %in% stat ) {parsr = hyp$type$maximizeL(hyp)}

    if ("Wald" %in% stat) {a1 = ncp.wald(hyp = hyp);names(a1)="Wald"}
    if ("LR" %in% stat) {a2 = ncp.lr(hyp = hyp,parsr=parsr);names(a2)="LR"}
    if ("Score" %in% stat) {a3 = ncp.score(hyp = hyp,parsr=parsr);names(a3)="Score"}

    re = c(a1,a2,a3)
  }
  return(re)
}



#' Calculate sampling based ncps (Wrapper)
#'
#' @param hyp Hypothesis Object created from the setup_hypothesis function
#' @param stat character vector containing the statistics to be calculated. Options are "Wald","LR","Score". Also, "all" (default) for all three.
#' @param n integer, number of persons that the ncp refers to (default 1)
#' @param n.pers integer, Number of persons in the artificial datasets
#' @param runs integer, Number of artificial datasets
#'
#' @return numerical vector containing the ncps
#' @export
#'
#' @examples
#' dat <- expand.table(LSAT7)
#' mirtfit <- mirt(dat,1)
#' pars = coef_short(mirtfit)
#' hyp <- setup_hypothesis(type = "1PLvs2PL", altpars = pars)
#' ncps <- ncp.sim(hyp=hyp)
ncp.sim = function(hyp,stat="all",n=1,n.pers=NULL,runs=3) {

  resmod = hyp$resmod
  unresmod = hyp$unresmod

  waldncp =  ncp.sim.wald(hyp,simbased.npers=n.pers) #  simbased.npers=n.pers*100 (previous simruns)

  res = list()
  for (i in 1:runs) {
    # paste("run",i,"of",runs,"runs") %>% print()
    datasetx = setup.data(hyp,n = n.pers)
    fittedx = mml.fit(hyp, data = datasetx)
    # print("fitted, now calculating stats")
    res[[i]] = c(lr_obs(fittedx),score_obs(fittedx))
    # res[[i]] = c(1,score_obs(fittedx))

  }
  re = do.call(rbind,res) %>% apply(.,2,function(x) get_ncp(x,df=nrow(resmod$Amat))$ncp) %>% c()
  re = re/n.pers
  re = c(waldncp,re)
  return(re*n)
}


#' Sampling-based ncp for the Wald statistic
#'
#' Uses a shortcut not available for LR and score (see paper)
#'
#' @param hyp Hypothesis Object created from the setup_hypothesis function
#' @param n integer, number of persons that the ncp refers to (default 1)
#' @param simbased.npers integer, Number of persons used in the sampling based approach
#'
#' @return
#' @export
#'
#' @examples
ncp.sim.wald = function(hyp,n=1,simbased.npers) {
  resmod = hyp$resmod
  unresmod = hyp$unresmod
  method="mirtOakessim"
  if (isTRUE(unresmod$multigroup)) {method="mirtOakesmultigroupsim"}

  A = resmod$Amat
  dif = A%*%unresmod$longcoef-resmod$cvec
  # print("calculating fisher expected infmat")
  sigma = infmat(unresmod,method=method,simbased.npers=simbased.npers) %>% solve()
  # print("finished calculating fisher expected infmat")
  re = t(dif) %*% solve(A%*% sigma %*% t(A)) %*% dif %>% c()
  return(re*n)
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

  method="mirtFisher"
  if (isTRUE(unresmod$multigroup)) {method="mirtFishermultigroup"}

  # if (isTRUE(unresmod$multigroup) & method=="mirtFisher") {method="mirtFishermultigroup"}
  # if (isTRUE(unresmod$multigroup) & method=="mirtOakessim") {method="mirtOakesmultigroupsim"}

  A = resmod$Amat
  dif = A%*%unresmod$longcoef-resmod$cvec
  # print("calculating fisher expected infmat")
  sigma = infmat(unresmod,method=method) %>% solve()
  # print("finished calculating fisher expected infmat")
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

    pars1 = unresmod[[1]][[1]]
    pars2 = unresmod[[1]][[2]]
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

    pars = unresmod
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

  if (isTRUE(unresmod$multigroup)) {

    method = "mirtFishermultigroup"

    pars1 = unresmod[[1]][[1]]
    pars2 = unresmod[[1]][[2]]
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

  } else {

    if(resmod$itemtype == "3PL") {method="customFisher"} else {method = "mirtFisher"}

    pars = unresmod
    n.kat = max(ncol(pars$d),2)
    patterns = lapply(1:resmod$n.items,function(x) 0:(n.kat-1)) %>% expand.grid() %>% as.matrix() %>% split(.,seq(nrow(.)))
    i=0
    ly = lapply(patterns,function(x) {ldot(x,parsr) * g(x,pars)})
    lx = do.call(rbind,ly) %>% colSums() %>% array(.,dim=c(length(.),1))
  }

  sigma = infmat(parsr,method=method) %>% solve()

  re = t(lx) %*% sigma %*% lx %>% c()

  return(re*n)
}




