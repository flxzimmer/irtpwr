

#' Calculate information matrix (expected / observed)
#'
#' @param pars Parameter Set
#' @param method string, method to be used, e.g. "Fisher" for the Fisher-expected information matrix calculated using the mirt package or "ApproxFisher" for its approximation
#' @param approx.npers integer, number of persons in the sampling based approach
#' @param data dataset for observed matrix methods
#' @param multigroup Boolean, multigroup model?
#' @param model Specific mirt model
#' @param itemtype String, itemtype
#'
#' @return
#' @export
#'
#' @examples
infmat = function(pars,method = "Fisher",approx.npers = 10^6,data=NULL,multigroup = FALSE,model = 1,itemtype="2PL",NCYCLES=5000,SE.type="Oakes") {

  if(is.null(pars$g)) {pars$g = rep(0,length(pars$a))}

  is.multi = "a2" %in% colnames(pars$a)


  if (method=="Fisher" & !multigroup) {
    dat = mirt::simdata(a = pars$a,d = pars$d,guess=pars$g,N =1000,itemtype = itemtype)
    synt = mirt::mirt(dat,1,itemtype = itemtype,technical = list(NCYCLES = 1),pars="values")
    apars = which(synt$name=="a1")
    synt$lbound[apars] =  synt$ubound[apars] = pars$a
    dpars = which(synt$name=="d")
    synt$lbound[dpars] = synt$ubound[dpars] = pars$d
    gpars = which(synt$name=="g")
    synt$lbound[gpars] = synt$ubound[gpars] = pars$g
    mml = mirt::mirt(dat,model=model,itemtype = itemtype,SE = TRUE,SE.type = "Fisher",technical = list(NCYCLES = NCYCLES),verbose = FALSE,pars=synt)
  }

  if (method=="Fisher" & multigroup) {
    pars1 = pars[[1]]
    pars2 = pars[[2]]
    df1 = mirt::simdata(a = pars1$a, d = pars1$d, N = 500, itemtype = itemtype)
    df2 = mirt::simdata(a = pars2$a, d = pars2$d, N = 500, itemtype = itemtype)
    group= rep(c("A","B"),each=500)
    dat =rbind(df1,df2)

    synt = mirt::multipleGroup(dat,group = group,itemtype = itemtype,model=model,technical = list(NCYCLES = 1),pars="values")
    apars = which(synt$name=="a1")
    dpars = which(synt$name=="d")
    synt$lbound[apars] =  synt$ubound[apars] = c(pars1$a,pars2$a)
    synt$lbound[dpars] = synt$ubound[dpars] = c(pars1$d,pars2$d)
    mml = mirt::multipleGroup(dat,model=model,itemtype = itemtype,SE = TRUE,SE.type = "Fisher",technical = list(NCYCLES = NCYCLES),verbose = FALSE,pars=synt,group=group)
  }

  if (method=="ApproxFisher"& !multigroup & itemtype!="gpcm" & !is.multi) {
    dat = mirt::simdata(a = pars$a,d = pars$d,guess=pars$g,N =approx.npers,itemtype = itemtype)
    synt = mirt::mirt(dat,1,itemtype = itemtype,technical = list(NCYCLES = 1),pars="values")
    apars = which(synt$name=="a1")
    synt$lbound[apars] =  synt$ubound[apars] = pars$a
    dpars = which(synt$name=="d")
    synt$lbound[dpars] = synt$ubound[dpars] = pars$d
    gpars = which(synt$name=="g")
    synt$lbound[gpars] = synt$ubound[gpars] = pars$g
    # if(pars$d[1]==1)browser()
    mml = mirt::mirt(dat,1,itemtype = itemtype,SE = TRUE,SE.type = SE.type,technical = list(NCYCLES = NCYCLES),verbose = FALSE,pars=synt)
  }

  if (method=="ApproxFisher"& !multigroup & itemtype!="gpcm" & is.multi) {
    dat = mirt::simdata(a = pars$a,d = pars$d,guess=pars$g,N =approx.npers,itemtype = itemtype)
    synt = mirt::mirt(dat,model,itemtype = itemtype,technical = list(NCYCLES = 1),pars="values")
    apars = grep("a",synt$name)
    synt$lbound[apars] =  synt$ubound[apars] = as.numeric(t(pars$a))
    dpars = which(synt$name=="d")
    synt$lbound[dpars] = synt$ubound[dpars] = pars$d
    gpars = which(synt$name=="g")
    synt$lbound[gpars] = synt$ubound[gpars] = pars$g
    mml = mirt::mirt(dat,model,itemtype = itemtype,SE = TRUE,SE.type = SE.type,technical = list(NCYCLES = NCYCLES),verbose = FALSE,pars=synt)
  }

  if (method=="ApproxFisher"& !multigroup & itemtype=="gpcm"& !is.multi) {
    # ggf funktioniert dieser Code auch für itemtype!="gpcm", ausprobieren!

    # browser()

    dat = mirt::simdata(a = pars$a,d = pars$d,guess=pars$g,N =approx.npers,itemtype = itemtype)
    synt = mirt::mirt(dat,1,itemtype = itemtype,technical = list(NCYCLES = 1),pars="values")
    apars = which(synt$name=="a1")
    synt$lbound[apars] =  synt$ubound[apars] = pars$a
    dpars = grep("d",synt$name)
    synt$lbound[dpars] = synt$ubound[dpars] = as.numeric(t(pars$d))
    gpars = which(synt$name=="g")
    synt$lbound[gpars] = synt$ubound[gpars] = pars$g
    mml = mirt::mirt(dat,1,itemtype = itemtype,SE = TRUE,SE.type = SE.type,technical = list(NCYCLES = NCYCLES),verbose = FALSE,pars=synt)
  }

  if (method=="ApproxFisher"& multigroup) {
    pars1 = pars[[1]]
    pars2 = pars[[2]]
    df1 = mirt::simdata(a = pars1$a, d = pars1$d, N = approx.npers%/%2, itemtype = itemtype)
    df2 = mirt::simdata(a = pars2$a, d = pars2$d, N = approx.npers%/%2, itemtype = itemtype)
    group= rep(c("A","B"),each=nrow(df1))
    dat =rbind(df1,df2)

    synt = mirt::multipleGroup(dat,group = group,itemtype = itemtype,model=model,technical = list(NCYCLES = 1),pars="values")
    apars = which(synt$name=="a1")
    dpars = which(synt$name=="d")
    synt$lbound[apars] =  synt$ubound[apars] = c(pars1$a,pars2$a)
    synt$lbound[dpars] = synt$ubound[dpars] = c(pars1$d,pars2$d)
    mml = mirt::multipleGroup(dat,model=model,itemtype = itemtype,SE = TRUE,SE.type = SE.type,technical = list(NCYCLES = NCYCLES),verbose = FALSE,pars=synt,group=group)
  }
    re = solve(mirt::vcov(mml))/nrow(dat)

  return(re)
}
