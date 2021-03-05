

#' Calculate information matrix (expected / observed)
#'
#' @param pars Parameter Set
#' @param method string, method to be used, e.g. "Fisher" for the Fisher-expected information matrix calculated using the mirt package
#' @param approx.npers integer, number of persons in the sampling based approach
#' @param data dataset for observed matrix methods
#'
#' @return
#' @export
#'
#' @examples
infmat = function(pars,method = "Fisher",approx.npers = 10^6,data=NULL,multigroup = FALSE,model = 1,itemtype="2PL") {

  if(is.null(pars$g)) {pars$g = rep(0,length(pars$a))}

  # j=0
  # success=FALSE
  # while((!success) & (j<10)) {
  #   j=j+1
  #   if(j>1) {print(paste("retry #",j-1))}

  if (method=="Fisher" & !multigroup) {
    dat = mirt::simdata(a = pars$a,d = pars$d,guess=pars$g,N =1000,itemtype = itemtype)
    synt = mirt::mirt(dat,1,itemtype = itemtype,technical = list(NCYCLES = 1),pars="values")
    apars = which(synt$name=="a1")
    synt$lbound[apars] =  synt$ubound[apars] = pars$a
    dpars = which(synt$name=="d")
    synt$lbound[dpars] = synt$ubound[dpars] = pars$d
    gpars = which(synt$name=="g")
    synt$lbound[gpars] = synt$ubound[gpars] = pars$g
    mml = mirt::mirt(dat,model=model,itemtype = itemtype,SE = TRUE,SE.type = "Fisher",technical = list(NCYCLES = 5000),verbose = FALSE,pars=synt)

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
    mml = mirt::multipleGroup(dat,model=model,itemtype = itemtype,SE = TRUE,SE.type = "Fisher",technical = list(NCYCLES = 5000),verbose = FALSE,pars=synt,group=group)
  }

  if (method=="ApproxFisher"& !multigroup) {

    dat = mirt::simdata(a = pars$a,d = pars$d,guess=pars$g,N =approx.npers,itemtype = itemtype)
    synt = mirt::mirt(dat,1,itemtype = itemtype,technical = list(NCYCLES = 1),pars="values")
    apars = which(synt$name=="a1")
    synt$lbound[apars] =  synt$ubound[apars] = pars$a
    dpars = which(synt$name=="d")
    synt$lbound[dpars] = synt$ubound[dpars] = pars$d
    gpars = which(synt$name=="g")
    synt$lbound[gpars] = synt$ubound[gpars] = pars$g
    mml = mirt::mirt(dat,1,itemtype = itemtype,SE = TRUE,SE.type = "Oakes",technical = list(NCYCLES = 5000),verbose = FALSE,pars=synt)

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
    mml = mirt::multipleGroup(dat,model=model,itemtype = itemtype,SE = TRUE,SE.type = "Oakes",technical = list(NCYCLES = 5000),verbose = FALSE,pars=synt,group=group)

  }

# delete/rework:

#   if (method=="customFisher") {
# # Fisher expected for the 3PL model.
#
#     patterns = as.matrix(expand.grid(lapply(1:length(pars$d),function(x) 0:1)))
#     load.functions(pars$itemtype)
#
#     res  = list()
#     for (i in 1:nrow(patterns)) {
#       l = ldot(patterns[i,],pars)
#       res[[i]] = l %*% t(l) * g(patterns[i,],pars)
#     }
#     re = Reduce('+', res)
#   }
#
#   if (method%in%c("mirtOakessim","mirtOakes")) {
#     if(method=="mirtOakessim"){df = mirt::simdata(a = pars$a,d = pars$d,guess=pars$g,N =simbased.npers,itemtype = itemtype)} else {df = data}
#
#     mml = mirt::mirt(df,1,itemtype = itemtype,SE = TRUE,SE.type = "Oakes",technical = list(NCYCLES = 5000),verbose = FALSE)
#     re = solve(mirt::vcov(mml))/nrow(df)
#   }
    re = solve(mirt::vcov(mml))/nrow(dat)
  #   success = !is.na(re[1,1])
  # }
  # if(j==2) {
  #   errorstring = paste(pars,method,multigroup,model)
  #   stop(errorstring)}

  return(re)
}
