
# Infmat ------------------------------------------------------------------

# Drop resmod argument, it's only needed for the number of items.


infmat = function(pars,method = "mirtFisher",resmod=NULL,simbased.npers = 200000,data=NULL) {
  if(is.null(pars$g)) {pars$g = rep(0,length(pars$a))}

  if (method=="mirtFishermultigroup") {

    if (!is.null(pars$a)){ # restricted model
      pars=list(list(pars,pars))}

    # Fisher Approach
    pars1 = pars[[1]][[1]]
    pars2 = pars[[1]][[2]]
    df1 = mirt::simdata(a = pars1$a, d = pars1$d, N = 5000, itemtype = pars1$itemtype)
    df2 = mirt::simdata(a = pars2$a, d = pars2$d, N = 5000, itemtype = pars2$itemtype)

    synt = mirt::mirt(df1,1,itemtype = c(pars$itemtype),technical = list(NCYCLES = 1),pars="values")
    apars = which(synt$name=="a1")
    dpars = which(synt$name=="d")
    synt$lbound[apars] =  synt$ubound[apars] = pars1$a
    synt$lbound[dpars] = synt$ubound[dpars] = pars1$d
    mml = mirt::mirt(df1,1,itemtype = c(pars$itemtype),SE = TRUE,SE.type = "Fisher",technical = list(NCYCLES = 5000),verbose = FALSE,pars=synt)
    re1 = solve(mirt::vcov(mml))/nrow(df1)

    synt$lbound[apars] =  synt$ubound[apars] = pars2$a
    synt$lbound[dpars] = synt$ubound[dpars] = pars2$d
    mml = mirt::mirt(df2,1,itemtype = c(pars$itemtype),SE = TRUE,SE.type = "Fisher",technical = list(NCYCLES = 5000),verbose = FALSE,pars=synt)
    re2 = solve(mirt::vcov(mml))/nrow(df2)

    re = matrix(0,nrow=nrow(re1)+2,ncol=ncol(re1)+2)
    re[1:nrow(re1),1:nrow(re1)]=re1
    re[(nrow(re1)+1):(nrow(re1)+2),(nrow(re1)+1):(nrow(re1)+2)]=re2[1:2,1:2]
    re[3:nrow(re1),(nrow(re1)+1):(nrow(re1)+2)]=re2[3:nrow(re2),1:2]
    re[(nrow(re1)+1):(nrow(re1)+2),3:nrow(re1)]=re2[1:2,3:nrow(re2)]

    re = re/2
    re[3:(nrow(re)-2),3:(nrow(re)-2)] = re[3:(nrow(re)-2),3:(nrow(re)-2)]*2

  }

  if (method%in%c("mirtOakesmultigroupsim","mirtOakesmultigroup")) {

    if (!is.null(pars$a)){ # restricted model
      pars=list(list(pars,pars))}

    # Oakes Approach
    pars1 = pars[[1]][[1]]
    pars2 = pars[[1]][[2]]

    if(method=="mirtOakesmultigroupsim"){
      df1 = mirt::simdata(a = pars1$a, d = pars1$d, N = simbased.npers/2, itemtype = pars1$itemtype)
      df2 = mirt::simdata(a = pars2$a, d = pars2$d, N = simbased.npers/2, itemtype = pars2$itemtype)
    } else {
      df1 = data[1:(nrow(data)/2),]
      df2 = data[(nrow(data)/2+1):nrow(data),]
    }

    mml = mirt::mirt(df1,1,itemtype = c(pars$itemtype),SE = TRUE,SE.type = "Oakes",technical = list(NCYCLES = 5000),verbose = FALSE)
    re1 = solve(mirt::vcov(mml))/nrow(df1)

    mml = mirt::mirt(df2,1,itemtype = c(pars$itemtype),SE = TRUE,SE.type = "Oakes",technical = list(NCYCLES = 5000),verbose = FALSE)
    re2 = solve(mirt::vcov(mml))/nrow(df2)

    re = matrix(0,nrow=nrow(re1)+2,ncol=ncol(re1)+2)
    re[1:nrow(re1),1:nrow(re1)]=re1
    re[(nrow(re1)+1):(nrow(re1)+2),(nrow(re1)+1):(nrow(re1)+2)]=re2[1:2,1:2]
    re[3:nrow(re1),(nrow(re1)+1):(nrow(re1)+2)]=re2[3:nrow(re2),1:2]
    re[(nrow(re1)+1):(nrow(re1)+2),3:nrow(re1)]=re2[1:2,3:nrow(re2)]

    re = re/2
    re[3:(nrow(re)-2),3:(nrow(re)-2)] = re[3:(nrow(re)-2),3:(nrow(re)-2)]*2

  }

  if (method=="mirtFisher") {
    df = mirt::simdata(a = pars$a,d = pars$d,guess=pars$g,N =1000,itemtype = pars$itemtype)
    synt = mirt::mirt(df,1,itemtype = c(pars$itemtype),technical = list(NCYCLES = 1),pars="values")
    apars = which(synt$name=="a1")
    synt$lbound[apars] =  synt$ubound[apars] = pars$a
    dpars = which(synt$name=="d")
    synt$lbound[dpars] = synt$ubound[dpars] = pars$d
    gpars = which(synt$name=="g")
    synt$lbound[gpars] = synt$ubound[gpars] = pars$g
    mml = mirt::mirt(df,1,itemtype = c(pars$itemtype),SE = TRUE,SE.type = "Fisher",technical = list(NCYCLES = 5000),verbose = FALSE,pars=synt)
    re = solve(mirt::vcov(mml))/nrow(df)

  }

  if (method=="customFisher") {

    patterns = as.matrix(expand.grid(lapply(1:resmod$n.items,function(x) 0:1)))
    load.functions(pars$itemtype)

    res  = list()
    for (i in 1:nrow(patterns)) {
      l = ldot(patterns[i,],pars)
      res[[i]] = l %*% t(l) * g(patterns[i,],pars)
    }
    re = Reduce('+', res)
  }


  if (method%in%c("mirtOakessim","mirtOakes")) {
    if(method=="mirtOakessim"){df = mirt::simdata(a = pars$a,d = pars$d,guess=pars$g,N =simbased.npers,itemtype = pars$itemtype)} else {df = data}
    mml = mirt::mirt(df,1,itemtype = c(pars$itemtype),SE = TRUE,SE.type = "Oakes",technical = list(NCYCLES = 5000),verbose = FALSE)
    re = solve(mirt::vcov(mml))/nrow(df)
  }

  if (method=="Thissen") {
    load.functions(pars$itemtype)
    re =  tmat(pars)
  }
  return(re)
}
