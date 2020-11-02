
library(spatstat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(mirt)
library(gridExtra)
library(simpackage)
library(parallel)
library(xtable)
library(Deriv)
library(digest)
library(ggfortify)
library(tidyr)
library(MASS)
library(grid)
library(stringr)
library(sn)
library(e1071)


# helper ------------------------------------------------------------------

coef.long = function(mirtfit=NULL,itemtype=NULL,shortcoef=NULL) {
  # extracts coefs from mml coef output
  if (is.null(shortcoef)) {
    if (itemtype=="2PL"& mirtfit@Call[[1]]=="multipleGroup") {
      pars = lapply(coef(mirtfit,simplify=TRUE),function(x) x$items[,1:2] %>% t() %>% as.numeric())
      re = c(pars$A,pars$B[1:2])
    } else if (itemtype=="2PL") {
      re = coef(mirtfit,simplify=TRUE)$items[,1:2] %>% t() %>% as.numeric()
    } else if (itemtype=="3PL") {
      re = coef(mirtfit,simplify=TRUE)$items[,1:3] %>% t() %>% as.numeric()
    } else if (itemtype=="gpcm") {
      nkatx = max(mirtfit@Data$data)
      a1 = coef(mirtfit,simplify=TRUE)$items
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


coef.short = function(mirtfit,itemtype) {
  # extracts coefs from mml coef output

  a = as.data.frame(coef(mirtfit,simplify=TRUE)$items)

  re =  list(a = a$a1, d  = a$d,g = a$g,itemtype=itemtype)
  if (itemtype=="gpcm") {
    re =  list(a = a$a1, d=cbind(rep(0,length(a$a1)),a$d1,a$d2),g = rep(0,length(a$a1)),itemtype=itemtype)
  }

  return(re)
}


load.functions = function (model) {
  # Loads model formula

  if (model=="2PL") {
    source('mmlpwr_func_2PL.R')
  }

  if (model=="2PLa") {
    source('mmlpwr_func_2PLa.R')
  }

  if (model=="3PL") {
    source('mmlpwr_func_3PL.R')
  }

  if (model=="3PL_nologit") {
    source('mmlpwr_func_3PL_nologit.R')
  }
  if (model=="3PL_logit") {
    source('mmlpwr_func_3PL_logit.R')
  }

  if (model=="gpcm") {
    source('mmlpwr_func_GPCM.R')
  }

}


maximizeL = function(null.hypothesis,alternative.hypothesis,bootstrap.start=TRUE) {
  #outputs pars that follow the restriction while exhibiting minimal KL distance
  print("maximizing L")
  if (null.hypothesis$preset=="1PLvs2PL") {

    pars = alternative.hypothesis
    load.functions(pars$itemtype)

    if(isTRUE(bootstrap.start)) {
      df = simdata(a = pars$a,d = pars$d,N =200000,itemtype = "2PL")
      mml = mirt(df,model = alternative.hypothesis$model,itemtype = null.hypothesis$itemtype,technical = list(NCYCLES = 5000),verbose = FALSE)
      startval=c(coef.short(mml,"2PL")$a[1],coef.short(mml,"2PL")$d)
    }else{
      startval = c(mean(pars$a),as.numeric(pars$d))
    }

    maxl2pre = maxl2preload(pars)
    optpar = optim(startval,function(x) {maxl2(x,pars,maxl2pre)})
    re = pars
    re$a = rep(optpar$par[1],length(pars$a))
    re$d = optpar$par[2:length(optpar$par)]
  }

  if (null.hypothesis$preset=="PCMvsGPCM") {

    pars = alternative.hypothesis
    load.functions(pars$itemtype)
    startval = c(mean(pars$a),as.numeric(pars$d[(null.hypothesis$n.items+1):length(as.numeric(pars$d))]))
    maxl2pre = maxl2preload(pars)
    optpar = optim(startval,function(x) {maxl2(x,pars,maxl2pre)})
    re = pars
    re$a = rep(optpar$par[1],length(pars$a))
    re$d = matrix(c(rep(0,null.hypothesis$n.items),optpar$par[2:length(optpar$par)]),ncol=ncol(pars$d))

  }

  if (null.hypothesis$preset=="3PLspec") {
    pars = alternative.hypothesis
    load.functions(pars$itemtype)
    re = pars
    re$g=rep(0.2,null.hypothesis$n.items)

    for (i in 1:length(pars$a)) {
      startval = c(re$a[i],re$d[i])
      optpar = optim(startval,function(x) {maxl(x,pars,i,.2)})
      re$a[i] = optpar$par[1]
      re$d[i] = optpar$par[2]
    }
  }

  if (null.hypothesis$preset=="DIF2PL") {

    pars1 = alternative.hypothesis[[1]][[1]]
    pars2 = alternative.hypothesis[[1]][[2]]

    load.functions(pars1$itemtype)
    re = pars1

    for (i in 1:length(pars1$a)) {
      startval = c(re$a[i],re$d[i])
      optpar = optim(startval,function(x) {maxl(x,pars1,pars2,i)})
      re$a[i] = optpar$par[1]
      re$d[i] = optpar$par[2]
    }

  }
  print("finished maximizing L")

  return(re)
}


get_ncp = function(chii,df) {
  ncp_x = mean(chii)-df
  chii[chii<=0] = .00000001
  if (ncp_x<0) {return(list(ncp=0,sd = NA))}
  chi_ncp <- fitdistr(chii,"chi-squared",start=list(ncp=0),
                      method="Brent",df=df,lower=0,upper=10000)
  return(list(ncp=chi_ncp$estimate,sd = chi_ncp$sd))
}


ssize = function(df,ncp,alpha,power) {
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
}

power = function(df,ncp,alpha,ssize) {
  crit = qchisq(1-alpha,df = df, ncp = 0) %>% as.numeric()
  1 - pchisq(q = crit, df = df, ncp = ncp*ssize)
}



# setup -------------------------------------------------------------------

setup.null.hypothesis = function(itemtype=NULL,n.items,model = 1,Amat=NULL, cvec=NULL, nullpars =NULL,preset=NULL) {
  #itemtype is the type of the unrestricted model
  # model is a mirt model for the unrestricted model (default 1, can also be multigroup)
  # Amat is the A matrix for the null hypothesis
  # cvec is the c vector for the null hypothesis
  # nullpars are the parameters under the null hypothesis
  # preset are presets

  if (preset == "1PLvs2PL") {
    itemtype = "2PL"
    Amat = c(1,0,-1,0,rep(0,(n.items-1)*2)) %>% rep(.,n.items-2) %>% c(.,1,0,-1,0) %>% matrix(.,ncol=n.items*2,byrow=TRUE)
    cvec = 0
    model = 1
    multigroup = FALSE
  }

  if (preset == "PCMvsGPCM") {
    itemtype = "gpcm"
    nkat = 3
    Amat = c(1,rep(0,nkat-1),-1,rep(0,nkat-1),rep(0,(n.items-1)*nkat)) %>% rep(.,n.items-2) %>% c(1,rep(0,nkat-1),-1,rep(0,nkat-1)) %>% matrix(.,ncol=n.items*nkat,byrow=TRUE)
    cvec = 0
    model = 1
    multigroup = FALSE
  }

  if (preset == "3PLspec") {
    itemtype = "3PL"
    Amat = c(0,0,1,rep(0,(n.items)*3)) %>% rep(.,n.items-1) %>% c(.,0,0,1) %>% matrix(.,ncol=n.items*3,byrow=TRUE)
    cvec = rep(.2,n.items)
    model = 1
    multigroup = FALSE
  }

  if (preset == "DIF2PL") {
    itemtype = "2PL"
    Amat = c(1,0,rep(0,(n.items-1)*2),-1,0) %>% c(.,c(0,1,rep(0,(n.items-1)*2),0,-1)) %>% matrix(.,ncol=(n.items+1)*2,byrow=TRUE)
    cvec = 0
    model = mirt.model(paste('F = 1-',n.items,'
                      CONSTRAINB = (2-',n.items,', d), (2-',n.items,', a1)'))
    multigroup = TRUE
  }

  re = list(itemtype=itemtype,
            n.items = n.items,
            model = model,
            Amat = Amat,
            cvec = cvec,
            preset = preset,
            multigroup = multigroup,
            df = nrow(Amat))

  return(re)
}


setup.alternative.hypothesis = function(null.hypothesis = NULL,esize= NULL,altpars= NULL,use.altpars = FALSE) {
  # altpars are the parameters under the alterNULLtive hypothesis

  n.items = null.hypothesis$n.items
  set.seed(123)

  if (is.null(esize)&!is.null(altpars)) {
    use.altpars = TRUE
    esize = "no" # setting esize temporarily
  }

  if (null.hypothesis$preset=="1PLvs2PL") {

  quantiles = (1:n.items)/(n.items+1)
  perm = sample(n.items)

  if (esize=="large") {
    apar = qlnorm(quantiles,sdlog = .1)
  } else if (esize=="small") {
    apar = qlnorm(quantiles,sdlog = .05)
  } else if (esize=="xsmall") {
    apar = qlnorm(quantiles,sdlog = .03)
  } else if (esize=="no") {
    apar = rep(1,n.items)}

  re = list(
    itemtype  = "2PL",
    a = apar,
    d = qnorm(quantiles)[perm],
    model = mirt.model(paste('F = 1-',n.items,'
                       CONSTRAIN = (1-',n.items,', a1)'))
  )
  re$longcoef = coef.long(shortcoef = re,itemtype=re$itemtype)

  } else if (null.hypothesis$preset=="PCMvsGPCM") {
    nkat= 3
    quantiles = (1:n.items)/(n.items+1)
    quantiles2 = (1:(n.items*(nkat-1)))/((n.items*(nkat-1))+1)
    perm = sample(n.items)

    if (esize=="large") {
      apar = qlnorm(quantiles,sdlog = .1)
    } else if (esize=="small") {
      apar = qlnorm(quantiles,sdlog = .05)
    } else if (esize=="no") {
      apar = rep(1,n.items)}

    re = list(
      itemtype = "gpcm",
      a = apar,
      d = cbind(rep(0,n.items),matrix(qnorm(quantiles2),ncol=(nkat-1))[perm,]),
      model = mirt.model(paste('F = 1-',n.items,'
                       CONSTRAIN = (1-',n.items,', a1)'))
    )
    re$longcoef = coef.long(shortcoef = re,itemtype=re$itemtype)

  } else if (null.hypothesis$preset=="3PLspec") {
    quantiles = (1:n.items)/(n.items+1)
    perm = sample(n.items)
    perm2 = sample(n.items)
    if (esize=="large") {mult = .4} else if (esize=="small") {mult=.2} else if (esize=="no") {mult=0}

    re = list(
      itemtype = "3PL",
      a = qlnorm(quantiles,sdlog = .05),
      d = qnorm(quantiles)[perm],
      g = qlnorm(quantiles,sdlog=.9)[perm2]^mult/5,
      model = mirt.model(paste('F = 1-',n.items,'
                       START = (1-',n.items,', g,.2),
                       FIXED = (1-',n.items,', g)')))
    re$longcoef = coef.long(shortcoef = re,itemtype=re$itemtype)

  } else if (null.hypothesis$preset=="DIF2PL") {
    quantiles = (1:n.items)/(n.items+1)
    perm = sample(n.items)
    if (esize=="large") {d = .4} else if (esize=="small") {d = .2} else if (esize=="no") {d=0}

    reA = reB = list(
      itemtype  = "2PL",
      a = qlnorm(quantiles,sdlog = .1),
      d = qnorm(quantiles)[perm]
    )
    reA$a[1] = 1+d/2
    reA$d[1] = d/2
    reB$a[1] = 1-d/2
    reB$d[1] = -d/2
    reA$longcoef = coef.long(shortcoef = reA,itemtype=reA$itemtype)
    reB$longcoef = coef.long(shortcoef = reB,itemtype=reB$itemtype)
    re = list(parsets = list(reA,reB), model = 1)
    re$longcoef = c(reA$longcoef,reB$longcoef[1:2])
  }

  if (isTRUE(use.altpars)) {
    re$a = altpars$a
    re$d = altpars$d
    re$g = altpars$g
    re$longcoef = coef.long(shortcoef = re,itemtype=re$itemtype)
  }
    return(re)
}



# simtools ----------------------------------------------------------------

setup.data = function(null.hypothesis, alternative.hypothesis, n,dist.type="norm") {
  if (dist.type=="norm") {
    distfun = function(x) {rnorm(x) %>% matrix(.,ncol=1) }
  } else if (dist.type =="unif"){
    distfun = function(x) {runif(x,min=-2,max=2) %>% matrix(.,ncol=1)}
  } else if (dist.type =="skewed"){
    distfun = function(x) {rsn(n=x, xi=0, omega=1, alpha=4)  %>% matrix(.,ncol=1)}

    # Skewed Dist Properties
    nmbrs = rsn(n=100000, xi=0, omega=1, alpha=4)
    nmbrs %>% summary()
    nmbrs %>% sd()
    nmbrs %>% skewness()


  } else {
    stop('dist.type unknown')
    }
  # distfun(10000) %>% as.numeric() %>% hist()
  # distfun(10000) %>% as.numeric() %>% mean()

  group = rep("A",n)

    if (isTRUE(null.hypothesis$multigroup)) {
      df = lapply(alternative.hypothesis$parsets,function(pars) simdata(a = pars$a,d = pars$d,Theta = distfun(n/2),itemtype = null.hypothesis$itemtype)) %>% do.call(rbind,.)
      group=rep(c("A","B"),each=n/2)
    } else if (null.hypothesis$itemtype %in% c("2PL","gpcm")) {
      df = simdata(a = alternative.hypothesis$a,d = alternative.hypothesis$d,Theta = distfun(n),itemtype = null.hypothesis$itemtype)
    } else if (null.hypothesis$itemtype=="3PL") {
      df = simdata(a = alternative.hypothesis$a,d = alternative.hypothesis$d,guess=alternative.hypothesis$g,Theta = distfun(n),itemtype = null.hypothesis$itemtype)
    }

  re=list(data = df,group=group)
  return(re)
}

sim.stats = function(null.hypothesis, alternative.hypothesis,runs, n.pers,include.score=TRUE,dist.type="norm") {

  res = list()
  for (i in 1:runs) {
    paste("started run",i,"of",runs) %>% print()
    datasetx = setup.data(null.hypothesis = null.hypothesis,alternative.hypothesis = alternative.hypothesis,n = n.pers,dist.type=dist.type)
    fittedx = mirt.fit(null.hypothesis = null.hypothesis,alternative.hypothesis = alternative.hypothesis, dataset = datasetx)
    res[[i]] = c(wald(fittedx,null.hypothesis),lr(fittedx))
    if (isTRUE(include.score)) {res[[i]] = c(res[[i]],score(fittedx,null.hypothesis))}
  }
  re = do.call(rbind,res) %>%  as.data.frame()
  names(re) =  c("Wald","LR","Score")[c(TRUE,TRUE,include.score)]

  return(re)
}

stats.simfunc = function(sim,dist.type="norm") {
  library(mirt)
  library(dplyr)
  library(spatstat)
  library(Deriv)
  library(digest)
  library(Matrix)
  library(MASS)
  library(sn)

  stats = sim.stats(sim$null.hypothesis, sim$alternative.hypothesis,runs=500, n.pers=sim$n.pers,dist.type=dist.type)
  return(list(sim=sim,stats=stats))
}

ncp.simfunc = function(sim) {

  library(mirt)
  library(dplyr)
  library(spatstat)
  library(Deriv)
  library(digest)
  library(Matrix)
  library(MASS)
  library(sn)

  if (sim$n.items<=15) {
    analytical = c(
      ncp.wald(sim$null.hypothesis,sim$alternative.hypothesis),
      ncp.lr(sim$null.hypothesis,sim$alternative.hypothesis),
      ncp.score(sim$null.hypothesis,sim$alternative.hypothesis)
    )
  } else {analytical = NULL}
  simbased  = ncp.sim(sim$null.hypothesis,sim$alternative.hypothesis)
  simbased.wald = ncp.wald(sim$null.hypothesis,sim$alternative.hypothesis ,method="mirtOakessim")

  return(list(sim=sim,analytical=analytical,simbased=simbased,simbased.wald=simbased.wald))
}

# Statistics --------------------------------------------------------------

mirt.fit = function(null.hypothesis,alternative.hypothesis,dataset=NULL,infmat.method = "Oakes") {

  if(is.null(dataset$data)) {
    temp= dataset
    dataset=list()
    dataset$data = temp}

  if (isFALSE(null.hypothesis$multigroup)) {
  unres =  mirt(dataset$data,model = null.hypothesis$model,itemtype = null.hypothesis$itemtype,SE = TRUE,SE.type = infmat.method,technical = list(NCYCLES = 5000),verbose = FALSE)
  } else {
  unres = multipleGroup(dataset$data,model = null.hypothesis$model,itemtype = null.hypothesis$itemtype,SE = TRUE,SE.type = infmat.method,technical = list(NCYCLES = 5000),verbose = FALSE,group=dataset$group)
  }
  res = mirt(dataset$data,model = alternative.hypothesis$model,itemtype = null.hypothesis$itemtype,technical = list(NCYCLES = 5000),verbose = FALSE)
  return(list(unres=unres,res=res))
}

wald = function(fitted,null.hypothesis) {
  pars=coef.long(mirtfit = fitted$unres,itemtype = null.hypothesis$itemtype)
  A = null.hypothesis$Amat
  dif = A%*%pars-null.hypothesis$cvec
  sigma=vcov(fitted$unres)
  re = t(dif) %*% solve(A%*% sigma %*% t(A)) %*% dif %>% as.numeric()
  return(re)
}

lr = function(fitted) {
  re = 2*(logLik(fitted$unres)-logLik(fitted$res))
  return(re)
}

score = function(fitted,null.hypothesis) {
  parsr = coef.short(fitted$res,itemtype = null.hypothesis$itemtype)
  load.functions(null.hypothesis$itemtype)
  patterns = fitted$res@Data$data
  hashs = apply(patterns,1,digest) %>% factor(.,levels=unique(.))
  rownames(patterns)=hashs
  up = unique(patterns)
  freq = table(hashs)

  if (isFALSE(null.hypothesis$multigroup)) {

    ly = lapply(rownames(up),function(i) ldot(up[i,],parsr)*freq[i] )
    lx = do.call(rbind,ly) %>% colSums() %>% array(.,dim=c(length(.),1))

    if (null.hypothesis$itemtype == "3PL") {method="customFisher"} else {method = "mirtOakes"}

  } else {

    freq1 = table(hashs[1:(nrow(patterns)/2)])
    freq2 = table(hashs[(nrow(patterns)/2+1):nrow(patterns)])

    ly = lapply(rownames(up),function(i) {l = ldot(up[i,],parsr); list(freq[i] * l,freq1[i] * l,freq2[i] * l)} )
    lx = lapply(ly,function(x) x[[1]]) %>% do.call(rbind,.) %>% colSums() %>% array(.,dim=c(length(.),1))
    lx1 = lapply(ly,function(x) x[[2]]) %>% do.call(rbind,.) %>% colSums() %>% array(.,dim=c(length(.),1))
    lx2 = lapply(ly,function(x) x[[3]]) %>% do.call(rbind,.) %>% colSums() %>% array(.,dim=c(length(.),1))

    lx = c(lx1[1:2,],lx[3:length(lx),],lx2[1:2,]) %>% array(.,dim=c(length(.),1))

    # method = "mirtFishermultigroup"
    method = "mirtOakesmultigroup"
  }
#     sigma = (infmat(parsr,method=method,null.hypothesis = null.hypothesis)*nrow(patterns)) %>% solve()

    sigma = (infmat(parsr,method=method,null.hypothesis = null.hypothesis,data=fitted$res@Data$data)*nrow(patterns)) %>% solve()


  re = t(lx) %*% sigma %*% lx %>% as.numeric()
  return(re)
}



# Infmat ------------------------------------------------------------------

infmat = function(pars,method = "mirtFisher",null.hypothesis=NULL,n.pers = 100000,data=NULL) {
  if(is.null(pars$g)) {pars$g = rep(0,length(pars$a))}

  if (method=="mirtFishermultigroup") {

      if (!is.null(pars$a)){ # restricted model
        pars=list(list(pars,pars))}

      # Fisher Approach
      pars1 = pars[[1]][[1]]
      pars2 = pars[[1]][[2]]
      n.pers = 5000
      df1 = simdata(a = pars1$a, d = pars1$d, N = n.pers, itemtype = pars1$itemtype)
      df2 = simdata(a = pars2$a, d = pars2$d, N = n.pers, itemtype = pars2$itemtype)

      synt = mirt(df1,1,itemtype = c(pars$itemtype),technical = list(NCYCLES = 1),pars="values")
      apars = which(synt$name=="a1")
      dpars = which(synt$name=="d")
      synt$lbound[apars] =  synt$ubound[apars] = pars1$a
      synt$lbound[dpars] = synt$ubound[dpars] = pars1$d
      mml = mirt(df1,1,itemtype = c(pars$itemtype),SE = TRUE,SE.type = "Fisher",technical = list(NCYCLES = 1000),verbose = FALSE,pars=synt)
      re1 = solve(vcov(mml))/nrow(df1)

      synt$lbound[apars] =  synt$ubound[apars] = pars2$a
      synt$lbound[dpars] = synt$ubound[dpars] = pars2$d
      mml = mirt(df2,1,itemtype = c(pars$itemtype),SE = TRUE,SE.type = "Fisher",technical = list(NCYCLES = 1000),verbose = FALSE,pars=synt)
      re2 = solve(vcov(mml))/nrow(df2)

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
      df1 = simdata(a = pars1$a, d = pars1$d, N = n.pers, itemtype = pars1$itemtype)
      df2 = simdata(a = pars2$a, d = pars2$d, N = n.pers, itemtype = pars2$itemtype)
    } else {
      df1 = data[1:(nrow(data)/2),]
      df2 = data[(nrow(data)/2+1):nrow(data),]
    }

    mml = mirt(df1,1,itemtype = c(pars$itemtype),SE = TRUE,SE.type = "Oakes",technical = list(NCYCLES = 1000),verbose = FALSE)
    re1 = solve(vcov(mml))/nrow(df1)

    mml = mirt(df2,1,itemtype = c(pars$itemtype),SE = TRUE,SE.type = "Oakes",technical = list(NCYCLES = 1000),verbose = FALSE)
    re2 = solve(vcov(mml))/nrow(df2)

    re = matrix(0,nrow=nrow(re1)+2,ncol=ncol(re1)+2)
    re[1:nrow(re1),1:nrow(re1)]=re1
    re[(nrow(re1)+1):(nrow(re1)+2),(nrow(re1)+1):(nrow(re1)+2)]=re2[1:2,1:2]
    re[3:nrow(re1),(nrow(re1)+1):(nrow(re1)+2)]=re2[3:nrow(re2),1:2]
    re[(nrow(re1)+1):(nrow(re1)+2),3:nrow(re1)]=re2[1:2,3:nrow(re2)]

    re = re/2
    re[3:(nrow(re)-2),3:(nrow(re)-2)] = re[3:(nrow(re)-2),3:(nrow(re)-2)]*2

  }

  if (method=="mirtFisher") {

    df = simdata(a = pars$a,d = pars$d,guess=pars$g,N =1000,itemtype = pars$itemtype)
    synt = mirt(df,1,itemtype = c(pars$itemtype),technical = list(NCYCLES = 1),pars="values")
    apars = which(synt$name=="a1")
    synt$lbound[apars] =  synt$ubound[apars] = pars$a
    dpars = which(synt$name=="d")
    synt$lbound[dpars] = synt$ubound[dpars] = pars$d
    gpars = which(synt$name=="g")
    synt$lbound[gpars] = synt$ubound[gpars] = pars$g
    mml = mirt(df,1,itemtype = c(pars$itemtype),SE = TRUE,SE.type = "Fisher",technical = list(NCYCLES = 1000),verbose = FALSE,pars=synt)
    re = solve(vcov(mml))/nrow(df)

  }

  if (method=="customFisher") {

      patterns = as.matrix(expand.grid(lapply(1:null.hypothesis$n.items,function(x) 0:1)))
      load.functions(pars$itemtype)

      res  = list()
      for (i in 1:nrow(patterns)) {
        l = ldot(patterns[i,],pars)
        res[[i]] = l %*% t(l) * g(patterns[i,],pars)
      }
      re = Reduce('+', res)
  }


  if (method%in%c("mirtOakessim","mirtOakes")) {
    if(method=="mirtOakessim"){df = simdata(a = pars$a,d = pars$d,guess=pars$g,N =n.pers,itemtype = pars$itemtype)} else {df = data}
    mml = mirt(df,1,itemtype = c(pars$itemtype),SE = TRUE,SE.type = "Oakes",technical = list(NCYCLES = 1000),verbose = FALSE)
    re = solve(vcov(mml))/nrow(df)
  }

  if (method=="Thissen") {
    load.functions(pars$itemtype)
    re =  tmat(pars)
  }

  return(re)
}


# Noncentrality -----------------------------------------------------------


ncp.sim = function(null.hypothesis, alternative.hypothesis, n=1,n.pers=15000,runs=3,include.score=TRUE) {

  if (null.hypothesis$n.items<=10) {
    n.pers=100000
    # n.pers=10000

    runs=3
  } else {
    n.pers=10000
    # n.pers=10000
    runs=3
  }

  res = list()
  for (i in 1:runs) {
    paste("run",i,"of",runs,"runs") %>% print()
    datasetx = setup.data(null.hypothesis = null.hypothesis,alternative.hypothesis = alternative.hypothesis,n = n.pers)
    fittedx = mirt.fit(null.hypothesis = null.hypothesis,alternative.hypothesis = alternative.hypothesis, dataset = datasetx)
    print("fitted, now calculating stats")
  res[[i]] = c(wald(fittedx,null.hypothesis),lr(fittedx))
  if (isTRUE(include.score)) {res[[i]] = c(res[[i]],score(fittedx,null.hypothesis))}
  }
  re = do.call(rbind,res) %>% apply(.,2,function(x) get_ncp(x,df=null.hypothesis$df)$ncp) %>% c()
  re = re/n.pers
  return(re*n)
}


ncp.wald = function(null.hypothesis,alternative.hypothesis,method="mirtFisher",n=1) {
  A = null.hypothesis$Amat
  dif = A%*%alternative.hypothesis$longcoef-null.hypothesis$cvec
  if (isTRUE(null.hypothesis$multigroup) & method=="mirtFisher") {method="mirtFishermultigroup"}
  if (isTRUE(null.hypothesis$multigroup) & method=="mirtOakessim") {method="mirtOakesmultigroupsim"}
  print("calculating fisher expected infmat")
  sigma = infmat(alternative.hypothesis,method=method,null.hypothesis = null.hypothesis) %>% solve()
  print("finished calculating fisher expected infmat")

  re = t(dif) %*% solve(A%*% sigma %*% t(A)) %*% dif %>% c()
  return(re*n)
  }


ncp.lr = function(null.hypothesis,alternative.hypothesis,method="mirtFisher",n=1) {

  parsr = maximizeL(null.hypothesis,alternative.hypothesis)
  res = c()
  unres = c()
  load.functions(null.hypothesis$itemtype)

  if (isFALSE(null.hypothesis$multigroup)) {
    pars = alternative.hypothesis
    n.kat = max(ncol(pars$d),2)
    patterns = as.matrix(expand.grid(lapply(1:null.hypothesis$n.items,function(x) 0:(n.kat-1))))

    for (i in 1:nrow(patterns)) {
      gr = g(patterns[i,],parsr)
      gu = g(patterns[i,],pars)
      res[i] = log(gr) * (gu)
      unres[i] = log(gu) * (gu)
    }

  } else {

    pars1 = alternative.hypothesis[[1]][[1]]
    pars2 = alternative.hypothesis[[1]][[2]]
    n.kat = max(ncol(pars1$d),2)
    patterns = as.matrix(expand.grid(lapply(1:null.hypothesis$n.items,function(x) 0:(n.kat-1))))

    for (i in 1:nrow(patterns)) {
      gr = g(patterns[i,],parsr)
      g1 = g(patterns[i,],pars1)
      g2 = g(patterns[i,],pars2)
      res[i] = log(gr) * (.5 * g1 + .5 * g2)
      unres[i] = .5 * log(g1) * g1 + .5 * log(g2) * g2
    }
  }

  res = sum(res)
  unres = sum(unres)

  re = 2*(unres-res) %>% c()

  return(re*n)

}


ncp.score = function(null.hypothesis,alternative.hypothesis,method="mirtFisher",n=1) {

  parsr = maximizeL(null.hypothesis,alternative.hypothesis)
  load.functions(null.hypothesis$itemtype)

  if (isFALSE(null.hypothesis$multigroup)) {
    pars = alternative.hypothesis
    n.kat = max(ncol(pars$d),2)
    patterns = lapply(1:null.hypothesis$n.items,function(x) 0:(n.kat-1)) %>% expand.grid() %>% as.matrix() %>% split(.,seq(nrow(.)))
    i=0
    ly = lapply(patterns,function(x) {ldot(x,parsr) * g(x,pars)})
    lx = do.call(rbind,ly) %>% colSums() %>% array(.,dim=c(length(.),1))
    if(null.hypothesis$itemtype == "3PL") {method="customFisher"} else {method = "mirtFisher"}

  } else {

    pars1 = alternative.hypothesis[[1]][[1]]
    pars2 = alternative.hypothesis[[1]][[2]]
    n.kat = max(ncol(pars1$d),2)
    patterns = lapply(1:null.hypothesis$n.items,function(x) 0:(n.kat-1)) %>% expand.grid() %>% as.matrix() %>% split(.,seq(nrow(.)))

    freq12 = lapply(patterns,function(x) c(g(x,pars1)/2,g(x,pars2)/2)) %>% do.call(rbind,.)
    freq = cbind(rowSums(freq12),freq12)
    ly = lapply(1:length(patterns),function(i) {paste(100*i/length(patterns)%>%round(.,1),"% finished") %>% print()
      l = ldot(patterns[[i]],parsr); list(freq[i,1] * l,freq[i,2] * l,freq[i,3] * l)} )
    lx = lapply(ly,function(x) x[[1]]) %>% do.call(rbind,.) %>% colSums() %>% array(.,dim=c(length(.),1))
    lx1 = lapply(ly,function(x) x[[2]]) dplyr::%>% do.call(rbind,.) dplyr::%>% colSums() dplyr::%>% array(.,dim=c(length(.),1))
    lx2 = lapply(ly,function(x) x[[3]]) dplyr::%>% do.call(rbind,.) dplyr::%>% colSums() dplyr::%>% array(.,dim=c(length(.),1))

    lx = c(lx1[1:2,],lx[3:length(lx),],lx2[1:2,]) dplyr::%>% array(.,dim=c(length(.),1))

    method = "mirtFishermultigroup"
  }
  # sigma = infmat(alternative.hypothesis,method=method,null.hypothesis = null.hypothesis) dplyr::%>% solve()
  sigma = infmat(parsr,method=method,null.hypothesis = null.hypothesis) %>% solve()

  re = t(lx) %*% sigma %*% lx %>% c()

  return(re*n)
}






# Plots -------------------------------------------------------------------


qq.plot = function(res.stat,res.ncp,preset) {

  setup.data = function(x) {

    df = x$sim$null.hypothesis$df

    n.pers = x$sim$n.pers
    if(!is.null(x$ncp$analytical)) {
      ncp.w = (x$ncp$analytical*n.pers)[1]
      ncp.l = (x$ncp$analytical*n.pers)[2]
      ncp.s = (x$ncp$analytical*n.pers)[3]
    } else {
      # ncp.w = (x$ncp$simbased*n.pers)[1]
      ncp.w = (x$ncp$simbased.wald*n.pers)
      ncp.l = (x$ncp$simbased*n.pers)[2]
      ncp.s = (x$ncp$simbased*n.pers)[3]

    }

    exp.w = qchisq((1:length(x$stats$Wald))/(length(x$stats$Wald)+1) ,df = df, ncp = ncp.w) %>% as.numeric()
    exp.l = qchisq((1:length(x$stats$LR))/(length(x$stats$LR)+1) ,df = df, ncp = ncp.l) %>% as.numeric()
    exp.s = qchisq((1:length(x$stats$Score))/(length(x$stats$Score)+1) ,df = df, ncp = ncp.s) %>% as.numeric()

    obs.w = as.numeric(quantile(x$stats$Wald,(1:length(x$stats$Wald))/(length(x$stats$Wald)+1)))
    obs.l = as.numeric(quantile(x$stats$LR,(1:length(x$stats$LR))/(length(x$stats$LR)+1)))
    obs.s = as.numeric(quantile(x$stats$Score,(1:length(x$stats$Score))/(length(x$stats$Score)+1)))

    ks.w = ks.test(x$stats$Wald, "pchisq",df=df,ncp=ncp.w)$p.value
    ks.l = ks.test(x$stats$LR, "pchisq",df=df,ncp=ncp.l)$p.value
    ks.s = ks.test(x$stats$Score, "pchisq",df=df,ncp=ncp.s)$p.value

    df = data.frame(obs=c(obs.w,obs.l,obs.s),exp=c(exp.w,exp.l,exp.s),kind = factor(c(rep("Wald",length(obs.w)),rep("LR",length(obs.l)),rep("Score",length(obs.s))),levels=c("Wald","LR","Score")),n = n.pers,esize = x$sim$esize)
    df2 = data.frame(ks=c(ks.w,ks.l,ks.s),kind = factor(c("Wald","LR","Score")),n = n.pers,esize = x$sim$esize)
    return(list(df,df2))
 }

  for (i in 1:length(res.stat)) {
    res.stat[[i]]$ncp = res.ncp[[i]]
  }

  indi = sapply(res.stat,function(x) isTRUE(x$sim$preset==preset))
  resx = res.stat[indi]


  datax = lapply(resx,FUN=setup.data)
  data = do.call(rbind,lapply(datax,function(x) x[[1]]))
  dataks = do.call(rbind,lapply(datax,function(x) x[[2]]))
  dataks$ks=dataks$ks*nrow(dataks)
  dataks$sig =""

  dataks$sig[dataks$ks<.05] = paste(dataks$kind[dataks$ks<.05],"*",sep="")


  sequence = seq(1,nrow(dataks),3)
  dataks$sig[sequence] = do.call(rbind,lapply(sequence,function(x) paste(dataks$sig[x:(x+2)],collapse=" ")))
  dataks$sig[c(sequence+1,sequence+2)] = ""

  # New facet label names for dose variable
  n.labs <- paste("n =",unique(data$n))
  names(n.labs) <- as.character(unique(data$n))

  # New facet label names for supp variable
  esize.labs <- c("No effect", "Small effect","Large effect")
  names(esize.labs) <- c("no", "small","large")

  re = ggplot(data, aes(x=exp, y=obs,group=kind)) +
    geom_point(aes(color=kind),alpha = 1) +
    xlab("Expected Quantiles")+  ylab("Observed Quantiles") + geom_abline(linetype="dashed", color = "black") + theme(legend.position = "bottom") + facet_grid(vars(esize),vars(n),scale="free", labeller = labeller(n = n.labs, esize = esize.labs)) + labs(color = "Statistic") + scale_color_grey(start=.0,end=.6) + theme_bw() +
    geom_text(dataks,mapping=aes(x=Inf,y=-Inf,label=sig),hjust=1.2,vjust=-1 )


  return(re)
}



ks.table = function(res.stat,res.ncp,preset) {
  #COLS: Hypothesis, Effect Size, Sample Size, Wald(D,p), LR(D,p), Score(D,p)
  # *significant after bonferroni correction

  setup.data = function(x) {

    df = x$sim$null.hypothesis$df

    n.pers = x$sim$n.pers

    if(!is.null(x$ncp$analytical)) {
      ncp.w = (x$ncp$analytical*n.pers)[1]
      ncp.l = (x$ncp$analytical*n.pers)[2]
      ncp.s = (x$ncp$analytical*n.pers)[3]
    } else {
      # ncp.w = (x$ncp$simbased*n.pers)[1]
      ncp.w = (x$ncp$simbased.wald*n.pers)
      ncp.l = (x$ncp$simbased*n.pers)[2]
      ncp.s = (x$ncp$simbased*n.pers)[3]

    }

    ks.wx = ks.test(x$stats$Wald, "pchisq",df=df,ncp=ncp.w)
    ks.lx = ks.test(x$stats$LR, "pchisq",df=df,ncp=ncp.l)
    ks.sx = ks.test(x$stats$Score, "pchisq",df=df,ncp=ncp.s)


    ds = format(round(c(ks.wx$statistic,ks.lx$statistic,ks.sx$statistic),3),nsmall=3)
    ps = paste(format(round(c(ks.wx$p.value,ks.lx$p.value,ks.sx$p.value),3),nsmall=3),c(if(ks.wx$p.value*27<.05) {"*"}else{" "},if(ks.lx$p.value*27<.05) {"*"}else{" "},if(ks.sx$p.value*27<.05) {"*"}else{" "}),sep="")

    df = c(
      hyp = as.character(x$sim$null.hypothesis$preset),
      items = x$sim$n.items,
      esize = as.character(x$sim$esize),
      n.pers = x$sim$n.pers,
      d1 = ds[1],p1 = ps[1],d1 = ds[2],p1 = ps[2],d1 = ds[3],p1 = ps[3]
    )

    names(df) = c("Hypothesis", "Items","Effect size", "Sample size","D","p","D","p","D","p")
    return(df)
  }

  for (i in 1:length(res.stat)) {
    res.stat[[i]]$ncp = res.ncp[[i]]
  }

  # indi = sapply(res.stat,function(x) isTRUE(x$sim$preset==preset))
  # resx = res.stat[indi]
  resx = res.stat

  re = as.data.frame(do.call(rbind,lapply(resx,FUN=setup.data)))

  return(re)
}



power.plot = function(res.stat,res.ncp,preset,alpha=.05,analytical.only=FALSE) {

  pw.sd = function(ncp,df,n.runs,ci=.95,alpha=.05) {

    crit = qchisq(1-alpha,df = df, ncp = 0) %>% as.numeric()

    res0 = c()

    for (i in 1:1000) {
      res0[i]= mean(rchisq(n.runs,df=df,ncp=ncp)>crit)
    }

    re = quantile(res0,c((1-ci)/2,1-(1-ci)/2))

    return(re)
  }

  pw.sd2 = function(pw,n.runs,ci=.95) {

    res = sapply(1:10000,function(x) {(runif(n.runs) < pw) %>% mean() })

    re = quantile(res,c((1-ci)/2,1-(1-ci)/2))

    return(re)
  }

  setup.data = function(x,simbased.only) {
    df = x$sim$null.hypothesis$df

    runs = length(x$stats$Wald)
    n.pers = x$sim$n.pers

    crit = qchisq(1-alpha,df = df, ncp = 0) %>% as.numeric()
    pwro.w = mean(x$stats$Wald>crit)
    pwro.l = mean(x$stats$LR>crit)
    pwro.s = mean(x$stats$Score>crit)

    # ncp.w1 = (x$ncp$simbased*n.pers)[1]
    ncp.w1 = (x$ncp$simbased.wald*n.pers)
    ncp.l1 = (x$ncp$simbased*n.pers)[2]
    ncp.s1 = (x$ncp$simbased*n.pers)[3]

    pwre.w1 = 1 - pchisq(q = crit, df = df, ncp = ncp.w1)
    pwre.l1 = 1 - pchisq(q = crit, df = df, ncp = ncp.l1)
    pwre.s1 = 1 - pchisq(q = crit, df = df, ncp = ncp.s1)

    # pw.w.sd1 = pw.sd(ncp=ncp.w1,df=df,n.runs=runs)
    # pw.l.sd1 = pw.sd(ncp=ncp.l1,df=df,n.runs=runs)
    # pw.s.sd1 = pw.sd(ncp=ncp.s1,df=df,n.runs=runs)
    pw.w.sd1 = pw.sd2(pwro.w,n.runs=runs)
    pw.l.sd1 = pw.sd2(pwro.s,n.runs=runs)
    pw.s.sd1 = pw.sd2(pwro.l,n.runs=runs)

    if (!simbased.only) {
    ncp.w = (x$ncp$analytical*n.pers)[1]
    ncp.l = (x$ncp$analytical*n.pers)[2]
    ncp.s = (x$ncp$analytical*n.pers)[3]

    pwre.w = 1 - pchisq(q = crit, df = df, ncp = ncp.w)
    pwre.l = 1 - pchisq(q = crit, df = df, ncp = ncp.l)
    pwre.s = 1 - pchisq(q = crit, df = df, ncp = ncp.s)

    # pw.w.sd = pw.sd(ncp=ncp.w,df=df,n.runs=runs)
    # pw.l.sd = pw.sd(ncp=ncp.l,df=df,n.runs=runs)
    # pw.s.sd = pw.sd(ncp=ncp.s,df=df,n.runs=runs)

    pw.w.sd = pw.w.sd1
    pw.l.sd = pw.l.sd1
    pw.s.sd = pw.s.sd1

    re = data.frame(
      stat = factor(c("Wald","LR","Score"),levels = c("Wald","LR","Score")),
      pw.exp= c(pwre.w ,pwre.l,pwre.s),
      pw.exp1= c(pwre.w1 ,pwre.l1,pwre.s1),
      pw.obs= c(pwro.w ,pwro.l,pwro.s),
      # pw.lower= c(pw.sdx[1],pw.sdx[1],pw.sdx[1]),
      # pw.upper= c(pw.sdx[2],pw.sdx[2],pw.sdx[2]),
      pw.lower= c(pw.w.sd[1],pw.l.sd[1],pw.s.sd[1]),
      pw.upper= c(pw.w.sd[2],pw.l.sd[2],pw.s.sd[2]),
      # pw.lower1= c(pw.sdx1[1],pw.sdx1[1],pw.sdx1[1]),
      # pw.upper1= c(pw.sdx1[2],pw.sdx1[2],pw.sdx1[2]))
      pw.lower1= c(pw.w.sd1[1],pw.l.sd1[1],pw.s.sd1[1]),
      pw.upper1= c(pw.w.sd1[2],pw.l.sd1[2],pw.s.sd1[2]))
    }


    if(simbased.only) {
    re = data.frame(
      stat = factor(c("Wald","LR","Score"),levels = c("Wald","LR","Score")),
      pw.exp1= c(pwre.w1 ,pwre.l1,pwre.s1),
      pw.obs= c(pwro.w ,pwro.l,pwro.s),
      pw.lower1= c(pw.w.sd1[1],pw.l.sd1[1],pw.s.sd1[1]),
      pw.upper1= c(pw.w.sd1[2],pw.l.sd1[2],pw.s.sd1[2]))
    }
    return(re)
  }

  for (i in 1:length(res.stat)) {
    res.stat[[i]]$ncp = res.ncp[[i]]
  }
  simbased.only = is.null(res.ncp[[1]]$analytical)

  indi = sapply(res.stat,function(x) isTRUE(x$sim$preset==preset))
  resx = res.stat[indi]

  dat = lapply(resx,function(x) setup.data(x,simbased.only)) %>% do.call(rbind,.) %>% as.data.frame()

  re = ggplot(dat, aes(pw.exp1, pw.obs, colour = stat)) +
    geom_point(size=1) +
    geom_abline(linetype="dashed", color = "black") +
    geom_errorbar(data=dat, mapping= aes(ymin = pw.lower1, ymax = pw.upper1),width=.05) +
    scale_color_grey(start=.0,end=.6) + labs(colour = "Statistic") + theme_bw() + scale_y_continuous("Observed Power",breaks=seq(0,1,.2),limits = c(0, 1)) + scale_x_continuous("Expected Power",breaks=seq(0,1,.2),limits = c(0, 1))

if(!simbased.only) {

    re2 = ggplot(dat, aes(pw.exp, pw.obs, colour = stat)) +
    geom_point(size=1) +
    geom_abline(linetype="dashed", color = "black") +
    geom_errorbar(data=dat, mapping= aes(ymin = pw.lower, ymax = pw.upper),width=.05) +
    scale_color_grey(start=.0,end=.6)  + xlab("Expected Power") + labs(colour = "Statistic") + theme_bw() + scale_y_continuous("Observed Power",breaks=seq(0,1,.2),limits = c(0, 1)) + scale_x_continuous("Expected Power",breaks=seq(0,1,.2),limits = c(0, 1))

    re = grid.arrange(re2,re,ncol=2)

    if(isTRUE(analytical.only)) {re = re2}
}
  return(re)
}


curve.plot.multi = function(ncpslist,df,alpha=.05,nolegend=FALSE,scales="free") {
  #ncps is a list here


  ncps = ncpslist[[1]]

  upper =  uniroot(f = function(x) .999 - power(df=df,ncp=ncps[1],alpha=alpha,ssize=x),interval = c(0, 1000000),tol = .Machine$double.eps ^ 0.5)$root

  ssizes = seq(0,upper,10)
  res = list()
  for (i in 1:length(ssizes)) {
    res[[i]] = sapply(ncps,function (x) power(df=df,ncp=x,alpha=alpha,ssize=ssizes[i]))
  }

  dat = data.frame(do.call(rbind,res))

  names(dat) = c("Wald","LR","Score")
  dat$ssize = ssizes
  dat1 = as.data.frame(pivot_longer(dat,cols=names(dat)[1:3]))


  ncps = ncpslist[[2]]

  upper =  uniroot(f = function(x) .999 - power(df=df,ncp=ncps[1],alpha=alpha,ssize=x),interval = c(0, 1000000),tol = .Machine$double.eps ^ 0.5)$root

  ssizes = seq(0,upper,10)
  res = list()
  for (i in 1:length(ssizes)) {
    res[[i]] = sapply(ncps,function (x) power(df=df,ncp=x,alpha=alpha,ssize=ssizes[i]))
  }

  dat = data.frame(do.call(rbind,res))
  names(dat) = c("Wald","LR","Score")
  dat$ssize = ssizes
  dat2 = as.data.frame(pivot_longer(dat,cols=names(dat)[1:3]))

  dat1$set = "LSAT6"
  dat2$set = "LSAT7"
  dat = rbind(dat1,dat2)
  dat$name = factor(dat$name,levels=c("Wald","LR","Score"))

  re = ggplot(dat, aes(x=ssize, y=value, group=name)) +
    geom_line(aes(color=name),size=1) + scale_color_grey(start=.0,end=.6) + labs(title = "")+  xlab("Sample size")+ ylab("Power") +labs(color = "Statistic")  + facet_wrap(~set,scales = scales) + scale_y_continuous("Power",breaks=seq(0,1,.2)) + theme_bw()


  return(re)
}



chisq.dens = function(ncp,df,n,alpha=.05,lim=30) {

  crit= qchisq(1-alpha,df)

  p = ggplot(data.frame(x = c(0, lim)), aes(x)) +

    stat_function(
      fun = dchisq,
      geom = "area",
      # alpha = .5,
      args = list(
        df=df
      ),n = 5000,aes(fill = "Null")
    ) +
    stat_function(
      fun = dchisq,
      geom = "area",
      # alpha = .5,
      xlim = c(crit,lim),
      args = list(
        df=df,
        ncp=ncp[1]*n
      ),aes(fill = "Wald")
    ) +
    stat_function(
      fun = dchisq,
      geom = "area",
      # alpha = .1,
      xlim = c(crit,lim),
      args = list(
        df=df,
        ncp=ncp[2]*n
      ),aes(fill = "LR")
    ) +
    stat_function(
      fun = dchisq,
      geom = "area",
      # alpha = .1,
      xlim = c(crit,lim),
      args = list(
        df=df,
        ncp=ncp[3]*n
      ),aes(fill = "Score")
    ) +
    stat_function(
      fun = dchisq,
      geom = "line",
      args = list(
        df=df
      ),n = 5000
    )+
    stat_function(
      fun = dchisq,
      geom = "line",
      args = list(
        df=df,
        ncp=ncp[1]*n
      )
    ) +
    stat_function(
      fun = dchisq,
      geom = "line",
      args = list(
        df=df,
        ncp=ncp[2]*n
      )
    ) +
    stat_function(
      fun = dchisq,
      geom = "line",
      args = list(
        df=df,
        ncp=ncp[3]*n
      )
    ) +
    labs(
      title = "",
      x = "Statistic",
      y = "Density"
    ) + scale_fill_manual("",breaks=c("Null","Wald","LR","Score"),values = c("grey90","black","grey30","grey60")) + theme_bw()

  return(p)
}


