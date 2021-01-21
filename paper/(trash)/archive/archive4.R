#' Specify the restricted model
#'
#' Corresponds to the null hypothesis
#'
#'
#' @param n.items Number of items
#' @param preset Specifiy a preset, currently available: \itemize{
#' \item 1PLvs2PL
#' \item etc.}
#' @param itemtype Type of the unrestricted model (used in mirt later)
#' @param model Mirt Model used for estimation of !NUll or ALT model??
#' @param Amat Hypothesis Design Matrix A in \eqn{A\beta=c}
#' @param cvec Hypothesis Design vector c in\eqn{A\beta=c}
#' @param nullpars (optional) Parameters under the Null hypothesis
#'
#' @return List specifying necessary null hypothesis parameters
#' @export
#'
#' @examples setup.null.hypothesis(n.items=5,preset="1PLvs2PL")
#'
restricted_model.o = function(n.items,preset=NULL,itemtype=NULL,model = 1,Amat=NULL, cvec=NULL, nullpars =NULL) {

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

#' Specify the unrestricted model
#'
#' Corresponds to the alternative hypothesis
#'
#' @param null.hypothesis List as generated from \code{setup.null.hypothesis}
#' @param esize
#' @param altpars Parameters under the alternative hypothesis
#' @param use.altpars
#'
#' @return
#' @export
#'
#' @examples
unrestricted_model = function(null.hypothesis = NULL,esize= NULL,altpars= NULL,use.altpars = FALSE) {

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
      model = 1
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


#' Specify the restricted model
#'
#' Corresponds to the null hypothesis
#'
#'
#' @param n.items Number of items
#' @param preset Specifiy a preset, currently available: \itemize{
#' \item 1PLvs2PL
#' \item etc.}
#' @param itemtype Type of the unrestricted model (used in mirt later)
#' @param model Mirt Model used for estimation of !NUll or ALT model??
#' @param Amat Hypothesis Design Matrix A in \eqn{A\beta=c}
#' @param cvec Hypothesis Design vector c in\eqn{A\beta=c}
#' @param nullpars (optional) Parameters under the Null hypothesis
#'
#' @return List specifying necessary null hypothesis parameters
#' @export
#'
#' @examples setup.null.hypothesis(n.items=5,preset="1PLvs2PL")
#'
restricted_model = function(n.items,preset=NULL,itemtype=NULL,model = NULL,Amat=NULL, cvec=NULL, multigroup=NULL, nullpars =NULL) {

  if (!is.null(preset)) {
    re = preset$res(n.items=n.items)
  } else {
    re = list(
      itemtype=itemtype,
      n.items = n.items,
      model = model,
      Amat = Amat,
      cvec = cvec,
      multigroup = multigroup,
      nullpars = nullpars
    )
  }
  return(re)
}


#' Specify the unrestricted model
#'
#' Corresponds to the alternative hypothesis
#'
#' @param null.hypothesis List as generated from \code{setup.null.hypothesis}
#' @param esize
#' @param altpars Parameters under the alternative hypothesis
#' @param use.altpars
#'
#' @return
#' @export
#'
#' @examples
unrestricted_model = function(resmod,preset = NULL, esize= NULL,altpars= NULL,use.altpars=FALSE) {

  n.items = resmod$n.items

  if (!is.null(preset)) {
    re = preset$unres(n.items=n.items,esize=esize)
  } else {
    re = list(
      itemtype,
      n.items = n.items,
      model = model,
      Amat = Amat,
      cvec = cvec,
      preset = preset,
      multigroup = multigroup
    )
  }
  #
  #   if (is.null(esize)&!is.null(altpars)) {
  #     use.altpars = TRUE
  #     esize = "no" # setting esize temporarily
  #   }


  if (isTRUE(use.altpars)) {
    re$a = altpars$a
    re$d = altpars$d
    re$g = altpars$g
    re$longcoef = coef.long(shortcoef = re,itemtype=re$itemtype)
  }

  return(re)
}


maximizeL = function(hyp,bootstrap.start=TRUE) {
  #outputs pars that follow the restriction while exhibiting minimal KL distance

  resmod = hyp$resmod
  unresmod = hyp$unresmod

  # print("maximizing L")
  if (resmod$preset=="1PLvs2PL") {

    pars = unresmod
    load.functions(pars$itemtype)

    if(isTRUE(bootstrap.start)) {
      df = mirt::simdata(a = pars$a,d = pars$d,N =200000,itemtype = "2PL")
      mml = mirt(df,model = unresmod$model,itemtype = resmod$itemtype,technical = list(NCYCLES = 5000),verbose = FALSE)
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

  if (resmod$preset=="PCMvsGPCM") {

    pars = unresmod
    load.functions(pars$itemtype)
    startval = c(mean(pars$a),as.numeric(pars$d[(resmod$n.items+1):length(as.numeric(pars$d))]))
    maxl2pre = maxl2preload(pars)
    optpar = optim(startval,function(x) {maxl2(x,pars,maxl2pre)})
    re = pars
    re$a = rep(optpar$par[1],length(pars$a))
    re$d = matrix(c(rep(0,resmod$n.items),optpar$par[2:length(optpar$par)]),ncol=ncol(pars$d))

  }

  if (resmod$preset=="3PLspec") {
    pars = unresmod
    load.functions(pars$itemtype)
    re = pars
    re$g=rep(0.2,resmod$n.items)

    for (i in 1:length(pars$a)) {
      startval = c(re$a[i],re$d[i])
      optpar = optim(startval,function(x) {maxl(x,pars,i,.2)})
      re$a[i] = optpar$par[1]
      re$d[i] = optpar$par[2]
    }
  }

  if (resmod$preset=="DIF2PL") {

    pars1 = unresmod[[1]][[1]]
    pars2 = unresmod[[1]][[2]]

    load.functions(pars1$itemtype)
    re = pars1

    for (i in 1:length(pars1$a)) {
      startval = c(re$a[i],re$d[i])
      optpar = optim(startval,function(x) {maxl(x,pars1,pars2,i)})
      re$a[i] = optpar$par[1]
      re$d[i] = optpar$par[2]
    }

  }
  # print("finished maximizing L")

  return(re)
}





## Basic Usage with observed data

First, we load an observed dataset and estimate the 2PL parameters with mirt.

```{r}
dat <- expand.table(LSAT6)
pars = mirt(dat,1)
```

We setup the 1PL vs 2PL hypothesis and test it.
```{r}

hyp = setup_hypothesis(type = h_1PLvs2PL, altpars = pars)

fitted = mml.fit(hyp = hyp, data = dat)

# m0 = restricted_model(n.items=5,preset=preset_1PLvs2PL)
# m1 = unrestricted_model(resmod = m0,altpars=pars6)
# fitted6 = mirt.fit(rm = m0, urm = m1, data = dat6)

stat6 = c(wald(fitted6,null.hyp),lr(fitted6),score(fitted6,null.hyp))

# Test results (p-values)
pvals6 = stat6 %>% pchisq(.,df=null.hyp$df,ncp=0,lower.tail=FALSE)


```


Now, we estimate the power of this test.

```{r}
ncps = c(ncp.wald(hyp = hyp),ncp.lr(hyp = hyp),ncp.score(hyp = hyp))
power(hyp=hyp,ncp=ncps,alpha=.05,ssize=500)
ssize(hyp=hyp,ncp=ncps,alpha=.05,power=.80)

#Power given parameters are true
ncp6 = c(ncp.wald(resmod = m0, unresmod = m1),ncp.lr(resmod = m0, unresmod = m1),ncp.score(resmod = m0, unresmod = m1))
pw6 = power(df=null.hyp$df,ncp6,alpha=.05,ssize=nrow(dat6))
ssize6 = ssize(null.hyp$df,ncp6 ,alpha=.05,power=.80)


```


load("C:/Users/felix/Google Drive/4 irt/mmlpwr/results1split.Rdata")
res1.split2 = lapply(res1.split,function(x) x[names(x)!="fits"])
res1.split = res1.split2

# glue all with the same condition together again
res1 = lapply(sim1,function(x) {x$obs = list();x = x[names(x)!="datasets"];return(x)})
for (i in res1.split) {
  index = which(sapply(sim1, function(x) digest(x$condition)==digest(i$condition))) %>%as.numeric()
  res1[[index]]$obs = c(res1[[index]]$obs,i$obs)
}
save(res1,file= "C:/Users/felix/Google Drive/4 irt/mmlpwr/results1.Rdata")


# Plot the expected and observed power for all the above conditions. Add a CI-Envelope (Variance of these values is n*p*q ?)

# p = .05
# q = 1-p
# nruns = 500
# a = sqrt(nruns*p*q)/nruns * qnorm(.975)
# c(p-a,p+a) %>% round(.,3)

# res =c()
# for (i in 1:100000) {
#   runs = 500
#   res[i] = mean((runif(runs)<=p))
# }
# quantile(res,c(.025,.975))



# # Power - Calculate ncps  -----------------------------------------
#
# cl <- makeCluster(10)
# clusterExport(cl, objects(), envir=environment())
#
# res2ncp =  parLapply(cl,X = sim2ncp,fun = function(x){
#
#   load.libraries()
#
#   analytical = NULL
#
#   # calculate ncps for the condition
#   if (x$n.items==10) {
#     analytical = calculate_ncps(x$hyp,simbased=FALSE,n=1)
#   }
#   simbased = calculate_ncps(x$hyp,simbased=TRUE,n=1,n.pers = NULL)
#
#   # Sim Condition without the datasets
#   condition = x[names(x)!="datasets"]
#
#   re = list(condition = condition,analytical = analytical,simbased=simbased)
#   return(re)
#
# })
#
# save(res2ncp,file= "results2ncp.Rdata")
# stopCluster(cl)
#
#
#
# # Power - Calculate sample size -------------------------------------
#
# # add the number of persons and the datasets to the simulation conditions
#
# powerset=c(.2,.4,.6,.8)
#
# # calculate the sample sizes according to the ncps
# sample_sizes = lapply(res2ncp, function(x) {
#
#   if (!is.null(x$analytical)) {
#     ncps = x$analytical
#   } else {
#     ncps = x$simbased
#   }
#
#   re = sapply(powerset,function(pw) {ssize(hyp=x$condition$hyp,ncp=ncps,power = pw,alpha=.05) %>% mean()})
#   return(re)
# })
#
# #add the sample sizes
#
# res2 = list()
#
# for (i in 1:length(sim2ncp)) {
#
#   ss = sample_sizes[i]
#   a = sim2ncp[i]
#   for (x in ss){
#     b = a
#     b$n.pers = x
#     res2 = c(res2,b)
#   }
# }
#
# #add the datasets
#
# sim2 = lapply(sim2,function(x){
#
#   x$datasets = list()
#   for (i in 1:x$simruns) {
#     x$datasets[[i]] = setup.data(hyp = x$hyp, n = x$n.pers)
#   }
#   return(x)
#
# })
#
#
# # Power - Calculate stats --------------------------------------------------
#
#
# cl <- makeCluster(10)
# clusterExport(cl, objects(), envir=environment())
#
#
# res2 =  parLapply(cl,X = sim2,fun = function(x){
#
#   load.libraries()
#
#   # fit the mirt model
#   fitted <- mml.fit(data = x$data,hyp = x$hyp)
#   # calculate statistics
#   stats_obs <- c(wald_obs(fitted),lr_obs(fitted),score_obs(fitted))
#
#   # Sim Condition without the datasets
#   condition = x[names(x)!="datasets"]
#
#   re = list(condition = condition,fitted=fitted,stats_obs = stats_obs)
#   return(re)
#
# })
# save(res2,file= "results2.Rdata")
# stopCluster(cl)
#
#
#
# # Misspecification: Calculate ncp --------------------------------------------------------
#
# # Builds on the Power condition 1PLvs2PL, 10 Items
#
# # The ncps are therefore included in "results2ncp.Rdata"
#
# indices = sapply(res2ncp,function(x) {
#   re = x$condition$type = "1PLvs2PL" & x$condition$n.items = 10
#   return(re)
# })
#
# res3ncp = res2ncp[indices]
#
# save(res3ncp,file= "results3ncp.Rdata")
#
#
# # add sample sizes to the data.
#
# # Setup the simulation condition - based on sim2. datasets have to be generated anew for the non-normal distributions
#
# indices = sapply(sim2,function(x) {
#   re = x$condition$type = "1PLvs2PL" & x$condition$n.items = 10
#   return(re)
# })
#
# sim3 = list()
#
# for (x in sim2[indices]) {
#   a = x
#   b = x
#
#   a$dist.fun = unifdist
#   b$dist.fun = skeweddist
#
#   sim3 = c(sim3,a,b)
# }
#
#
# sim3 = lapply(sim3,function(x){
#
#   x$datasets = list()
#   for (i in 1:x$simruns) {
#     x$datasets[[i]] = setup.data(hyp = x$hyp, n = x$n.pers,dist.type=x$dist.fun)
#   }
#   return(x)
#
# })
#
#
#
# # Misspecification: Calculate stats -----------------------------------------
#
#
# cl <- makeCluster(10)
# clusterExport(cl, objects(), envir=environment())
#
#
# res3 =  parLapply(cl,X = sim3,fun = function(x){
#
#   load.libraries()
#
#   # fit the mirt model
#   fitted <- mml.fit(data = x$data,hyp = x$hyp)
#   # calculate statistics
#   stats_obs <- c(wald_obs(fitted),lr_obs(fitted),score_obs(fitted))
#
#   # Sim Condition without the datasets
#   condition = x[names(x)!="datasets"]
#
#   re = list(condition = condition,fitted=fitted,stats_obs = stats_obs)
#   return(re)
#
# })
# save(res3,file= "results3.Rdata")
# stopCluster(cl)
#
#
#
# # save updated sim files --------------------------------------------------
#
# save(sim2,sim3,file= "C:/Users/felix/Google Drive/4 irt/mmlpwr/simsetup2.Rdata")



# Power ------------------------------------------------------------

sim2 = expand.grid(type=c("1PLvs2PL","DIF2PL"),n.items=c(10,50),esize ="",simruns = 50)
sim2$esize[sim2$type =="1PLvs2PL"& sim2$n.items == 10] = "large"
sim2$esize[sim2$type =="DIF2PL"& sim2$n.items == 10] = "large"
sim2$esize[sim2$type =="1PLvs2PL"& sim2$n.items == 50] = "xsmall"
sim2$esize[sim2$type =="DIF2PL"& sim2$n.items == 50] = "small"

sim2 = sim2[order(sim2$type,sim2$n.items,sim2$esize,sim2$n.pers),]
sim2$condition = list(type= sim2$type,n.items = sim2$n.items,esize = sim2$esize)
sim2 = sim2  %>% split(., seq(nrow(.))) %>% lapply(.,as.list)


# Add the hypothesis object to the list according to the specifications of the conditions.
sim2ncp = lapply(sim2,function(x) {

  n.items = x$n.items
  esize = x$esize
  quantiles = (1:n.items)/(n.items+1)
  perm = sample(n.items)

  if (x$type == "1PLvs2PL") {

    if (esize=="large") {d = .1} else if (esize=="small") {d = .05} else if (esize=="xsmall") {d=0.03} else if (esize=="no") {d=0}

    x$altpars = list(
      a = qlnorm(quantiles,sdlog = d),
      d = qnorm(quantiles)[perm]
    )
  }

  if (x$type == "DIF2PL") {

    if (esize=="large") {d = .4} else if (esize=="small") {d = .2} else if (esize=="no") {d=0}

    group1 = group2 <- list(
      a = qlnorm(quantiles,sdlog = .1),
      d = qnorm(quantiles)[perm]
    )
    group1$a[1] = 1+d/2
    group1$d[1] = d/2
    group2$a[1] = 1-d/2
    group2$d[1] = -d/2

    x$altpars <- list(group1,group2)

  }

  x$hyp = setup_hypothesis(type = x$type, altpars = x$altpars)

  return(x)
})



# h_1PLvs2PL --------------------------------------------------------------

# Alternative A matrix to test for invariance - Item 1 minus all other


h_1PLvs2PLa = list(

  res = function(altpars,nullpars = NULL) {

    n.items = length(altpars[[1]])

    init = matrix(0,n.items-1,n.items*2)
    init[,1] = 1
    s = seq(3,n.items*2-1,2)
    for ( i in 1:nrow(init)) {
      init[i,s[i]] = -1
    }
    Amat = init

    re = list(
      n.items = n.items,
      itemtype = "2PL",
      Amat = Amat,
      cvec = 0,
      model = mirt::mirt.model(paste('F = 1-',n.items,'
                           CONSTRAIN = (1-',n.items,', a1)'))
    )
    return(re)
  },

  unres = function(altpars) {

    if (!is.null(altpars)) {
      re = altpars
      re$model = 1
      re$itemtype = "2PL"
    }

    re$longcoef = coef.long(shortcoef = re,itemtype=re$itemtype) # check if this is really needed.
    return(re)
  },

  maximizeL = function(hyp,bootstrap.start=TRUE) {
    # Hypothesis-specific algorithm to find the maximum likelihood restricted parameter set

    resmod = hyp$resmod
    unresmod = hyp$unresmod

    pars = unresmod
    load.functions(pars$itemtype)

    if(isTRUE(bootstrap.start)) {
      set.seed(1234)
      df = mirt::simdata(a = pars$a,d = pars$d,N =200000,itemtype = "2PL")
      mml = mirt(df,model = unresmod$model,itemtype = resmod$itemtype,technical = list(NCYCLES = 5000),verbose = FALSE)
      startval=c(coef_short(mml,"2PL")$a[1],coef_short(mml,"2PL")$d)
    }else{
      startval = c(mean(pars$a),as.numeric(pars$d))
    }

    maxl2pre = maxl2preload(pars)
    optpar = optim(startval,function(x) {maxl2(x,pars,maxl2pre)})
    re = pars
    re$a = rep(optpar$par[1],length(pars$a))
    re$d = optpar$par[2:length(optpar$par)]

    return(re)
  }
)


# Thissen Mat

ft = function(th,a,d,x) log(f(th,a,d,x))
ftd = Deriv::Deriv(ft,c("a","d"),nderiv=2)


tmat = function(pars) {
  res = list()
  for (i in 1:length(pars$a)) {
    res[[i]] = spatstat::gauss.hermite(function(th) {
      -matrix(ftd(th,pars$a[i],pars$d[i],0),2,2)
      # -matrix(ftd(th,pars$a[i],pars$d[i],p(th,pars$a[i],pars$d[i])),2,2)
    },order=30)
  }
  as.matrix(do.call(bdiag,res))
}



We can test it:
  ```{r}
altpars <- list(
  a = rlnorm(6,sdlog = .4),
  d = rnorm(6)
)

hyp <- setup_hypothesis(type = h_1PLvs2PL, altpars = altpars)

ptm <- proc.time() #start timer
ncps <- calculate_ncps(hyp=hyp)
time = proc.time() - ptm #stop timer

time = as.numeric(time[1])
time

```


# New facet label names for supp variable
# esize.labs <- c(expression("H"[0]*": No effect"), expression("H"[1]*": Small effect"),expression("H"[1]*": Large effect"))
# esize.labs <- function(st) {
#   c(expression("H"[0]*": No effect"), expression("H"[1]*": Small effect"),expression("H"[1]*": Large effect"))  }




