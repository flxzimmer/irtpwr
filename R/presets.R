


# h_1PLvs2PL --------------------------------------------------------------

h_1PLvs2PL = list(

  res = function(altpars,nullpars = NULL) {

    n.items = length(altpars[[1]])

    re = list(
      n.items = n.items,
      itemtype = "2PL",
      Amat = c(1,0,-1,0,rep(0,(n.items-1)*2)) %>% rep(.,n.items-2) %>% c(.,1,0,-1,0) %>% matrix(.,ncol=n.items*2,byrow=TRUE),
      cvec = 0,
      model = mirt::mirt.model(paste('F = 1-',n.items,'
                           CONSTRAIN = (1-',n.items,', a1)'))
    )
    return(re)
  },

  unres = function(esize=NULL,altpars = NULL) {

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
      df = simdata(a = pars$a,d = pars$d,N =200000,itemtype = "2PL")
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



# h_DIF2PL ----------------------------------------------------------------

h_DIF2PL = list(

  res = function(altpars,nullpars = NULL) {

    n.items = length(altpars[[1]][[1]])

    re = list(
      n.items = n.items,
      itemtype = "2PL",
      Amat = c(1,0,rep(0,(n.items-1)*2),-1,0) %>% c(.,c(0,1,rep(0,(n.items-1)*2),0,-1)) %>% matrix(.,ncol=(n.items+1)*2,byrow=TRUE),
      cvec = 0,
      model = mirt::mirt.model(paste('F = 1-',n.items,'
                       CONSTRAINB = (2-',n.items,', d), (2-',n.items,', a1)'))
    )
    return(re)
  },

  unres = function(esize=NULL,altpars = NULL) {

    if (!is.null(altpars)) {

      reA = altpars[[1]]
      reB = altpars[[2]]

      reA$itemtype = "2PL"
      reB$itemtype = "2PL"

      reA$longcoef = coef.long(shortcoef = reA,itemtype=reA$itemtype)
      reB$longcoef = coef.long(shortcoef = reB,itemtype=reB$itemtype)
      re = list(parsets = list(reA,reB), model = 1)

      }

    re$longcoef = c(reA$longcoef,reB$longcoef[1:2])
    re$multigroup = TRUE
    re$itemtype = "2PL"

    return(re)
  },

  maximizeL = function(hyp,bootstrap.start=TRUE) {
    # Hypothesis-specific algorithm to find the maximum likelihood restricted parameter set

    resmod = hyp$resmod
    unresmod = hyp$unresmod

    pars = unresmod
    load.functions(pars$itemtype)

    pars1 = pars[[1]][[1]]
    pars2 = pars[[1]][[2]]

    load.functions(pars1$itemtype)
    re = pars1

    for (i in 1:length(pars1$a)) {
      startval = c(re$a[i],re$d[i])
      optpar = optim(startval,function(x) {maxl(x,pars1,pars2,i)})
      re$a[i] = optpar$par[1]
      re$d[i] = optpar$par[2]
    }


    return(re)
  }
)



# Other -------------------------------------------------------------------


# setup.null.hypothesis = function(itemtype=NULL,n.items,model = 1,Amat=NULL, cvec=NULL, nullpars =NULL,preset=NULL) {
#   #itemtype is the type of the unrestricted model
#   # model is a mirt model for the unrestricted model (default 1, can also be multigroup)
#   # Amat is the A matrix for the null hypothesis
#   # cvec is the c vector for the null hypothesis
#   # nullpars are the parameters under the null hypothesis
#   # preset are presets
#
#   if (preset == "1PLvs2PL") {
#     itemtype = "2PL"
#     Amat = c(1,0,-1,0,rep(0,(n.items-1)*2)) %>% rep(.,n.items-2) %>% c(.,1,0,-1,0) %>% matrix(.,ncol=n.items*2,byrow=TRUE)
#     cvec = 0
#     model = 1
#     multigroup = FALSE
#   }
#
#   if (preset == "PCMvsGPCM") {
#     itemtype = "gpcm"
#     nkat = 3
#     Amat = c(1,rep(0,nkat-1),-1,rep(0,nkat-1),rep(0,(n.items-1)*nkat)) %>% rep(.,n.items-2) %>% c(1,rep(0,nkat-1),-1,rep(0,nkat-1)) %>% matrix(.,ncol=n.items*nkat,byrow=TRUE)
#     cvec = 0
#     model = 1
#     multigroup = FALSE
#   }
#
#   if (preset == "3PLspec") {
#     itemtype = "3PL"
#     Amat = c(0,0,1,rep(0,(n.items)*3)) %>% rep(.,n.items-1) %>% c(.,0,0,1) %>% matrix(.,ncol=n.items*3,byrow=TRUE)
#     cvec = rep(.2,n.items)
#     model = 1
#     multigroup = FALSE
#   }
#
#   if (preset == "DIF2PL") {
#     itemtype = "2PL"
#     Amat = c(1,0,rep(0,(n.items-1)*2),-1,0) %>% c(.,c(0,1,rep(0,(n.items-1)*2),0,-1)) %>% matrix(.,ncol=(n.items+1)*2,byrow=TRUE)
#     cvec = 0
#     model = mirt.model(paste('F = 1-',n.items,'
#                       CONSTRAINB = (2-',n.items,', d), (2-',n.items,', a1)'))
#     multigroup = TRUE
#   }
#
#   re = list(itemtype=itemtype,
#             n.items = n.items,
#             model = model,
#             Amat = Amat,
#             cvec = cvec,
#             preset = preset,
#             multigroup = multigroup,
#             df = nrow(Amat))
#
#   return(re)
# }


# setup.alternative.hypothesis = function(null.hypothesis = NULL,esize= NULL,altpars= NULL,use.altpars = FALSE) {
#   # altpars are the parameters under the alterNULLtive hypothesis
#
#   n.items = null.hypothesis$n.items
#   set.seed(123)
#
#   if (is.null(esize)&!is.null(altpars)) {
#     use.altpars = TRUE
#     esize = "no" # setting esize temporarily
#   }
#
#   if (null.hypothesis$preset=="1PLvs2PL") {
#
#     quantiles = (1:n.items)/(n.items+1)
#     perm = sample(n.items)
#
#     if (esize=="large") {
#       apar = qlnorm(quantiles,sdlog = .1)
#     } else if (esize=="small") {
#       apar = qlnorm(quantiles,sdlog = .05)
#     } else if (esize=="xsmall") {
#       apar = qlnorm(quantiles,sdlog = .03)
#     } else if (esize=="no") {
#       apar = rep(1,n.items)}
#
#     re = list(
#       itemtype  = "2PL",
#       a = apar,
#       d = qnorm(quantiles)[perm],
#       model = mirt.model(paste('F = 1-',n.items,'
#                        CONSTRAIN = (1-',n.items,', a1)'))
#     )
#     re$longcoef = coef.long(shortcoef = re,itemtype=re$itemtype)
#
#   } else if (null.hypothesis$preset=="PCMvsGPCM") {
#     nkat= 3
#     quantiles = (1:n.items)/(n.items+1)
#     quantiles2 = (1:(n.items*(nkat-1)))/((n.items*(nkat-1))+1)
#     perm = sample(n.items)
#
#     if (esize=="large") {
#       apar = qlnorm(quantiles,sdlog = .1)
#     } else if (esize=="small") {
#       apar = qlnorm(quantiles,sdlog = .05)
#     } else if (esize=="no") {
#       apar = rep(1,n.items)}
#
#     re = list(
#       itemtype = "gpcm",
#       a = apar,
#       d = cbind(rep(0,n.items),matrix(qnorm(quantiles2),ncol=(nkat-1))[perm,]),
#       model = mirt.model(paste('F = 1-',n.items,'
#                        CONSTRAIN = (1-',n.items,', a1)'))
#     )
#     re$longcoef = coef.long(shortcoef = re,itemtype=re$itemtype)
#
#   } else if (null.hypothesis$preset=="3PLspec") {
#     quantiles = (1:n.items)/(n.items+1)
#     perm = sample(n.items)
#     perm2 = sample(n.items)
#     if (esize=="large") {mult = .4} else if (esize=="small") {mult=.2} else if (esize=="no") {mult=0}
#
#     re = list(
#       itemtype = "3PL",
#       a = qlnorm(quantiles,sdlog = .05),
#       d = qnorm(quantiles)[perm],
#       g = qlnorm(quantiles,sdlog=.9)[perm2]^mult/5,
#       model = mirt.model(paste('F = 1-',n.items,'
#                        START = (1-',n.items,', g,.2),
#                        FIXED = (1-',n.items,', g)')))
#     re$longcoef = coef.long(shortcoef = re,itemtype=re$itemtype)
#
#   } else if (null.hypothesis$preset=="DIF2PL") {
#     quantiles = (1:n.items)/(n.items+1)
#     perm = sample(n.items)
#     if (esize=="large") {d = .4} else if (esize=="small") {d = .2} else if (esize=="no") {d=0}
#
#     reA = reB = list(
#       itemtype  = "2PL",
#       a = qlnorm(quantiles,sdlog = .1),
#       d = qnorm(quantiles)[perm]
#     )
#     reA$a[1] = 1+d/2
#     reA$d[1] = d/2
#     reB$a[1] = 1-d/2
#     reB$d[1] = -d/2
#     reA$longcoef = coef.long(shortcoef = reA,itemtype=reA$itemtype)
#     reB$longcoef = coef.long(shortcoef = reB,itemtype=reB$itemtype)
#     re = list(parsets = list(reA,reB), model = 1)
#     re$longcoef = c(reA$longcoef,reB$longcoef[1:2])
#   }
#
#   if (isTRUE(use.altpars)) {
#     re$a = altpars$a
#     re$d = altpars$d
#     re$g = altpars$g
#     re$longcoef = coef.long(shortcoef = re,itemtype=re$itemtype)
#   }
#   return(re)
# }



# maximizeL = function(null.hypothesis,alternative.hypothesis,bootstrap.start=TRUE) {
#   #outputs pars that follow the restriction while exhibiting minimal KL distance
#   print("maximizing L")
#   if (null.hypothesis$preset=="1PLvs2PL") {
#
#     pars = alternative.hypothesis
#     load.functions(pars$itemtype)
#
#     if(isTRUE(bootstrap.start)) {
#       df = simdata(a = pars$a,d = pars$d,N =200000,itemtype = "2PL")
#       mml = mirt(df,model = alternative.hypothesis$model,itemtype = null.hypothesis$itemtype,technical = list(NCYCLES = 5000),verbose = FALSE)
#       startval=c(coef_short(mml,"2PL")$a[1],coef_short(mml,"2PL")$d)
#     }else{
#       startval = c(mean(pars$a),as.numeric(pars$d))
#     }
#
#     maxl2pre = maxl2preload(pars)
#     optpar = optim(startval,function(x) {maxl2(x,pars,maxl2pre)})
#     re = pars
#     re$a = rep(optpar$par[1],length(pars$a))
#     re$d = optpar$par[2:length(optpar$par)]
#   }
#
#   if (null.hypothesis$preset=="PCMvsGPCM") {
#
#     pars = alternative.hypothesis
#     load.functions(pars$itemtype)
#     startval = c(mean(pars$a),as.numeric(pars$d[(null.hypothesis$n.items+1):length(as.numeric(pars$d))]))
#     maxl2pre = maxl2preload(pars)
#     optpar = optim(startval,function(x) {maxl2(x,pars,maxl2pre)})
#     re = pars
#     re$a = rep(optpar$par[1],length(pars$a))
#     re$d = matrix(c(rep(0,null.hypothesis$n.items),optpar$par[2:length(optpar$par)]),ncol=ncol(pars$d))
#
#   }
#
#   if (null.hypothesis$preset=="3PLspec") {
#     pars = alternative.hypothesis
#     load.functions(pars$itemtype)
#     re = pars
#     re$g=rep(0.2,null.hypothesis$n.items)
#
#     for (i in 1:length(pars$a)) {
#       startval = c(re$a[i],re$d[i])
#       optpar = optim(startval,function(x) {maxl(x,pars,i,.2)})
#       re$a[i] = optpar$par[1]
#       re$d[i] = optpar$par[2]
#     }
#   }
#
#   if (null.hypothesis$preset=="DIF2PL") {
#
#     pars1 = alternative.hypothesis[[1]][[1]]
#     pars2 = alternative.hypothesis[[1]][[2]]
#
#     load.functions(pars1$itemtype)
#     re = pars1
#
#     for (i in 1:length(pars1$a)) {
#       startval = c(re$a[i],re$d[i])
#       optpar = optim(startval,function(x) {maxl(x,pars1,pars2,i)})
#       re$a[i] = optpar$par[1]
#       re$d[i] = optpar$par[2]
#     }
#
#   }
#   print("finished maximizing L")
#
#   return(re)
# }



