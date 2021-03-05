

# Presets containing:

# - Functions to generate a restricted and unrestricted model
# - Function to maximize the Likelihood
# currently implemented: 1PLvs2PL, DIF2PL


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

  unres = function(altpars) {

    re = list(
      parsets = altpars,
      model = 1,
      itemtype = "2PL",
      longpars = pars.long(pars = altpars,itemtype="2PL")
    )

    return(re)
  },

  maximizeL = function(hyp) {
    # Hypothesis-specific algorithm to find the maximum likelihood restricted parameter set


    maxlpreload= function(pars) {
      # returns the density for each response pattern under the model parameters pars

      patterns = as.matrix(expand.grid(lapply(1:length(pars$a),function(x) c(0,1))))

      pre = c()
      for (i in 1:nrow(patterns)) {
        pre[i] = g(patterns[i,],pars)
      }

      return(pre)
    }


    maxl = function(x,pars,pre) {
      # calculates the likelihood of parameters x given model "pars"

      patterns = as.matrix(expand.grid(lapply(1:length(pars$a),function(x) c(0,1))))
      x = list(a=rep(x[1],length(pars$a)),d=x[2:length(x)])

      res  = c()
      for (i in 1:nrow(patterns)) {
        px = pre[i]
        qx = g(patterns[i,],x)
        res[i] =  {px*log(qx)}
      }
      re = -sum(res)
    }
    resmod = hyp$resmod
    unresmod = hyp$unresmod

    pars = unresmod$parsets
    load.functions(unresmod$itemtype)

    startval = c(mean(pars$a),as.numeric(pars$d))

    maxlpre = maxlpreload(pars)
    optpar = optim(startval,function(x) {maxl(x,pars,maxlpre)},method = "BFGS")
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
                       CONSTRAINB = (1-',n.items,', d), (1-',n.items,', a1)')),
      multigroup = TRUE
    )

    return(re)
  },

  unres = function(altpars) {

    n.items = length(altpars[[1]][[1]])

    reA = altpars[[1]]
    reB = altpars[[2]]

    reA$itemtype = reB$itemtype = "2PL"

    reA$longpars = pars.long(pars = reA,itemtype="2PL")
    reB$longpars = pars.long(pars = reB,itemtype="2PL")
    re = list(parsets = list(reA,reB),
              model = mirt::mirt.model(paste('F = 1-',n.items,'
                     CONSTRAINB = (2-',n.items,', d), (2-',n.items,', a1)')),
              longpars = c(reA$longpars,reB$longpars[1:2]),
              multigroup = TRUE,
              itemtype = "2PL"
    )

    return(re)
  },

  maximizeL = function(hyp) {
    # Hypothesis-specific algorithm to find the maximum likelihood restricted parameter set
    # L Optimizer

    maxl = function(x,pars1,pars2,i) {

      px1 = function(th) {f(th,pars1$a[i],pars1$d[i],1)}
      px2 = function(th) {f(th,pars2$a[i],pars2$d[i],1)}
      qx = function(th) {f(th,x[1],x[2],1)}
      kl = function(th) {px1(th)*log(qx(th))+(1-px1(th))*log(1-qx(th)) + px2(th)*log(qx(th))+(1-px2(th))*log((1-qx(th)))
      }
      re = -spatstat::gauss.hermite(kl,order=20)
    }

    resmod = hyp$resmod
    unresmod = hyp$unresmod

    pars = unresmod$parsets
    load.functions(unresmod$itemtype)

    pars1 = pars[[1]]
    pars2 = pars[[2]]

    load.functions(pars1$itemtype)
    re = pars1

    for (i in 1:length(pars1$a)) {
      startval = c(re$a[i],re$d[i])
      optpar = optim(startval,function(x) {maxl(x,pars1,pars2,i)},method = "BFGS")
      re$a[i] = optpar$par[1]
      re$d[i] = optpar$par[2]
    }

    return(re)
  }
)


# h_basic --------------------------------------------------------------

#hypothesis that the first item has difficulty 0

h_basic = list(

  res = function(altpars,nullpars = NULL) {

    n.items = length(altpars[[1]])

    re = list(
      n.items = n.items,
      itemtype = "2PL",
      Amat = c(0,1,rep(0,(n.items-1)*2)) %>% matrix(.,ncol=n.items*2,byrow=TRUE),
      cvec = 0,
      model = mirt::mirt.model(paste('F = 1-',n.items,'
                           FIXED = (1, d)
                           START = (1,d,0)'))
    )
    return(re)
  },

  unres = function(altpars) {

    re = list(
      parsets = altpars,
      model = 1,
      itemtype = "2PL",
      longpars = pars.long(pars = altpars,itemtype="2PL")
    )

    return(re)
  },

  maximizeL = function(hyp) {
    # Hypothesis-specific algorithm to find the maximum likelihood restricted parameter set


    maxlpreload= function(pars) {
      # returns the density for each response pattern under the model parameters pars

      patterns = as.matrix(expand.grid(lapply(1:length(pars$a),function(x) c(0,1))))

      pre = c()
      for (i in 1:nrow(patterns)) {
        pre[i] = g(patterns[i,],pars)
      }

      return(pre)
    }


    maxl = function(x,pars,pre) {
      # calculates the likelihood of parameters x given model "pars"
      patterns = as.matrix(expand.grid(lapply(1:length(pars$a),function(x) c(0,1))))

      x = list(
        a=c(x,pars$a[2:length(pars$a)]),
        d=c(0,pars$d[2:length(pars$d)])
        )

      res  = c()
      for (i in 1:nrow(patterns)) {
        px = pre[i]
        qx = g(patterns[i,],x)
        res[i] =  {px*log(qx)}
      }
      re = -sum(res)
    }
    resmod = hyp$resmod
    unresmod = hyp$unresmod

    pars = unresmod$parsets
    load.functions(unresmod$itemtype)

    startval = pars$a[1]

    maxlpre = maxlpreload(pars)

    optpar = optim(startval,function(x) {maxl(x,pars,maxlpre)},method = "BFGS")
    re = pars
    re$a = c(optpar$par[1],pars$a[2:length(pars$a)])
    re$d = c(0,pars$d[2:length(pars$d)])

    return(re)
  }
)
