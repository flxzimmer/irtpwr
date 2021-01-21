

# Presets containing:

# - Functions to generate a restricted and unrestricted model
# - Function to maximize the Likelihood
# - Funciton to calculate runtime
#
# currently implemented: 1PLvs2PL, DIF2PL
#
# Tutorial on how to create new presets is upcoming

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


    maxl2preload= function(pars) {
      # returns the density for each response pattern under the model parameters pars

      patterns = as.matrix(expand.grid(lapply(1:length(pars$a),function(x) c(0,1))))

      pre = c()
      for (i in 1:nrow(patterns)) {
        pre[i] = g(patterns[i,],pars)
      }

      return(pre)
    }


    maxl2 = function(x,pars,pre) {
      # calculates the likelihood of parameters x given model "pars"

      patterns = as.matrix(expand.grid(lapply(1:length(pars$a),function(x) c(0,1))))
      #patterns = split(patterns, rep(1:nrow(patterns), each = ncol(patterns)))

      x = list(a=rep(x[1],length(pars$a)),d=x[2:length(x)])

      res  = c()
      for (i in 1:nrow(patterns)) {
        px = pre[i]
        qx = g(patterns[i,],x)
        res[i] =  {px*log(qx)}
      }

      # px = pre
      # qx = sapply(patterns,function(y) g(as.numeric(y),x))
      # re = -sum(px*log(qx))
      re = -sum(res)
    }

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
    set.seed(1234)

    maxl2pre = maxl2preload(pars)
    optpar = optim(startval,function(x) {maxl2(x,pars,maxl2pre)})
    re = pars
    re$a = rep(optpar$par[1],length(pars$a))
    re$d = optpar$par[2:length(optpar$par)]

    return(re)
  },
  calctime = function(n.items) {

    altpars <- list(
      a = rlnorm(3,sdlog = .4),
      d = rnorm(3)
    )

    hyp <- setup_hypothesis(type = h_1PLvs2PL, altpars = altpars)

    ptm <- proc.time() #start timer
    ncps <- calculate_ncps(hyp=hyp)
    time = proc.time() - ptm #stop timer

    time = as.numeric(time[1])

    # time is t * 2^3
    re = time *2^(n.items-3)

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
      # model = mirt::mirt.model(paste('F = 1-',n.items,'
      #                  CONSTRAINB = (2-',n.items,', d), (2-',n.items,', a1)'))
      model =1
    )
    return(re)
  },

  unres = function(altpars) {

    n.items = length(altpars[[1]][[1]])

    if (!is.null(altpars)) {

      reA = altpars[[1]]
      reB = altpars[[2]]

      reA$itemtype = "2PL"
      reB$itemtype = "2PL"

      reA$longcoef = coef.long(shortcoef = reA,itemtype=reA$itemtype)
      reB$longcoef = coef.long(shortcoef = reB,itemtype=reB$itemtype)
      re = list(parsets = list(reA,reB),
                model = mirt::mirt.model(paste('F = 1-',n.items,'
                       CONSTRAINB = (2-',n.items,', d), (2-',n.items,', a1)')))
      }

    re$longcoef = c(reA$longcoef,reB$longcoef[1:2])
    re$multigroup = TRUE
    re$itemtype = "2PL"

    return(re)
  },

  maximizeL = function(hyp,bootstrap.start=TRUE) {
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
  },
  calctime = function(n.items) {

    group1 = group2 <- list(
      a = rlnorm(8,sdlog = .2),
      d = rnorm(8)
    )

    group2$a[1] = (group2$a[1])^2
    group2$d[1] = group2$d[1] + .5

    altpars <- list(group1,group2)

    hyp <- setup_hypothesis(type = h_DIF2PL, altpars = altpars)

    ptm <- proc.time() #start timer
    ncps <- calculate_ncps(hyp=hyp)
    time = proc.time() - ptm #stop timer

    time = as.numeric(time[1])

    # time is t * 2^3
    re = time *2^(n.items-8)

    return(re)

  }
)


