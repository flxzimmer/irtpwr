load.functions = function (model) {
  # Loads model formula

  # 2PL ---------------------------------------------------------------------

  if (model=="2PL") {

    p = function(th,a,d) 1/(1+exp(-(a*th+d)))
    f = function(th,a,d,x) {c(1-p(th,a,d),p(th,a,d))[x+1]}
    f1 =  Deriv::Deriv(f,c("a"))
    f2 =  Deriv::Deriv(f,c("d"))
    fvec = function(pattern,th,pars) {
      re = c()
      for (i in 1:length(pars$a)) {
        re[i] = f(th,pars$a[i],pars$d[i],pattern[i])
      }
      return(re)
    }
    g = function(pattern, pars) {
      re = spatstat::gauss.hermite(function(th) {
        prod(fvec(pattern, th, pars))
      }, order = 30)
      return(re)
    }
    fdot = function(pattern,th,pars) {
      res = fvec(pattern,th,pars)
      re =c()
      for (i in 1:length(pars$a)) {
        temp1 = res
        temp2 = res
        temp1[i] = sum(f1(th,pars$a[i],pars$d[i],pattern[i]),0,na.rm=T)
        temp2[i] = sum(f2(th,pars$a[i],pars$d[i],pattern[i]),0,na.rm=T)
        re = c(re,c(prod(temp1),prod(temp2)))
      }
      return(re)
    }
    gdot = function(pattern, pars) {
      re = spatstat::gauss.hermite(function(th) {
        fdot(pattern, th, pars)
      }, order = 30)
      return(re)
    }
    ldot = function(pattern,pars) {
      re=gdot(pattern,pars)/g(pattern,pars)
      return(re)
    }

    # # Thissen Mat
    #
    # ft = function(th,a,d,x) log(f(th,a,d,x))
    # ftd = Deriv::Deriv(ft,c("a","d"),nderiv=2)
    #
    #
    # tmat = function(pars) {
    #   res = list()
    #   for (i in 1:length(pars$a)) {
    #     res[[i]] = spatstat::gauss.hermite(function(th) {
    #       -matrix(ftd(th,pars$a[i],pars$d[i],0),2,2)
    #       # -matrix(ftd(th,pars$a[i],pars$d[i],p(th,pars$a[i],pars$d[i])),2,2)
    #     },order=30)
    #   }
    #   as.matrix(do.call(bdiag,res))
    # }


    # # L Optimizer
    #
    # maxl = function(x,pars1,pars2,i) {
    #   px1 = function(th) {f(th,pars1$a[i],pars1$d[i],1)}
    #   px2 = function(th) {f(th,pars2$a[i],pars2$d[i],1)}
    #   qx = function(th) {f(th,x[1],x[2],1)}
    #   kl = function(th) {px1(th)*log(qx(th))+(1-px1(th))*log(1-qx(th)) + px2(th)*log(qx(th))+(1-px2(th))*log((1-qx(th)))
    #   }
    #   re = -spatstat::gauss.hermite(kl,order=20)
    # }


  #   maxl2preload= function(pars) {
  # # returns the density for each response pattern under the model parameters pars
  #
  #     patterns = as.matrix(expand.grid(lapply(1:length(pars$a),function(x) c(0,1))))
  #
  #     pre = c()
  #     for (i in 1:nrow(patterns)) {
  #       pre[i] = g(patterns[i,],pars)
  #     }
  #
  #     return(pre)
  #
  #   }
  #
  #
  #   maxl2 = function(x,pars,pre) {
  #
  #     patterns = as.matrix(expand.grid(lapply(1:length(pars$a),function(x) c(0,1))))
  #     #patterns = split(patterns, rep(1:nrow(patterns), each = ncol(patterns)))
  #
  #     x = list(a=rep(x[1],length(pars$a)),d=x[2:length(x)])
  #
  #     res  = c()
  #     for (i in 1:nrow(patterns)) {
  #       px = pre[i]
  #       qx = g(patterns[i,],x)
  #       res[i] =  {px*log(qx)}
  #     }
  #
  #     # px = pre
  #     # qx = sapply(patterns,function(y) g(as.numeric(y),x))
  #     # re = -sum(px*log(qx))
  #     re = -sum(res)
  #   }

    # funlist=list(f=f,fvec=fvec,g=g,fdot=fdot,gdot=gdot,ldot=ldot,tmat=tmat,maxl=maxl,maxl2preload=maxl2preload,maxl2=maxl2)
    funlist=list(f=f,fvec=fvec,g=g,fdot=fdot,gdot=gdot,ldot=ldot)


  }

  # if (model=="2PLa") {
  #   source('mmlpwr_func_2PLa.R')
  # }
  #
  # if (model=="3PL") {
  #   source('mmlpwr_func_3PL.R')
  # }
  #
  # if (model=="3PL_nologit") {
  #   source('mmlpwr_func_3PL_nologit.R')
  # }
  # if (model=="3PL_logit") {
  #   source('mmlpwr_func_3PL_logit.R')
  # }
  #
  # if (model=="gpcm") {
  #   source('mmlpwr_func_GPCM.R')
  # }

  attach(funlist,warn.conflicts = FALSE)

}
