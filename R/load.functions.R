

#' Load IRT model functions and derivatives
#'
#' This is a helper function used to generate custom hypotheses. See the "adding_hypotheses" vignette.
#'
#' @param model character, Desired Model (2PL,3PL, GPCM)
#' @param multi logical, multidimensional model if TRUE (available for 2PL)
#'
#' @return nothing
#' @export
#'
#' @examples
#'
#' funs = load.functions("2PL")
#'
load.functions = function (model,multi=FALSE) {


# 2PL ---------------------------------------------------------------------


  # 2PL
  if (model=="2PL"&!multi) {

    f_2pl = function(th,a,d,x) exp(x*(a*th+d))/(1+exp(a*th+d))
    f1 =  Deriv::Deriv(f_2pl,c("a"))
    f2 =  Deriv::Deriv(f_2pl,c("d"))

    g = function(pattern, pars) {
      re = spatstat.random::gauss.hermite(function(th) {
        prod(f_2pl(th,pars$a,pars$d,pattern))
      }, order = 30)
      return(re)
    }
    fdot_2pl = function(th,a,d,x) {
      res = f_2pl(th,a,d,x)
      temp1 = f1(th,a,d,x)
      temp2 = f2(th,a,d,x)

      a = lapply(seq_len(length(x)),function(i) {
       prod(res[-i])* c(temp1[i],temp2[i])
      })
      as.numeric(unlist(a))
    }
    gdot = function(pattern, pars) {
      re = spatstat.random::gauss.hermite(function(th) {
        fdot_2pl(th,pars$a,pars$d,pattern)
      }, order = 30)
      return(re)
    }
    ldot = function(pattern,pars) {
      re=gdot(pattern,pars)/g(pattern,pars)
      return(re)
    }

    funlist=list(f=f_2pl,g=g,fdot=fdot_2pl,gdot=gdot,ldot=ldot)

  }


# 2PL multi ---------------------------------------------------------------


  # 2PL
  if (model=="2PL"& multi) {

    f_2pl_multi = function(th1,th2,a1,a2,d,x) exp(x*(a1*th1+a2*th2+d))/(1+exp(a1*th1+a2*th2+d))
    f1 =  Deriv::Deriv(f_2pl_multi,c("a1"))
    f2 =  Deriv::Deriv(f_2pl_multi,c("a2"))
    f3 =  Deriv::Deriv(f_2pl_multi,c("d"))

    g = function(pattern, pars) {
      re = spatstat.random::gauss.hermite(function(th1) {
      re = spatstat.random::gauss.hermite(function(th2) {
        prod(f_2pl_multi(th1,th2,pars$a[,1],pars$a[,2],pars$d,pattern))
      }, order = 30)},order=30)
      return(re)
    }
    fdot_2pl_multi = function(th1,th2,a1,a2,d,x) {
      res = f_2pl_multi(th1,th2,a1,a2,d,x)
      temp1 = f1(th1,th2,a1,a2,d,x)
      temp2 = f2(th1,th2,a1,a2,d,x)
      temp3 = f3(th1,th2,a1,a2,d,x)

      a = lapply(seq_len(length(x)),function(i) {
        prod(res[-i])* c(temp1[i],temp2[i],temp3[i])
      })
      as.numeric(unlist(a))
    }
    gdot = function(pattern, pars) {
      re = spatstat.random::gauss.hermite(function(th1) {
        re = spatstat.random::gauss.hermite(function(th2) {
          fdot_2pl_multi(th1,th2,pars$a[,1],pars$a[,2],pars$d,pattern)
        }, order = 30)},order=30)
      return(re)
    }
    ldot = function(pattern,pars) {
      re=gdot(pattern,pars)/g(pattern,pars)
      return(re)
    }

    funlist=list(f=f_2pl_multi,g=g,fdot=fdot_2pl_multi,gdot=gdot,ldot=ldot)

  }


# 3PL ---------------------------------------------------------------------

  # 3PL
  if (model=="3PL"&!multi) {

    f_3pl = function(th,a,d,g,x) (1-x) + (2 * x - 1) * logitinv(g) + (2 * x - 1) * (1-logitinv(g)) /(1+exp(-(a*th+d)))

    f1 =  Deriv::Deriv(f_3pl,c("a"))
    f2 =  Deriv::Deriv(f_3pl,c("d"))
    f3 =  Deriv::Deriv(f_3pl,c("g"))

    g = function(pattern, pars) {
      re = spatstat.random::gauss.hermite(function(th) {
        prod(f_3pl(th,pars$a,pars$d,pars$g,pattern))
      }, order = 30)
      return(re)
    }
    fdot_3pl = function(th,a,d,g,x) {
      res = f_3pl(th,a,d,g,x)
      temp1 = f1(th,a,d,g,x)
      temp2 = f2(th,a,d,g,x)
      temp3 = f3(th,a,d,g,x)

      a = lapply(seq_len(length(x)),function(i) {
        prod(res[-i])* c(temp1[i],temp2[i],temp3[i])
      })
      as.numeric(unlist(a))
    }
    gdot = function(pattern, pars) {
      re = spatstat.random::gauss.hermite(function(th) {
        fdot_3pl(th,pars$a,pars$d,pars$g,pattern)
      }, order = 30)
      return(re)
    }
    ldot = function(pattern,pars) {
      re=gdot(pattern,pars)/g(pattern,pars)
      return(re)
    }

    funlist=list(f=f_3pl,g=g,fdot=fdot_3pl,gdot=gdot,ldot=ldot)

  }



# GPCM --------------------------------------------------------------------


  # GPCM
  if (model=="gpcm") {

    f_gpcm = function(th,a,d1,d2,x) exp(x*((2-x)*(a*th+d1)+(x-1)/2*(2*a*th+d2)))/(1+exp(a*th+d1)+exp(2*a*th+d2))
    f1 =  Deriv::Deriv(f_gpcm,c("a"))
    f2 =  Deriv::Deriv(f_gpcm,c("d1"))
    f3 =  Deriv::Deriv(f_gpcm,c("d2"))

    g = function(pattern, pars) {
      re = spatstat.random::gauss.hermite(function(th) {
        prod(f_gpcm(th,pars$a,pars$d[,2],pars$d[,3],pattern))
      }, order = 50)
      return(re)
    }

    fdot_gpcm = function(th,a,d1,d2,x) {
      res = f_gpcm(th,a,d1,d2,x)
      temp1 = f1(th,a,d1,d2,x)
      temp2 = f2(th,a,d1,d2,x)
      temp3 = f3(th,a,d1,d2,x)

      a = lapply(seq_len(length(x)),function(i) {
        prod(res[-i])* c(temp1[i],temp2[i],temp3[i])
      })
      as.numeric(unlist(a))
    }

    gdot = function(pattern, pars) {
      re = spatstat.random::gauss.hermite(function(th) {
        fdot_gpcm(th,pars$a,pars$d[,2],pars$d[,3],pattern)
      }, order = 50)
      return(re)
    }
    ldot = function(pattern,pars) {
      re=gdot(pattern,pars)/g(pattern,pars)
      return(re)
    }

    funlist=list(f=f_gpcm,g=g,fdot=fdot_gpcm,gdot=gdot,ldot=ldot)
  }


  return(funlist)

}


