

#' Load IRT model functions and derivatives
#'
#' @param model Desired Model (2PL only at the moment)
#'
#' @return
#' @export
#'
#' @examples
load.functions = function (model,multi=FALSE) {



# 2PL ---------------------------------------------------------------------


  # 2PL
  if (model=="2PL"&!multi) {

    # p = function(th,a,d) 1/(1+exp(-(a*th+d)))
    f = function(th,a,d,x) exp(x*(a*th+d))/(1+exp(a*th+d))
    f1 =  Deriv::Deriv(f,c("a"))
    f2 =  Deriv::Deriv(f,c("d"))

    g = function(pattern, pars) {
      re = spatstat.core::gauss.hermite(function(th) {
        prod(f(th,pars$a,pars$d,pattern))
      }, order = 30)
      return(re)
    }
    fdot = function(th,a,d,x) {
      res = f(th,a,d,x)
      temp1 = f1(th,a,d,x)
      temp2 = f2(th,a,d,x)
      # clear NAs? can also do later and not in each function eval?

      a = lapply(1:length(x),function(i) {
       prod(res[-i])* c(temp1[i],temp2[i])
      })
      as.numeric(unlist(a))
    }
    gdot = function(pattern, pars) {
      re = spatstat.core::gauss.hermite(function(th) {
        fdot(th,pars$a,pars$d,pattern)
      }, order = 30)
      return(re)
    }
    ldot = function(pattern,pars) {
      re=gdot(pattern,pars)/g(pattern,pars)
      return(re)
    }

    funlist=list(f=f,g=g,fdot=fdot,gdot=gdot,ldot=ldot)

  }


# 2PL multi ---------------------------------------------------------------


  # 2PL
  if (model=="2PL"& multi) {

    # p = function(th,a,d) 1/(1+exp(-(a*th+d)))
    f = function(th1,th2,a1,a2,d,x) exp(x*(a1*th1+a2*th2+d))/(1+exp(a1*th1+a2*th2+d))
    f1 =  Deriv::Deriv(f,c("a1"))
    f2 =  Deriv::Deriv(f,c("a2"))
    f3 =  Deriv::Deriv(f,c("d"))

    g = function(pattern, pars) {
      re = spatstat.core::gauss.hermite(function(th1) {
      re = spatstat.core::gauss.hermite(function(th2) {
        prod(f(th1,th2,pars$a[,1],pars$a[,2],pars$d,pattern))
      }, order = 5)},order=5)
      return(re)
    }
    fdot = function(th1,th2,a1,a2,d,x) {
      res = f(th1,th2,a1,a2,d,x)
      temp1 = f1(th1,th2,a1,a2,d,x)
      temp2 = f2(th1,th2,a1,a2,d,x)
      temp3 = f3(th1,th2,a1,a2,d,x)

      # clear NAs? can also do later and not in each function eval?

      a = lapply(1:length(x),function(i) {
        prod(res[-i])* c(temp1[i],temp2[i],temp3[i])
      })
      as.numeric(unlist(a))
    }
    gdot = function(pattern, pars) {
      re = spatstat.core::gauss.hermite(function(th1) {
        re = spatstat.core::gauss.hermite(function(th2) {
          fdot(th1,th2,pars$a[,1],pars$a[,2],pars$d,pattern)
        }, order = 5)},order=5)
      return(re)
    }
    ldot = function(pattern,pars) {
      re=gdot(pattern,pars)/g(pattern,pars)
      return(re)
    }

    funlist=list(f=f,g=g,fdot=fdot,gdot=gdot,ldot=ldot)

  }


# 3PL ---------------------------------------------------------------------

  # Problem ist dass f(x=0) + f(x=1) != 1 ist?

  # 3PL
  if (model=="3PL"&!multi) {

    # f = function(th,a,d,g,x) g + (1-g) * exp(x*(a*th+d))/(1+exp(a*th+d))
    # f = function(th,a,d,g,x) (2 * x - 1) * g + (1-g) * exp(x*(a*th+d))/(1+exp(a*th+d))
    f = function(th,a,d,g,x) (1-x) + (2 * x - 1) * logitinv(g) + (2 * x - 1) * (1-logitinv(g)) /(1+exp(-(a*th+d)))

    f1 =  Deriv::Deriv(f,c("a"))
    f2 =  Deriv::Deriv(f,c("d"))
    f3 =  Deriv::Deriv(f,c("g"))

    g = function(pattern, pars) {
      re = spatstat.core::gauss.hermite(function(th) {
        prod(f(th,pars$a,pars$d,pars$g,pattern))
      }, order = 30)
      return(re)
    }
    fdot = function(th,a,d,g,x) {
      res = f(th,a,d,g,x)
      temp1 = f1(th,a,d,g,x)
      temp2 = f2(th,a,d,g,x)
      temp3 = f3(th,a,d,g,x)

      a = lapply(1:length(x),function(i) {
        prod(res[-i])* c(temp1[i],temp2[i],temp3[i])
      })
      as.numeric(unlist(a))
    }
    gdot = function(pattern, pars) {
      re = spatstat.core::gauss.hermite(function(th) {
        fdot(th,pars$a,pars$d,pars$g,pattern)
      }, order = 30)
      return(re)
    }
    ldot = function(pattern,pars) {
      re=gdot(pattern,pars)/g(pattern,pars)
      return(re)
    }

    funlist=list(f=f,g=g,fdot=fdot,gdot=gdot,ldot=ldot)

  }



# GPCM --------------------------------------------------------------------


  # GPCM
  if (model=="gpcm") {

    f = function(th,a,d1,d2,x) exp(x*((2-x)*(a*th+d1)+(x-1)/2*(2*a*th+d2)))/(1+exp(a*th+d1)+exp(2*a*th+d2))
    f1 =  Deriv::Deriv(f,c("a"))
    f2 =  Deriv::Deriv(f,c("d1"))
    f3 =  Deriv::Deriv(f,c("d2"))

    g = function(pattern, pars) {
      re = spatstat.core::gauss.hermite(function(th) {
        prod(f(th,pars$a,pars$d[,2],pars$d[,3],pattern))
      }, order = 50)
      return(re)
    }

    fdot = function(th,a,d1,d2,x) {
      res = f(th,a,d1,d2,x)
      temp1 = f1(th,a,d1,d2,x)
      temp2 = f2(th,a,d1,d2,x)
      temp3 = f3(th,a,d1,d2,x)

      a = lapply(1:length(x),function(i) {
        prod(res[-i])* c(temp1[i],temp2[i],temp3[i])
      })
      as.numeric(unlist(a))
    }

    gdot = function(pattern, pars) {
      re = spatstat.core::gauss.hermite(function(th) {
        fdot(th,pars$a,pars$d[,2],pars$d[,3],pattern)
      }, order = 50)
      return(re)
    }
    ldot = function(pattern,pars) {
      re=gdot(pattern,pars)/g(pattern,pars)
      return(re)
    }

    funlist=list(f=f,g=g,fdot=fdot,gdot=gdot,ldot=ldot)
  }


  attach(funlist,warn.conflicts = FALSE)




}


