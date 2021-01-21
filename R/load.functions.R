

#' Load IRT model functions and derivatives
#'
#' @param model Desired Model (2PL only at the moment)
#'
#' @return
#' @export
#'
#' @examples
load.functions = function (model) {

  # 2PL
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

    funlist=list(f=f,fvec=fvec,g=g,fdot=fdot,gdot=gdot,ldot=ldot)

  }


  attach(funlist,warn.conflicts = FALSE)

}



