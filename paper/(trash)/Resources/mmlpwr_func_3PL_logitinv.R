
p = function(th,a,d,g) logitinv(g)+(1-logitinv(g))/(1+exp(-(a*th+d)))
#f = function(th,a,d,g,x) p(th,a,d,g)^x * (1-p(th,a,d,g))^(1-x)
# f = function(th,a,d,g,x) {1-c(p(th,a,d,g),p(th,a,d,g))[x+1]}
f = function(th,a,d,g,x) {c(1-p(th,a,d,g),p(th,a,d,g))[x+1]}

f1 =  Deriv(f,c("a"))
f2 =  Deriv(f,c("d"))
f3 =  Deriv(f,c("g"))

fvec = function(pattern,th,pars) {
  re = c()
  for (i in 1:length(pars$a)) {
    re[i] = f(th,pars$a[i],pars$d[i],pars$g[i],pattern[i])
  }
  return(re)
}

g = function(pattern, pars) {
  re = gauss.hermite(function(th) {
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
    temp3 = res
    temp1[i] = f1(th,pars$a[i],pars$d[i],pars$g[i],pattern[i])
    temp2[i] = f2(th,pars$a[i],pars$d[i],pars$g[i],pattern[i])
    temp3[i] = f3(th,pars$a[i],pars$d[i],pars$g[i],pattern[i])
    re = c(re,c(prod(temp1),prod(temp2),prod(temp3)))
  }
  return(re)
}

gdot = function(pattern, pars) {
  re = gauss.hermite(function(th) {
    fdot(pattern, th, pars)
  }, order = 30)
  return(re)
}


ldot = function(pattern,pars) {
  re=gdot(pattern,pars)/g(pattern,pars)
  return(re)
}



# L Optimizer ------------------------------------------------------------



maxl = function(x,pars,i,restriction) {
  if (restriction=="g0.2") {gval= .2} else
    if (restriction=="g0") {gval= 0}
  p = function(th,a,d,g) g+(1-g)/(1+exp(-(a*th+d)))
  f = function(th,a,d,g,x) p(th,a,d,g)^x * (1-p(th,a,d,g))^(1-x)
  px = function(th) {f(th,pars$a[i],pars$d[i],pars$g[i],1)}
  qx = function(th) {f(th,x[1],x[2],gval,1)}
  kl = function(th) {px(th)*log(qx(th))+(1-px(th))*log(1-qx(th)) }
  re = -gauss.hermite(kl,order=20)
}

