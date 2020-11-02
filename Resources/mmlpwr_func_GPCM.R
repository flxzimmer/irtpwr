
# Basics ------------------------------------------------------------------


# p1  = function(th,a,d1,d2) exp(a*th+d1)
# p2 = function(th,a,d1,d2) exp(2*a*th+d2)
# 
# f.0 = function(th,a,d1,d2) {
#   1/(1+p1(th,a,d1,d2)+p2(th,a,d1,d2))
# }
# f.1 = function(th,a,d1,d2) {
#   p1(th,a,d1,d2)/(1+p1(th,a,d1,d2)+p2(th,a,d1,d2))
# }
# f.2 = function(th,a,d1,d2) {
#   p2(th,a,d1,d2)/(1+p1(th,a,d1,d2)+p2(th,a,d1,d2))
# }
# 
# f = function(th,a,d1,d2,x) {
#   y = (x-1)^2
#   alt=f.0(th,a,d1,d2)^(1-x/2)*f.2(th,a,d1,d2)^(x/2)
#   alt^y * f.1(th,a,d1,d2)^(1-y)
# }

p = function(th,a,d1,d2,x) c(1,exp(a*th+d1),exp(2*a*th+d2))[x+1]
f = function(th,a,d1,d2,x) p(th,a,d1,d2,x)/(1+p(th,a,d1,d2,1)+p(th,a,d1,d2,2))

fderiv = Deriv(f,c("a","d1","d2"))


fvec = function(pattern,th,pars) {
  re = c()
  for (i in 1:length(pars$a)) {
    re[i] = f(th,pars$a[i],pars$d[i,2],pars$d[i,3],pattern[i])
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
    temp1[i] = fderiv(th,pars$a[i],pars$d[i,2],pars$d[i,3],pattern[i])[1]
    temp2[i] = fderiv(th,pars$a[i],pars$d[i,2],pars$d[i,3],pattern[i])[2]
    temp3[i] = fderiv(th,pars$a[i],pars$d[i,2],pars$d[i,3],pattern[i])[3]
    
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


# Thissen Mat -------------------------------------------------------------

ft = function(th,a,d1,d2,x) log(f(th,a,d1,d2,x))
ftd = Deriv(ft,c("a","d1","d2"),nderiv=2)


# ftd = function (th, a, d1, d2, x) 
# {
#   .e1 <- x/2
#   .e2 <- p2(th, a, d1, d2)
#   .e3 <- p1(th, a, d1, d2)
#   .e4 <- a * th
#   .e6 <- 1 + .e3 + .e2
#   .e7 <- (x - 1)^2
#   .e8 <- f.0(th, a, d1, d2)
#   .e9 <- f.2(th, a, d1, d2)
#   .e10 <- 1 - .e1
#   .e11 <- .e8^.e10
#   .e12 <- .e9^.e1
#   .e15 <- exp(2 * .e4 + d2)
#   .e17 <- exp(.e4 + d1)
#   .e18 <- 2 * .e15
#   .e19 <- .e18 + .e17
#   .e20 <- .e11 * .e12
#   .e21 <- .e1 - 1
#   .e22 <- .e8^.e1
#   .e23 <- f.1(th, a, d1, d2)
#   .e24 <- .e9^.e21
#   .e25 <- 1 - .e7
#   .e26 <- .e7 - 1
#   .e27 <- .e20^.e26
#   .e28 <- .e6 * .e22
#   .e29 <- .e20^.e7
#   .e30 <- .e23^.e7
#   .e31 <- .e10 * .e12
#   .e32 <- .e6^2
#   .e33 <- .e2/.e6
#   .e36 <- .e18 - .e19 * .e2/.e6
#   .e37 <- 1 - .e33
#   .e38 <- .e23^.e25
#   .e39 <- .e3/.e6
#   .e40 <- x * .e37
#   .e41 <- x * .e36
#   .e42 <- f(th, a, d1, d2, x)
#   .e43 <- .e10 * .e19
#   .e46 <- 1 - .e39
#   .e47 <- .e17 - .e19 * .e3/.e6
#   .e52 <- .e31/.e22 + x * .e11 * .e24 * .e2/2
#   .e56 <- .e40 * .e11 * .e24/2 - .e31/.e28
#   .e60 <- .e41 * .e11 * .e24/2 - .e43 * .e12/.e28
#   .e61 <- .e19/.e6
#   .e62 <- .e32 * .e22
#   .e64 <- .e27 * .e38 * .e7
#   .e65 <- 2 * .e61
#   .e66 <- .e6 * .e30
#   .e67 <- .e20^(.e7 - 2)
#   .e71 <- .e8^(0.5 * x + 1)
#   .e72 <- .e23^(.e7 + 1)
#   .e73 <- .e9^(.e1 - 2)
#   .e74 <- .e25 * .e46
#   .e75 <- .e52 * .e27
#   .e76 <- .e25 * .e47
#   .e77 <- .e6 * .e72
#   .e78 <- 2 * .e6
#   .e84 <- .e74 * .e29/.e30 - .e75 * .e38 * .e7/.e6
#   .e87 <- .e76 * .e29/.e30 + .e64 * .e60
#   .e92 <- .e64 * .e56 - .e25 * .e29 * .e3/.e66
#   .e93 <- 1 - .e65
#   .e94 <- 2 - .e65
#   .e95 <- 4 * .e15
#   .e96 <- x * .e19
#   .e98 <- .e52 * .e26 * .e67
#   .e99 <- .e6 * .e42
#   .e101 <- .e32 * .e30
#   .e103 <- th^2
#   .e105 <- .e96 * .e12/(2 * (.e32 * .e71))
#   .e106 <- .e52 * .e25
#   .e107 <- .e26 * .e67
#   .e108 <- .e10 * .e36
#   .e109 <- 1 - 2 * (.e15/.e6)
#   .e110 <- 1 - 2 * (.e17/.e6)
#   .e111 <- 2 * .e71
#   .e112 <- .e95 + .e17
#   .e114 <- .e84 * .e87/.e42
#   .e116 <- .e84 * .e92/.e42
#   .e118 <- .e87 * .e92/.e42
#   .e119 <- (.e93 * .e2 + .e18) * .e11
#   .e121 <- .e98 * .e56/.e6
#   .e123 <- .e98 * .e60/.e6
#   .e124 <- .e52 * .e17
#   .e125 <- (.e94 * .e3 + .e17) * .e29
#   .e128 <- .e107 * .e56 * .e60/.e6
#   .e129 <- .e74 * .e27
#   .e130 <- .e93 * .e12
#   .e134 <- .e46 * .e47 * .e29 * .e7/.e77
#   .e136 <- .e46 * .e27 * .e7
#   .e140 <- .e46 * .e29 * .e3 * .e7/.e77
#   .e141 <- .e37 * .e10
#   .e146 <- .e37 * .e36 * .e11 * .e73 * .e21/.e6
#   .e151 <- .e37 * .e11 * .e73 * .e2 * .e21/.e6
#   .e153 <- .e32 * .e42
#   .e155 <- .e6^3 * .e30
#   .e156 <- .e94 * .e12
#   .e161 <- .e36 * .e11 * .e73 * .e2 * .e21/.e6
#   .e163 <- .e47 * .e27 * .e7
#   .e167 <- .e47 * .e29 * .e3 * .e7/.e77
#   .e169 <- 2 * .e32
#   .e170 <- 2 * .e39
#   .e171 <- 2 * .e33
#   .e172 <- 2 * .e12
#   .e174 <- .e41 * .e24/.e78
#   .e176 <- x * .e12/(2 * (.e6 * .e71))
#   c(a = c(a = (((.e107 * (th * .e60/.e6)^2 + (0.5 * (x * (.e11 * 
#                                                             .e73 * (th * .e36/.e6)^2 * .e21 + .e103 * ((.e95 - (.e112 * 
#                                                                                                                   .e2 + 2 * (.e36 * .e19))/.e6) * .e11 - .e108 * .e19/.e62) * 
#                                                             .e24/.e6)) - .e10 * (.e103 * ((.e112 - 2 * (.e19^2/.e6)) * 
#                                                                                             .e12 + .e41 * .e19 * .e24/.e78)/.e62 + x * (-(th * .e19/.e32))^2 * 
#                                                                                    .e12/.e111)) * .e27) * .e38 + .e103 * .e25 * .e47 * .e27 * 
#                   .e60/.e101) * .e7 + .e25 * (.e103 * ((.e17 - (.e112 * 
#                                                                   .e3 + 2 * (.e19 * .e47))/.e6) * .e29 + .e163 * .e60/.e6)/.e66 - 
#                                                 .e29 * (th * .e47/.e6)^2 * .e7/.e72) - (th * .e87/.e6)^2/.e42)/.e42, 
#           d1 = th * (((((0.5 * (x * ((.e43 * .e2/.e62 - .e119) * 
#                                        .e24 - .e161)) - ((.e130 + .e174)/.e22 + .e105) * 
#                            .e10) * .e27 - .e123) * .e38 - .e106 * .e47 * .e27/.e66) * 
#                         .e7 - .e114)/.e6 + (((1 - ((1 - .e61) * .e3 + .e46 * 
#                                                      .e19 + .e17)/.e6) * .e29 + .e136 * .e60/.e6)/.e30 - 
#                                               .e134) * .e25) * .e17/.e99, d2 = th * (((.e128 + 
#                                                                                          (0.5 * (x * (((.e37 * (2 - .e61) - .e36/.e6) * .e11 - 
#                                                                                                          .e141 * .e19/.e62) * .e24 + .e146)) - ((.e156 + 
#                                                                                                                                                    .e174)/.e22 + .e105) * .e10/.e6) * .e27) * .e38 + 
#                                                                                         .e76 * .e27 * .e56/.e66) * .e7 + ((.e167 - (.e125 + 
#                                                                                                                                       .e27 * .e3 * .e7 * .e60/.e6)/.e30) * .e25 - .e118)/.e6) * 
#             .e15/.e99), d1 = c(a = th * (((((((.e96 * .e24 * 
#                                                  .e2/.e169 - .e130)/.e22 - .e105) * .e10 - 0.5 * (x * 
#                                                                                                     ((.e119 + .e108/.e28) * .e24 + .e161))) * .e27 - .e123) * 
#                                              .e38 + .e129 * .e60/.e30) * .e7 - .e114)/.e6 + (((1 - 
#                                                                                                  (.e93 * .e3 + .e18 + 2 * .e17)/.e6) * .e29 - .e52 * .e47 * 
#                                                                                                 .e27 * .e7/.e32)/.e30 - .e134) * .e25) * .e17/.e99, d1 = ((((-(.e124/.e32))^2 * 
#                                                                                                                                                               .e26 * .e67 + (.e10 * (.e17 * (x * .e17 * .e24 * .e2/.e169 - 
#                                                                                                                                                                                                .e110 * .e12)/.e62 - x * (-(.e17/.e32))^2 * .e12/.e111) + 
#                                                                                                                                                                                0.5 * (x * ((-(.e17 * .e2/.e32))^2 * .e11 * .e73 * .e21 + 
#                                                                                                                                                                                              (.e10 * .e17/.e62 - .e110 * .e11) * .e17 * .e24 * 
#                                                                                                                                                                                              .e2/.e32))) * .e27) * .e38 - .e106 * .e46 * .e17^2 * 
#                                                                                                                                                              .e27/.e155) * .e7 + ((.e110 * .e29 - .e124 * .e27 * .e7/.e32) * 
#                                                                                                                                                                                     .e46 * .e17/.e66 - (.e46 * .e17/.e6)^2 * .e29 * .e7/.e72) * 
#                                                                                                                                                             .e25 - (.e84 * .e17/.e6)^2/.e42)/.e42, d2 = ((((((.e172 + 
#                                                                                                                                                                                                                 x * .e24 * .e2/.e78)/.e22 - .e176) * .e10/.e6 + 0.5 * 
#                                                                                                                                                                                                               (x * (((.e171 - 1) * .e11 - .e141/.e28) * .e24 - .e151))) * 
#                                                                                                                                                                                                              .e27 - .e121) * .e38 + .e129 * .e56/.e30) * .e7 + ((.e75 * 
#                                                                                                                                                                                                                                                                    .e3 * .e7/.e32 - (1 - .e170) * .e29)/.e30 + .e140) * 
#                                                                                                                                                                                                            .e25 - .e116) * .e15 * .e17/.e153), d2 = c(a = th * ((((.e163 * 
#                                                                                                                                                                                                                                                                      .e56 - .e125)/.e30 + .e167) * .e25 - .e118)/.e6 + ((.e128 + 
#                                                                                                                                                                                                                                                                                                                            (0.5 * (x * (((2 - (.e94 * .e2 + .e95 + .e17)/.e6) * 
#                                                                                                                                                                                                                                                                                                                                            .e11 - .e108/.e62) * .e24 + .e146)) - ((.e156 + .e40 * 
#                                                                                                                                                                                                                                                                                                                                                                                      .e19 * .e24/.e78)/.e22 + .e105) * .e10/.e6) * .e27) * 
#                                                                                                                                                                                                                                                                                                                           .e38 - .e25 * .e27 * .e3 * .e60/.e101) * .e7) * .e15/.e99, 
#                                                                                                                                                                                                                                                       d1 = ((((((.e172 - .e40 * .e24/2)/.e22 - .e176) * .e10/.e6 + 
#                                                                                                                                                                                                                                                                  0.5 * (x * ((.e10 * .e2/.e62 - (1 - .e171) * .e11) * 
#                                                                                                                                                                                                                                                                                .e24 - .e151))) * .e27 - .e121) * .e38 + .e106 * 
#                                                                                                                                                                                                                                                                .e27 * .e3/.e101) * .e7 + ((.e136 * .e56 + (.e170 - 
#                                                                                                                                                                                                                                                                                                              1) * .e29)/.e30 + .e140) * .e25 - .e116) * .e15 * 
#                                                                                                                                                                                                                                                         .e17/.e153, d2 = (((.e26 * (.e15 * .e56/.e6)^2 * 
#                                                                                                                                                                                                                                                                               .e67 + (0.5 * (x * ((.e109 * .e11 - .e10 * .e15/.e62) * 
#                                                                                                                                                                                                                                                                                                     .e37 * .e15 * .e24/.e6 + (.e37 * .e15/.e6)^2 * .e11 * 
#                                                                                                                                                                                                                                                                                                     .e73 * .e21)) - ((.e109 * .e12 + .e40 * .e15 * .e24/.e78) * 
#                                                                                                                                                                                                                                                                                                                        .e15/.e62 + x * (-(.e15/.e32))^2 * .e12/.e111) * 
#                                                                                                                                                                                                                                                                                         .e10) * .e27) * .e38 - .e25 * .e15^2 * .e27 * .e3 * 
#                                                                                                                                                                                                                                                                              .e56/.e155) * .e7 - (((-(.e15 * .e3/.e32))^2 * .e29 * 
#                                                                                                                                                                                                                                                                                                      .e7/.e72 + (.e109 * .e29 + .e15 * .e27 * .e7 * .e56/.e6) * 
#                                                                                                                                                                                                                                                                                                      .e15 * .e3/.e101) * .e25 + (.e92 * .e15/.e6)^2/.e42))/.e42))
# }


tmat = function(pars) {
  res = list()
  for (i in 1:length(pars$a)) {
    res[[i]] = gauss.hermite(function(th) {
      likelihoods = f(th,pars$a[i],pars$d[i,2],pars$d[i,3],0:2)
      sderiv = Reduce('+',lapply(0:2,function(x) likelihoods[x+1]*ftd(th,pars$a[i],pars$d[i,2],pars$d[i,3],x)))
      -matrix(sderiv,3,3) 
    },order=30)
  }
  as.matrix(do.call(bdiag,res))
}


# L Optimizer ------------------------------------------------------------



maxl2preload= function(pars) {
  
  n.items = length(pars$a)
  n.kat = max(ncol(pars$d),2)
  patterns = as.matrix(expand.grid(lapply(1:n.items,function(x) 0:(n.kat-1))))
  
  pre = c()
  for (i in 1:nrow(patterns)) {
    pre[i] = g(patterns[i,],pars)
  }
  
  return(pre)
  
}



maxl2 = function(x,pars,pre) {
  
  n.items = length(pars$a)
  n.kat = max(ncol(pars$d),2)
  patterns = as.matrix(expand.grid(lapply(1:n.items,function(x) 0:(n.kat-1))))
  x = list(a=rep(x[1],n.items),d=matrix(c(rep(0,n.items),x[2:length(x)]),ncol=ncol(pars$d)))
  
  res  = c()
  for (i in 1:nrow(patterns)) {
    px = pre[i]
    qx = g(patterns[i,],x)
    res[i] =  {px*log(qx)}
  }
  
  re = -sum(res)
  
}






