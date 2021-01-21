
library(mmlpwrpackage)
library(mirt)
library(dplyr)
library(parallel)

set.seed(1234)
altpars <- list(
  a = rlnorm(10,sdlog = .1),
  d = rnorm(10)
)
hyp <- setup_hypothesis(type = "1PLvs2PL", altpars = altpars)

cond = lapply(1:100,function(x) list(hyp=hyp))

cl <- makeCluster(120)

res =  parLapply(cl,X = cond,fun = function(x){

  library(mmlpwrpackage)
  library(mirt)
  library(dplyr)
  library(parallel)

  k = 20
  # k = 2

  ptm <- proc.time()
  ncps1 <- calculate_ncps(hyp=x$hyp,simbased = TRUE,n.pers = 1000*k,runs=1)
  t1 = proc.time() - ptm
  ptm <- proc.time()
  ncps2 <- calculate_ncps(hyp=x$hyp,simbased = TRUE,n.pers = 1000,runs = k)
  t2 = proc.time() - ptm

  re = list(ncps1,ncps2,t1,t2)

  return(re)
})


save(res,file = "restemp123.RData")
stopCluster(cl)
