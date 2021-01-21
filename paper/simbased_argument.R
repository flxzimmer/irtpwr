
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


# analyse -----------------------------------------------------------------

load(file= "C:/Users/admin/Google Drive/4 irt/mmlpwr/restemp123.Rdata")

# time comparisons
ts1=sapply(res,function(x) x[[3]]["elapsed"])
ts2=sapply(res,function(x) x[[4]]["elapsed"])

mean(ts1);mean(ts2)



# same sd?
ncps1=lapply(res,function(x) x[[1]]) %>% do.call(rbind,.)
ncps2=lapply(res,function(x) x[[2]]) %>% do.call(rbind,.)
sds1 = ncps1 %>% apply(.,2,sd)
sds2 = ncps2 %>% apply(.,2,sd)

sds1/sds2

# same mean?
means1 = ncps1 %>% apply(.,2,mean)
means2 = ncps2 %>% apply(.,2,mean)

means1;means2


# comparison to analytical ncp
ncpsana <- calculate_ncps(hyp=hyp,simbased = FALSE)

# can the get_ncp function be improved so that it does not come out as 0 that often?


# old  --------------------------------------------------------------------


# altpars <- list(
#   a = rlnorm(10,sdlog = .1),
#   d = rnorm(10)
# )
#
# hyp <- setup_hypothesis(type = "1PLvs2PL", altpars = altpars)
#
# k = 100
#
# res = list()
#
#
# for (i in 1:20) {
#   print(paste("run",i))
#   ptm <- proc.time()
#   ncps1 <- calculate_ncps(hyp=hyp,simbased = TRUE,n.pers = 1000*k,runs=1)
#   t1 = proc.time() - ptm
#   ptm <- proc.time()
#   ncps2 <- calculate_ncps(hyp=hyp,simbased = TRUE,n.pers = 1000,runs = k)
#   t2 = proc.time() - ptm
#
#
#   res[[i]] = list(ncps1,ncps2,t1,t2)
# }
#
# save(res,file = "restemp123.RData")
#
# # time comparisons -> boosted n is faster
# ts1=sapply(res,function(x) x[[3]]["elapsed"])
# ts2=sapply(res,function(x) x[[4]]["elapsed"])
#
# mean(ts1);mean(ts2)
#
#
# # same sd?
# ncps1=lapply(res,function(x) x[[1]]) %>% do.call(rbind,.)
# ncps2=lapply(res,function(x) x[[2]]) %>% do.call(rbind,.)
# sds1 = ncps1 %>% apply(.,2,sd)
# sds2 = ncps2 %>% apply(.,2,sd)
#
# sds1/sds2
#
# # same mean?
#
#


