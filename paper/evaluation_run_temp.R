###########################################################################
# Run simulation on the server
###########################################################################

# prep --------------------------------------------------------------------

library(mmlpwrpackage)
library(sn)
library(rlist)
library(digest)
library(parallel)
library(dplyr)


load(file= "simsetupa.Rdata")
# load(file= "C:/Users/felix/Google Drive/4 irt/mmlpwr/simsetup.Rdata")


# Distribution: Calculate Stats ---------------------------------------------------------

sim1.split = list()
setsize = 5
for (i in sim1) {
  for (j in 1:(i$simruns/setsize)) {
    a = i
    a$datasets = i$datasets[((j-1)*setsize+1):(j*setsize)]
    sim1.split = a %>% list.append(sim1.split,.)
  }
}

cl <- makeCluster(120)

res1.split =  parLapply(cl,X = sim1.split,fun = function(x){

  library(mmlpwrpackage)
  library(rlist)
  library(sn)
  library(mirt)
  library(dplyr)
  library(spatstat)
  library(Deriv)
  library(digest)
  library(Matrix)
  library(MASS)
  library(sn)

  obs = list()

  for (i in x$datasets) {
    # fit the mirt model
    fitted = mml.fit(data = i,hyp = x$hyp)
    # calculate statistics
    obs = c(wald_obs(fitted),lr_obs(fitted),score_obs(fitted)) %>% list.append(obs,.)
  }
  re = list(condition = x$condition,obs = obs)
  return(re)
})

save(res1.split,file= "results1splita.Rdata")
stopCluster(cl)


# glue all with the same condition together again
res1 = lapply(sim1,function(x) {x$obs = list();x = x[names(x)!="datasets"];return(x)})
for (i in res1.split) {
  index = which(sapply(sim1, function(x) digest(x$condition)==digest(i$condition))) %>%as.numeric()
  res1[[index]]$obs = c(res1[[index]]$obs,i$obs)
}
save(res1,file= "results1a.Rdata")


# Distribution: Calculate ncps ----------------------------------------------------------

cl <- makeCluster(10)
#clusterExport(cl, objects(), envir=environment())

res1ncp =  parLapply(cl,X = sim1ncp,fun = function(x){

  library(mmlpwrpackage)
  library(rlist)
  library(sn)
  library(mirt)
  library(dplyr)
  library(spatstat)
  library(Deriv)
  library(digest)
  library(Matrix)
  library(MASS)
  library(sn)

  # calculate ncps for the condition
  if (isFALSE(x$condition$simbased.only)) {
    analytical = calculate_ncps(x$hyp,simbased=FALSE)
  } else {
    analytical = NULL
  }
  simbased = calculate_ncps(x$hyp,simbased=TRUE,n.pers = x$condition$simbased.pers,runs =x$condition$simbased.runs)

  re = list(condition = x$condition,analytical = analytical,simbased=simbased)
  return(re)
})
save(res1ncp,file= "results1ncpa.Rdata")
stopCluster(cl)
