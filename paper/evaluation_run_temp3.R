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


load(file= "simsetup.Rdata")
# load(file= "C:/Users/felix/Google Drive/4 irt/mmlpwr/simsetup.Rdata")

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
save(res1ncp,file= "results1ncp.Rdata")
stopCluster(cl)


# Misspecification: Distribution --------------------------------------------------------

sim2.split = list()
setsize = 5
for (i in sim2) {
  for (j in 1:(i$simruns/setsize)) {
    a = i
    a$datasets = i$datasets[((j-1)*setsize+1):(j*setsize)]
    sim2.split = a %>% list.append(sim2.split,.)
  }
}

cl <- makeCluster(60)

res2.split =  parLapply(cl,X = sim2.split,fun = function(x){

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

save(res2.split,file= "results2split.Rdata")
stopCluster(cl)


# glue all with the same condition together again
res2 = lapply(sim2,function(x) {x$obs = list();x = x[names(x)!="datasets"];return(x)})
for (i in res2.split) {
  index = which(sapply(sim2, function(x) digest(x$condition)==digest(i$condition))) %>%as.numeric()
  res2[[index]]$obs = c(res2[[index]]$obs,i$obs)
}
save(res2,file= "results2.Rdata")

