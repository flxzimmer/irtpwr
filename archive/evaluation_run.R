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


#load(file= "simsetup.Rdata")
load(file= "C:/Users/felix/Google Drive/4 irt/mmlpwr/simsetup.Rdata")


# Distribution: Calculate Stats ---------------------------------------------------------

sim1.split = list()
setsize = 20
for (i in sim1) {
  for (j in 1:(i$simruns/setsize)) {
    a = i
    a$datasets = i$datasets[((j-1)*setsize+1):(j*setsize)]
    sim1.split = a %>% list.append(sim1.split,.)
  }
}

cl <- makeCluster(120)
#clusterExport(cl, "load.libraries", envir=environment())

res1.split =  parLapply(cl,X = sim1.split,fun = function(x){

  #  load.libraries()

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

  # fits = list()
  obs = list()

  for (i in x$datasets) {
    # fit the mirt model
    fitted = mml.fit(data = i,hyp = x$hyp)
    # fits = fitted %>% list.append(fits,.)
    # calculate statistics
    obs = c(wald_obs(fitted),lr_obs(fitted),score_obs(fitted)) %>% list.append(obs,.)
  }
  # re = list(condition = x$condition,fits=fits,obs = obs)
  re = list(condition = x$condition,obs = obs)
  return(re)
})

save(res1.split,file= "results1split.Rdata")
stopCluster(cl)


# glue all with the same condition together again
res1 = lapply(sim1,function(x) {x$obs = list();x = x[names(x)!="datasets"];return(x)})
for (i in res1.split) {
  index = which(sapply(sim1, function(x) digest(x$condition)==digest(i$condition))) %>%as.numeric()
  res1[[index]]$obs = c(res1[[index]]$obs,i$obs)
}
save(res1,file= "results1.Rdata")



# Distribution: Calculate ncps ----------------------------------------------------------

cl <- makeCluster(30)
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
  if (x$n.items==10) {
    analytical = calculate_ncps(x$hyp,simbased=FALSE,n=1)
  } else {
    analytical = NULL
  }
  simbased = calculate_ncps(x$hyp,simbased=TRUE,n=1,n.pers = NULL)

  # # Sim Condition without the datasets
  # condition = x[names(x)!="datasets"]

  re = list(condition = condition,analytical = analytical,simbased=simbased)
  return(re)
  })
save(res1ncp,file= "results1ncp.Rdata")
stopCluster(cl)



# Misspecification --------------------------------------------------------

