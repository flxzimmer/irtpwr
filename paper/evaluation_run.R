###########################################################################
# Run simulation on the server
###########################################################################

# prep --------------------------------------------------------------------

library(mmlpwrpackage)
library(sn)
library(rlist)

set.seed(456)

load(file= "C:/Users/felix/Google Drive/4 irt/mmlpwr/simsetup.Rdata")

# Distribution: Calculate Stats ---------------------------------------------------------

sim1.split = list()
setsize = 50
for (i in sim1) {
  for (j in 1:i$simruns/setsize) {
    a = list(condition = i$condition,datasets = i$datasets[((j-1)*setsize+1):(j*setsize)])
    a %>% list.append(sim1.split,.)
  }
}

cl <- makeCluster(10)
clusterExport(cl, objects(), envir=environment())

res1.split =  parLapply(cl,X = sim1.split,fun = function(x){

  load.libraries()

  fits = list()
  obs = list()

  for (i in x$datasets) {
    # fit the mirt model
    mml.fit(data = x$data,hyp = x$hyp) %>% list.append(fits,.)
    # calculate statistics
    c(wald_obs(fitted),lr_obs(fitted),score_obs(fitted)) %>% list.append(obs,.)
  }

  re = list(condition = x$condition,fits=fits,obs = obs)
  return(re)

})

# glue all with the same condition together again
res1 = lapply(sim1,function(x) {x$fits=list();x$obs = list();x = x[names(x)!="datasets"]})
for (i in res1.split) {
  i$fits %>% list.append(x[x$condition==i$condition]$fits,.)
  i$obs %>% list.append(x[x$condition==i$condition]$obs,.)

}

save(res1,file= "results1.Rdata")
stopCluster(cl)



# Distribution: Calculate ncps ----------------------------------------------------------

cl <- makeCluster(10)
clusterExport(cl, objects(), envir=environment())

res1ncp =  parLapply(cl,X = sim1ncp,fun = function(x){

  load.libraries()

  analytical = NULL

  # calculate ncps for the condition
  if (x$n.items==10) {
    analytical = calculate_ncps(x$hyp,simbased=FALSE,n=1)
  }
  simbased = calculate_ncps(x$hyp,simbased=TRUE,n=1,n.pers = NULL)

  # Sim Condition without the datasets
  condition = x[names(x)!="datasets"]

  re = list(condition = condition,analytical = analytical,simbased=simbased)
  return(re)
  })
save(res1ncp,file= "results1ncp.Rdata")
stopCluster(cl)



# Power - Calculate ncps  -----------------------------------------

cl <- makeCluster(10)
clusterExport(cl, objects(), envir=environment())

res2ncp =  parLapply(cl,X = sim2ncp,fun = function(x){

  load.libraries()

  analytical = NULL

  # calculate ncps for the condition
  if (x$n.items==10) {
    analytical = calculate_ncps(x$hyp,simbased=FALSE,n=1)
  }
  simbased = calculate_ncps(x$hyp,simbased=TRUE,n=1,n.pers = NULL)

  # Sim Condition without the datasets
  condition = x[names(x)!="datasets"]

  re = list(condition = condition,analytical = analytical,simbased=simbased)
  return(re)

})

save(res2ncp,file= "results2ncp.Rdata")
stopCluster(cl)



# Power - Calculate sample size -------------------------------------

# add the number of persons and the datasets to the simulation conditions

powerset=c(.2,.4,.6,.8)

# calculate the sample sizes according to the ncps
sample_sizes = lapply(res2ncp, function(x) {

  if (!is.null(x$analytical)) {
    ncps = x$analytical
  } else {
    ncps = x$simbased
  }

  re = sapply(powerset,function(pw) {ssize(hyp=x$condition$hyp,ncp=ncps,power = pw,alpha=.05) %>% mean()})
  return(re)
})

#add the sample sizes

res2 = list()

for (i in 1:length(sim2ncp)) {

  ss = sample_sizes[i]
  a = sim2ncp[i]
  for (x in ss){
    b = a
    b$n.pers = x
    res2 = c(res2,b)
  }
}

#add the datasets

sim2 = lapply(sim2,function(x){

  x$datasets = list()
  for (i in 1:x$simruns) {
    x$datasets[[i]] = setup.data(hyp = x$hyp, n = x$n.pers)
  }
  return(x)

})


# Power - Calculate stats --------------------------------------------------


cl <- makeCluster(10)
clusterExport(cl, objects(), envir=environment())


res2 =  parLapply(cl,X = sim2,fun = function(x){

  load.libraries()

  # fit the mirt model
  fitted <- mml.fit(data = x$data,hyp = x$hyp)
  # calculate statistics
  stats_obs <- c(wald_obs(fitted),lr_obs(fitted),score_obs(fitted))

  # Sim Condition without the datasets
  condition = x[names(x)!="datasets"]

  re = list(condition = condition,fitted=fitted,stats_obs = stats_obs)
  return(re)

})
save(res2,file= "results2.Rdata")
stopCluster(cl)



# Misspecification: Calculate ncp --------------------------------------------------------

# Builds on the Power condition 1PLvs2PL, 10 Items

# The ncps are therefore included in "results2ncp.Rdata"

indices = sapply(res2ncp,function(x) {
  re = x$condition$type = "1PLvs2PL" & x$condition$n.items = 10
  return(re)
})

res3ncp = res2ncp[indices]

save(res3ncp,file= "results3ncp.Rdata")


# add sample sizes to the data.

# Setup the simulation condition - based on sim2. datasets have to be generated anew for the non-normal distributions

indices = sapply(sim2,function(x) {
  re = x$condition$type = "1PLvs2PL" & x$condition$n.items = 10
  return(re)
})

sim3 = list()

for (x in sim2[indices]) {
  a = x
  b = x

  a$dist.fun = unifdist
  b$dist.fun = skeweddist

  sim3 = c(sim3,a,b)
}


sim3 = lapply(sim3,function(x){

  x$datasets = list()
  for (i in 1:x$simruns) {
    x$datasets[[i]] = setup.data(hyp = x$hyp, n = x$n.pers,dist.type=x$dist.fun)
  }
  return(x)

})



# Misspecification: Calculate stats -----------------------------------------


cl <- makeCluster(10)
clusterExport(cl, objects(), envir=environment())


res3 =  parLapply(cl,X = sim3,fun = function(x){

  load.libraries()

  # fit the mirt model
  fitted <- mml.fit(data = x$data,hyp = x$hyp)
  # calculate statistics
  stats_obs <- c(wald_obs(fitted),lr_obs(fitted),score_obs(fitted))

  # Sim Condition without the datasets
  condition = x[names(x)!="datasets"]

  re = list(condition = condition,fitted=fitted,stats_obs = stats_obs)
  return(re)

})
save(res3,file= "results3.Rdata")
stopCluster(cl)



# save updated sim files --------------------------------------------------

save(sim2,sim3,file= "C:/Users/felix/Google Drive/4 irt/mmlpwr/simsetup2.Rdata")

