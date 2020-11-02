


# prep --------------------------------------------------------------------

source('mmlpwr_func5.R')

load(file= "sim.Rdata")

# sim3 =sim[37:44]
sim3 =sim[37:40]



# exec --------------------------------------------------------------------


cl <- makeCluster(10)
clusterExport(cl, objects(), envir=environment())

# misspecification

res1 =  parLapply(cl,X = sim3,fun = function(x){
  list(stats=stats.simfunc(x,dist.type="unif"))
})
save(res1,file= "results1.Rdata")

res2 =  parLapply(cl,X = sim3,fun = function(x){
  list(stats=stats.simfunc(x,dist.type="skewed"))
})
save(res2,file= "results2.Rdata")


stopCluster(cl)


