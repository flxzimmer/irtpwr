# 
# 
# 
# # prep --------------------------------------------------------------------
# 
# source('mmlpwr_func5.R')
# 
# load(file= "sim - Copy (4).Rdata")
# 
# # exec --------------------------------------------------------------------
# 
# 
# cl <- makeCluster(60)
# clusterExport(cl, objects(), envir=environment())
# 
# 
# res =  parLapply(cl,X = sim,fun = function(x){
#   list(stats=stats.simfunc(x),ncps=ncp.simfunc(x))
#             })
# save(res,file= "results_addon.Rdata")
# stopCluster(cl)


source('mmlpwr_func5.R')

load(file= "sim - Copy (5).Rdata")
sim = sim[45:48]


# exec --------------------------------------------------------------------


cl <- makeCluster(10)
clusterExport(cl, objects(), envir=environment())


res =  parLapply(cl,X = sim,fun = function(x){
  list(stats=stats.simfunc(x),ncps=ncp.simfunc(x))
})
save(res,file= "results_addon.Rdata")
stopCluster(cl)

