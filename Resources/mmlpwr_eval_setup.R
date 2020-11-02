


# prep --------------------------------------------------------------------

source('mmlpwr_func5.R')

# setup basic  ------------------------------------------------------------


sim1 = expand.grid(preset=c("1PLvs2PL","DIF2PL"),n.items=10,esize=c("no","small","large"),n.pers=c(500,1000,3000))
sim1 = sim1[order(sim1$preset,sim1$n.items,sim1$esize,sim1$n.pers),]
sim1 = sim1  %>% split(., seq(nrow(.))) %>% lapply(.,as.list)


sim1 = lapply(sim1,function(x) {
  x$null.hypothesis =setup.null.hypothesis(n.items=x$n.items,preset=x$preset)
  x$alternative.hypothesis = setup.alternative.hypothesis(null.hypothesis = x$null.hypothesis,esize=x$esize)
  return(x)
})


sim2 = expand.grid(preset=c("1PLvs2PL","DIF2PL"),n.items=50,esize=c("no","small","large"),n.pers=c(500,1000,3000))
sim2 = sim2[order(sim2$preset,sim2$n.items,sim2$esize,sim2$n.pers),]
sim2 = sim2  %>% split(., seq(nrow(.))) %>% lapply(.,as.list)


sim2 = lapply(sim2,function(x) {
  x$null.hypothesis =setup.null.hypothesis(n.items=x$n.items,preset=x$preset)
  x$alternative.hypothesis = setup.alternative.hypothesis(null.hypothesis = x$null.hypothesis,esize=x$esize)
  return(x)
})



# setup power -------------------------------------------------------------

null.hypothesis = setup.null.hypothesis(n.items=10,preset="1PLvs2PL")
alternative.hypothesis = setup.alternative.hypothesis(null.hypothesis = null.hypothesis,esize="large")
mean.ncp = c(ncp.wald(null.hypothesis,alternative.hypothesis),
             ncp.lr(null.hypothesis,alternative.hypothesis),
             ncp.score(null.hypothesis,alternative.hypothesis)) %>% mean()
powerset = c(.2,.4,.6,.8)
sample.sizes1 = sapply(powerset, function(x) ssize(df=null.hypothesis$df,ncp=mean.ncp,alpha=.05,power=x))


null.hypothesis = setup.null.hypothesis(n.items=10,preset="DIF2PL")
alternative.hypothesis = setup.alternative.hypothesis(null.hypothesis = null.hypothesis,esize="large")
mean.ncp = c(ncp.wald(null.hypothesis,alternative.hypothesis),
             ncp.lr(null.hypothesis,alternative.hypothesis),
             ncp.score(null.hypothesis,alternative.hypothesis)) %>% mean()
powerset = c(.2,.4,.6,.8)
sample.sizes2 = sapply(powerset, function(x) ssize(df=null.hypothesis$df,ncp=mean.ncp,alpha=.05,power=x))


sim3 = rbind(expand.grid(preset=c("1PLvs2PL"),n.items=10,esize=c("large"),n.pers=sample.sizes1),expand.grid(preset=c("DIF2PL"),n.items=10,esize=c("large"),n.pers=sample.sizes2))
sim3 = sim3[order(sim3$preset,sim3$n.items,sim3$esize,sim3$n.pers),]
sim3 = sim3  %>% split(., seq(nrow(.))) %>% lapply(.,as.list)

sim3 = lapply(sim3,function(x) {
  x$null.hypothesis =setup.null.hypothesis(n.items=x$n.items,preset=x$preset)
  x$alternative.hypothesis = setup.alternative.hypothesis(null.hypothesis = x$null.hypothesis,esize=x$esize)
  return(x)
})


# setup power 50 items -------------------------------------------------------------

null.hypothesis = setup.null.hypothesis(n.items=50,preset="1PLvs2PL")
alternative.hypothesis = setup.alternative.hypothesis(null.hypothesis = null.hypothesis,esize="xsmall")
mean.ncp = ncp.wald(null.hypothesis,alternative.hypothesis ,method="mirtOakessim")
powerset = c(.2,.4,.6,.8)
sample.sizes1 = sapply(powerset, function(x) ssize(df=null.hypothesis$df,ncp=mean.ncp,alpha=.05,power=x))

null.hypothesis = setup.null.hypothesis(n.items=50,preset="DIF2PL")
alternative.hypothesis = setup.alternative.hypothesis(null.hypothesis = null.hypothesis,esize="small")
mean.ncp = ncp.wald(null.hypothesis,alternative.hypothesis ,method="mirtOakessim")
powerset = c(.2,.4,.6,.8)
sample.sizes2 = sapply(powerset, function(x) ssize(df=null.hypothesis$df,ncp=mean.ncp,alpha=.05,power=x))


sim4 = rbind(expand.grid(preset=c("1PLvs2PL"),n.items=50,esize=c("xsmall"),n.pers=sample.sizes1),expand.grid(preset=c("DIF2PL"),n.items=50,esize=c("small"),n.pers=sample.sizes2))
sim4 = sim4[order(sim4$preset,sim4$n.items,sim4$esize,sim4$n.pers),]
sim4 = sim4  %>% split(., seq(nrow(.))) %>% lapply(.,as.list)
sim4 = lapply(sim4,function(x) {
  x$null.hypothesis =setup.null.hypothesis(n.items=x$n.items,preset=x$preset)
  x$alternative.hypothesis = setup.alternative.hypothesis(null.hypothesis = x$null.hypothesis,esize=x$esize)
  return(x)
})


sim = c(sim1,sim2,sim3,sim4)
save(sim,file= "sim.Rdata")




# # prep --------------------------------------------------------------------
# 
# source('mmlpwr_func5.R')
# 
# # setup basic  ------------------------------------------------------------
# 
# 
# sim1 = expand.grid(preset=c("1PLvs2PL","DIF2PL"),n.items=5,esize=c("no","small","large"),n.pers=c(500,1000,3000))
# sim1 = sim1[order(sim1$preset,sim1$n.items,sim1$esize,sim1$n.pers),] 
# sim1 = sim1  %>% split(., seq(nrow(.))) %>% lapply(.,as.list)
# 
# 
# sim1 = lapply(sim1,function(x) {
#   x$null.hypothesis =setup.null.hypothesis(n.items=x$n.items,preset=x$preset)
#   x$alternative.hypothesis = setup.alternative.hypothesis(null.hypothesis = x$null.hypothesis,esize=x$esize)
#   return(x)
#   })
# 
# 
# sim2 = expand.grid(preset=c("1PLvs2PL","DIF2PL"),n.items=10,esize=c("no","small","large"),n.pers=c(500,1000,3000))
# sim2 = sim2[order(sim2$preset,sim2$n.items,sim2$esize,sim2$n.pers),]
# sim2 = sim2  %>% split(., seq(nrow(.))) %>% lapply(.,as.list)
# 
# 
# sim2 = lapply(sim2,function(x) {
#   x$null.hypothesis =setup.null.hypothesis(n.items=x$n.items,preset=x$preset)
#   x$alternative.hypothesis = setup.alternative.hypothesis(null.hypothesis = x$null.hypothesis,esize=x$esize)
#   return(x)
# })
# 
# 
# 
# # setup power -------------------------------------------------------------
# 
# null.hypothesis = setup.null.hypothesis(n.items=5,preset="1PLvs2PL")
# alternative.hypothesis = setup.alternative.hypothesis(null.hypothesis = null.hypothesis,esize="large")
# mean.ncp = c(ncp.wald(null.hypothesis,alternative.hypothesis),
#              ncp.lr(null.hypothesis,alternative.hypothesis),
#              ncp.score(null.hypothesis,alternative.hypothesis)) %>% mean()
# powerset = c(.2,.4,.6,.8)
# sample.sizes1 = sapply(powerset, function(x) ssize(df=null.hypothesis$df,ncp=mean.ncp,alpha=.05,power=x)) 
# 
# 
# null.hypothesis = setup.null.hypothesis(n.items=5,preset="DIF2PL")
# alternative.hypothesis = setup.alternative.hypothesis(null.hypothesis = null.hypothesis,esize="large")
# mean.ncp = c(ncp.wald(null.hypothesis,alternative.hypothesis),
#              ncp.lr(null.hypothesis,alternative.hypothesis),
#              ncp.score(null.hypothesis,alternative.hypothesis)) %>% mean()
# powerset = c(.2,.4,.6,.8)
# sample.sizes2 = sapply(powerset, function(x) ssize(df=null.hypothesis$df,ncp=mean.ncp,alpha=.05,power=x)) 
# 
# 
# sim3 = rbind(expand.grid(preset=c("1PLvs2PL"),n.items=5,esize=c("large"),n.pers=sample.sizes1),expand.grid(preset=c("DIF2PL"),n.items=5,esize=c("large"),n.pers=sample.sizes2))
# sim3 = sim3[order(sim3$preset,sim3$n.items,sim3$esize,sim3$n.pers),] 
# sim3 = sim3  %>% split(., seq(nrow(.))) %>% lapply(.,as.list)
# 
# sim3 = lapply(sim3,function(x) {
#   x$null.hypothesis =setup.null.hypothesis(n.items=x$n.items,preset=x$preset)
#   x$alternative.hypothesis = setup.alternative.hypothesis(null.hypothesis = x$null.hypothesis,esize=x$esize)
#   return(x)
# })
# 
# 
# # setup power 50 items -------------------------------------------------------------
# 
# null.hypothesis = setup.null.hypothesis(n.items=10,preset="1PLvs2PL")
# # alternative.hypothesis = setup.alternative.hypothesis(null.hypothesis = null.hypothesis,esize="xsmall")
# alternative.hypothesis = setup.alternative.hypothesis(null.hypothesis = null.hypothesis,esize="large")
# mean.ncp = ncp.sim(null.hypothesis,alternative.hypothesis) %>% mean()
# powerset = c(.2,.4,.6,.8)
# sample.sizes1 = sapply(powerset, function(x) ssize(df=null.hypothesis$df,ncp=mean.ncp,alpha=.05,power=x))
# 
# null.hypothesis = setup.null.hypothesis(n.items=10,preset="DIF2PL")
# alternative.hypothesis = setup.alternative.hypothesis(null.hypothesis = null.hypothesis,esize="small")
# mean.ncp = ncp.sim(null.hypothesis,alternative.hypothesis) %>% mean()
# powerset = c(.2,.4,.6,.8)
# sample.sizes2 = sapply(powerset, function(x) ssize(df=null.hypothesis$df,ncp=mean.ncp,alpha=.05,power=x))
# 
# 
# sim4 = rbind(expand.grid(preset=c("1PLvs2PL"),n.items=10,esize=c("small"),n.pers=sample.sizes1),expand.grid(preset=c("DIF2PL"),n.items=10,esize=c("small"),n.pers=sample.sizes2))
# sim4 = sim4[order(sim4$preset,sim4$n.items,sim4$esize,sim4$n.pers),]
# sim4 = sim4  %>% split(., seq(nrow(.))) %>% lapply(.,as.list)
# sim4 = lapply(sim4,function(x) {
#   x$null.hypothesis =setup.null.hypothesis(n.items=x$n.items,preset=x$preset)
#   x$alternative.hypothesis = setup.alternative.hypothesis(null.hypothesis = x$null.hypothesis,esize=x$esize)
#   return(x)
# })
# 
# 
# sim = c(sim1,sim2,sim3,sim4)
# save(sim,file= "sim.Rdata")
