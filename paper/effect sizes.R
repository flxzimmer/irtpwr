

library(mirt)
library(dplyr)
library(mmlpwrpackage)
library(TAM)


# load simulation data

load(file= "C:/Users/admin/Google Drive/4 irt/mmlpwr/simsetup.Rdata")

# 1PL vs 2PL: Calculate sd of a parameters

df = sim1[sapply(sim1,function(x) x$type=="1PLvs2PL")] # 1PLvs2PL conditions only
rm(sim1)
rm(sim2)

sd.apar = sapply(df, function(x) {

  re = c()
  for (i in 1:500) {
    pars = mirt(as.data.frame(x$datasets[[i]]),1,verbose = FALSE) %>% coef_short(.)
    re=c(re,sd(pars$a))
  }
  return(mean(re))
})


# Load datasets
dat6 <- expand.table(LSAT6)
dat7 <- expand.table(LSAT7)
dat8 <- expand.table(deAyala)
data(data.geiser)
dat9 <- data.geiser

data(data.numeracy)
data(data.janssen2)

dat10 <-data.numeracy

# Estimate pars
pars6 = mirt(dat6,1,verbose = FALSE) %>% coef_short(.)
pars7 = mirt(dat7,1,verbose = FALSE) %>% coef_short(.)
pars8 = mirt(dat8,1,verbose = FALSE) %>% coef_short(.)
pars9 = mirt(dat9,1,verbose = FALSE) %>% coef_short(.)
pars10= mirt(dat10,1,verbose = FALSE) %>% coef_short(.)

sd(pars6$a)
sd(pars7$a)
sd(pars8$a)
sd(pars9$a)
sd(pars10$a)


# DIF: Calculate ETS classification

ddif = function(x)
  {
  data = x$data
  group = x$group
  n.a1 = sum(data[group=="A",1] == 1)
  n.a0 = sum(data[group=="A",1] == 0)
  n.b1 = sum(data[group=="B",1] == 1)
  n.b0 = sum(data[group=="B",1] == 0)

  alpha = n.a1*n.b0/(n.a0*n.b1)

  re = -2.35 * log(alpha)

  }


df = sim1[sapply(sim1,function(x) x$type=="DIF2PL")] # DIF2PL conditions only

sapply(df, function(x) {

  res = c()

  for (dt in x$datasets) {
  res=c(res,ddif(dt))
  }

  return(mean(abs(res)))
})

# sapply(df, function(x) {
#
# print(x$esize)
#
# })



