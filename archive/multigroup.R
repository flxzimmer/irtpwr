
# Testing the application of mirt multigroup command
# Could it be specified differently? Did we make an implicit decision that needs to be reported?

library(mmlpwrpackage)
library(mirt)
library(dplyr)


group1 = group2 <- list(
  a = rlnorm(5,sdlog = .2),
  d = rnorm(5)
)

group2$a[1] = (group2$a[1])^2
group2$d[1] = group2$d[1] + .5

altpars <- list(group1,group2)

hyp <- setup_hypothesis(type = "DIF2PL", altpars = altpars)

data <- setup.data(hyp=hyp,n=500)


group = data$group
data = data$data
n.items = ncol(data)

model = mirt::mirt.model(paste('F = 1-',n.items,'
                       CONSTRAINB = (2-',n.items,', d), (2-',n.items,', a1)'))

# as in paper
res = mirt::multipleGroup(data,model = model,itemtype = "2PL",technical = list(NCYCLES = 5000),verbose = FALSE,group=group)
logLik(res)
coef(res)

# free means
res1 = mirt::multipleGroup(data,model = model,itemtype = "2PL",technical = list(NCYCLES = 5000),verbose = FALSE,group=group,invariance = "free_means")
logLik(res1)
coef(res1)

# + free vars
res1 = mirt::multipleGroup(data,model = model,itemtype = "2PL",technical = list(NCYCLES = 5000),verbose = FALSE,group=group,invariance = c("free_means","free_var"))
logLik(res1)
coef(res1)



# simulation --------------------------------------------------------------

#simulation: wird der itemparameter-unterschied vom Personenparameter Mittelwert "geschluckt"?
# resultat: Unterschied ist größer wenn gruppen-personenparameter frei geschätzt werden.

re=list()

for (i in 1:100) {

  data <- setup.data(hyp=hyp,n=500)
  group = data$group
  data = data$data

  # as in paper
  res = mirt::multipleGroup(data,model = model,itemtype = "2PL",technical = list(NCYCLES = 5000),verbose = FALSE,group=group)

  # free means
  res1 = mirt::multipleGroup(data,model = model,itemtype = "2PL",technical = list(NCYCLES = 5000),verbose = FALSE,group=group,invariance = "free_means")

  a = abs(coef(res)$A$Item_1[2]-coef(res)$B$Item_1[2])
  b = abs(coef(res1)$A$Item_1[2]-coef(res1)$B$Item_1[2])

  re[[i]] = c(a,b)

  }

a = do.call(rbind,re)

sum(a[,1]-a[,2])



# Using Package Code and Tests --------------------------------------------

library(mmlpwrpackage)
library(mirt)
library(dplyr)

group1 = group2 <- list(
  a = rlnorm(5,sdlog = .2),
  d = rnorm(5)
)

group2$a[1] = (group2$a[1])^2
group2$d[1] = group2$d[1] + .5

altpars <- list(group1,group2)

hyp <- setup_hypothesis(type = "DIF2PL", altpars = altpars)

data <- setup.data(hyp=hyp,n=500)

fitted <- mml.fit(data = data,hyp = hyp,free_mean = FALSE)
wald_obs(fitted)

fitted <- mml.fit(data = data,hyp = hyp,free_mean = TRUE)
wald_obs(fitted)





