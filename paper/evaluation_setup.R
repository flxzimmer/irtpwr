###########################################################################
# Prepare Simulation Conditions
###########################################################################

# 4 conditions: (1PLvs2PL,DIF2PL) x (10Items,50Items)
# For each condition: Agreement of Distribution, Agreement of Power

# Calculations of observed stats and expected ncps
# Simulation Conditions are reduced for the calculation of ncps - for conditions differing only in the n.pers argument, only one estimation run is necessary. therefore naming: sim1 and sim1ncp respectively.


## Agreement of Distribution:

# Conditions are given, we will estimate expected ncps and observed stats

## Agreement of Power:

# Conditions are not fully determined, sample size will be determined according to ncps. One simulation "prerun" will be to determine the sample sizes, which are then added to the conditions for a second run analogous to the agreement of distribution.


# prep --------------------------------------------------------------------

library(mmlpwrpackage)
library(dplyr)
library(sn)

set.seed(123)


# Distribution ------------------------------------------------------------

# Make a list of all conditions.
sim1 = expand.grid(type=c("1PLvs2PL","DIF2PL"),n.items=c(10,50),esize=c("no","small","large"),n.pers=c(500,1000,3000),simruns = 500)
sim1 = sim1[order(sim1$type,sim1$n.items,sim1$esize,sim1$n.pers),]
sim1$condition = list(type = sim1$type,n.items = sim1$n.items,esize = sim1$esize)
sim1 = sim1  %>% split(., seq(nrow(.))) %>% lapply(.,as.list)

# Add the hypothesis object and the datasets to the list according to the specifications of the conditions.
sim1 = lapply(sim1,function(x) {

  n.items = x$n.items
  esize = x$esize
  quantiles = (1:n.items)/(n.items+1)
  perm = sample(n.items)

  if (x$type == "1PLvs2PL") {

    if (esize=="large") {d = .1} else if (esize=="small") {d = .05} else if (esize=="xsmall") {d=0.03} else if (esize=="no") {d=0}

    x$altpars = list(
      a = qlnorm(quantiles,sdlog = d),
      d = qnorm(quantiles)[perm]
    )
  }

  if (x$type == "DIF2PL") {

    if (esize=="large") {d = .4} else if (esize=="small") {d = .2} else if (esize=="no") {d=0}

    group1 = group2 <- list(
      a = qlnorm(quantiles,sdlog = .1),
      d = qnorm(quantiles)[perm]
    )
    group1$a[1] = 1+d/2
    group1$d[1] = d/2
    group2$a[1] = 1-d/2
    group2$d[1] = -d/2

    x$altpars <- list(group1,group2)

  }

  x$hyp = setup_hypothesis(type = x$type, altpars = x$altpars)

  x$datasets = list()
  for (i in 1:x$simruns) {
    x$datasets[[i]] = setup.data(hyp = x$hyp, n = x$n.pers)
  }

  return(x)
})


# The ncps will be identical for conditions with identical alternative hypotheses, dropping unnecessary npers conditions
sim1ncp = sim1[sapply(sim1,function(x) x$n.pers==sim1[[1]]$n.pers)]


# Power ------------------------------------------------------------

sim2 = expand.grid(type=c("1PLvs2PL","DIF2PL"),n.items=c(10,50),esize ="",simruns = 500)
sim2$esize[sim2$type =="1PLvs2PL"& sim2$n.items == 10] = "large"
sim2$esize[sim2$type =="DIF2PL"& sim2$n.items == 10] = "large"
sim2$esize[sim2$type =="1PLvs2PL"& sim2$n.items == 50] = "xsmall"
sim2$esize[sim2$type =="DIF2PL"& sim2$n.items == 50] = "small"

sim2 = sim2[order(sim2$type,sim2$n.items,sim2$esize,sim2$n.pers),]
sim2$condition = list(type= sim2$type,n.items = sim2$n.items,esize = sim2$esize)
sim2 = sim2  %>% split(., seq(nrow(.))) %>% lapply(.,as.list)


# Add the hypothesis object to the list according to the specifications of the conditions.
sim2ncp = lapply(sim2,function(x) {

  n.items = x$n.items
  esize = x$esize
  quantiles = (1:n.items)/(n.items+1)
  perm = sample(n.items)

  if (x$type == "1PLvs2PL") {

    if (esize=="large") {d = .1} else if (esize=="small") {d = .05} else if (esize=="xsmall") {d=0.03} else if (esize=="no") {d=0}

    x$altpars = list(
      a = qlnorm(quantiles,sdlog = d),
      d = qnorm(quantiles)[perm]
    )
  }

  if (x$type == "DIF2PL") {

    if (esize=="large") {d = .4} else if (esize=="small") {d = .2} else if (esize=="no") {d=0}

    group1 = group2 <- list(
      a = qlnorm(quantiles,sdlog = .1),
      d = qnorm(quantiles)[perm]
    )
    group1$a[1] = 1+d/2
    group1$d[1] = d/2
    group2$a[1] = 1-d/2
    group2$d[1] = -d/2

    x$altpars <- list(group1,group2)

  }

  x$hyp = setup_hypothesis(type = x$type, altpars = x$altpars)

  return(x)
})


# Misspecification --------------------------------------------------------

# setting up functions to generate the non-normal distributions.

unifdist = function(x) runif(x,min=-sqrt(3),max=sqrt(3))

skeweddist = function(x) {
  #skeweddist: alpha = 0 entspricht normalverteilung.
  # Wähle alpha = 4 und passe alle anderen Werte an, damit mean~0 und sd~1 gegeben ist.
  alpha = 4
  delta = alpha/sqrt(1+alpha^2)
  omega = sqrt(1/(1-2*delta^2/pi))
  xi = -omega*delta*sqrt(2/pi)
  return(rsn(n=x, xi=xi, omega=omega, alpha=alpha))
  }

# check the properties of the distributions (should be mean~0 and sd~1)
a= skeweddist(10^6);mean(a);sd(a);hist(a)
a= unifdist(10^6);mean(a);sd(a)


# save all ----------------------------------------------------------------

save(sim1, sim1ncp, sim2ncp, unifdist, skeweddist,file= "C:/Users/felix/Google Drive/4 irt/mmlpwr/simsetup.Rdata")


