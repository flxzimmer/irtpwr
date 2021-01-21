###########################################################################
# Prepare Simulation Conditions
###########################################################################

# 4 conditions: (1PLvs2PL,DIF2PL) x (10Items,50Items)
# For each condition: Agreement of Distribution, Agreement of Power

# Calculations of observed stats and expected ncps
# Simulation Conditions are reduced for the calculation of ncps - for conditions differing only in the n.pers argument, only one estimation run is necessary. therefore naming: sim1 and sim1ncp respectively.


## Agreement of Distribution:

# Conditions are given, we will estimate expected ncps and observed stats


# Demo Parameters: reduced simruns, reduced simbased.pers

# prep --------------------------------------------------------------------

library(mmlpwrpackage)
library(dplyr)
library(sn)

set.seed(1234)


# Distribution ------------------------------------------------------------

# Setup a dataframe specifying the condition
sim1 = expand.grid(type=c("1PLvs2PLa1","1PLvs2PLa2"),n.items=c(50),esize=c("no","small","large"),n.pers=c(500,1000,3000),simruns = 500)

# Order the dataframe
sim1 = sim1[order(sim1$type,sim1$n.items,sim1$esize,sim1$n.pers),]

# Convert the conditions to a list
sim1 = sim1  %>% split(., seq(nrow(.))) %>% lapply(.,as.list)

# Create a list of condition specifications for each condition
# Add the settings for the simbased noncentrality parameter estimation
sim1 = lapply(sim1,function(sim1) {
  simbased.runs=1

  if (sim1$n.items==10) {simbased.pers=10^6;simbased.only=FALSE} else {simbased.pers=10^5;simbased.only=TRUE}

  sim1$condition=list(type= sim1$type,n.items = sim1$n.items,esize = sim1$esize,n.pers=sim1$n.pers,simbased.runs = simbased.runs, simbased.pers = simbased.pers,simbased.only = simbased.only)
  return(sim1)
  })

# Add the hypothesis object and the datasets to the list according to the specifications of the conditions.
sim1 = lapply(sim1,function(x) {

  n.items = x$condition$n.items
  esize = x$condition$esize
  type =  x$condition$type
  quantiles = (1:n.items)/(n.items+1)
  perm = sample(n.items)

    #esize for 10 Items
    if (esize=="large") {d1 = .15} else if (esize=="small") {d1 = .1} else if (esize=="no") {d1=0}
    #esize for 50 Items
    if (esize=="large") {d2 = .07} else if (esize=="small") {d2 = .04} else if (esize=="no") {d2=0}

    if(n.items==10){d = d1} else {d = d2}

    x$altpars = list(
      a = qlnorm(quantiles,sdlog = d),
      d = qnorm(quantiles)[perm]
    )

  x$hyp = setup_hypothesis(type = type, altpars = x$altpars)

  x$datasets = list()
  for (i in 1:x$simruns) {
     x$datasets[[i]] = setup.data(hyp = x$hyp, n = x$n.pers)
  }

  return(x)
})

# For the calculation of ncps: Drop datasets for the conditions
sim1ncp = lapply(sim1,function(x) x =x[names(x)!="datasets"])
# The ncps will be identical for conditions with identical alternative hypotheses, dropping unnecessary conditions that vary only in npers
sim1ncp = sim1ncp[sapply(sim1ncp,function(x) x$n.pers==sim1ncp[[1]]$n.pers)]
# Dropping conditions with esize = "no"
sim1ncp = sim1ncp[sapply(sim1ncp,function(x) x$condition$esize!="no")]



# save all ----------------------------------------------------------------

save(sim1, sim1ncp,file= "simsetupa.Rdata")
#save(sim1, sim1ncp, sim2,file= "C:/Users/felix/Google Drive/4 irt/mmlpwr/simsetup.Rdata")


