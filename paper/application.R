
###########################################################################
# Application Example
###########################################################################

source('C:/Users/admin/Google Drive/4 irt/mmlpwrpackage/paper/plot_functions.R')

library(mmlpwrpackage)
library(dplyr)
library(mirt)
library(tidyr)
library(ggplot2)


# Test the original hypothesis --------------------------------------------

# Load datasets
dat6 <- expand.table(LSAT6)
dat7 <- expand.table(LSAT7)

# Estimate pars
pars6 = mirt(dat6,1,verbose = FALSE) %>% coef_short(.)
pars7 = mirt(dat7,1,verbose = FALSE) %>% coef_short(.)

# Setting up hypotheses
hyp6 <- setup_hypothesis(type = "1PLvs2PL", altpars = pars6)
hyp7 <- setup_hypothesis(type = "1PLvs2PL", altpars = pars7)

# Fitting mirt models
fitted6 = mml.fit(data = dat6, hyp = hyp6)
fitted7 = mml.fit(data = dat7, hyp = hyp7)

# Estimated statistics
stat6 <- c(wald_obs(fitted6),lr_obs(fitted6),score_obs(fitted6))
stat7 <- c(wald_obs(fitted7),lr_obs(fitted7),score_obs(fitted7))

# Test results (p-values)
pvals6 <- pchisq(stat6,df=nrow(hyp6$resmod$Amat),ncp=0,lower.tail=FALSE)
pvals7 <- pchisq(stat7,df=nrow(hyp7$resmod$Amat),ncp=0,lower.tail=FALSE)

#Power given parameters are true
ncps6 <- calculate_ncps(hyp=hyp6)
ncps7 <- calculate_ncps(hyp=hyp7)

power6 = power(hyp=hyp6,ncp=ncps6,alpha=.05,ssize=1000)
power7 = power(hyp=hyp7,ncp=ncps7,alpha=.05,ssize=1000)

round(power6,3)
round(power7,3)

#Sample Size needed to find this effect with a power of 80%
ssize6 = ssize(hyp=hyp6,ncp=ncps6,alpha=.05,power=.80)
ssize7 = ssize(hyp=hyp7,ncp=ncps7,alpha=.05,power=.80)


# Power Curve Plot --------------------------------------------------------

p1 = curve.plot.multi(list(ncps6,ncps7),hyp = hyp6)
pdf("paper/figures/application_power.pdf",height=6,width=8);p1;dev.off()


# Density plot for Introduction -------------------------------------------------------

p2=chisq.dens(ncps7,df=nrow(hyp6$resmod$Amat),n=nrow(dat7),lim=35);p2
pdf("paper/figures/introduction_density.pdf",height=6,width=8);p2;dev.off()


