library(mmlpwrpackage)
library(mirt)
library(dplyr)


altpars <- list(
  a = rlnorm(50,sdlog = .4),
  d = rnorm(50)
)

# hyp <- setup_hypothesis(type = "1PLvs2PL", altpars = altpars)
# hyp <- setup_hypothesis(type = h_1PLvs2PLa1, altpars = altpars)
hyp <- setup_hypothesis(type = h_1PLvs2PLa2, altpars = altpars)


data <- setup.data(hyp=hyp,n=500)
fitted <- mml.fit(data = data,hyp = hyp)
stats_obs <- c(wald_obs(fitted),lr_obs(fitted),score_obs(fitted))
pvals <- pchisq(stats_obs,df=nrow(hyp$resmod$Amat),ncp=0,lower.tail=FALSE)


stats_obs
pvals



