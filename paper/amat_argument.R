
# A mat comp --------------------------------------------------------------


altpars <- list(
  a = rlnorm(5,sdlog = .4),
  d = rnorm(5)
)

hyp <- setup_hypothesis(type = "1PLvs2PL", altpars = altpars)
hypa <- setup_hypothesis(type = h_1PLvs2PLa, altpars = altpars)

data <- setup.data(hyp=hyp,n=500)

## Performing Hypothesis Tests

fitted <- mml.fit(data = data,hyp = hyp)
fitteda <- mml.fit(data = data,hyp = hypa)

c(wald_obs(fitted),lr_obs(fitted),score_obs(fitted))
c(wald_obs(fitteda),lr_obs(fitteda),score_obs(fitteda))


ncps <- calculate_ncps(hyp=hyp)
power(hyp=hyp,ncp=ncps,alpha=.05,ssize=500)
ssize(hyp=hyp,ncp=ncps,alpha=.05,power=.80)


ncpsa <- calculate_ncps(hyp=hypa)
power(hyp=hypa,ncp=ncpsa,alpha=.05,ssize=500)
ssize(hyp=hypa,ncp=ncpsa,alpha=.05,power=.80)



