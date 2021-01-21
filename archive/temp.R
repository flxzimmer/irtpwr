## Usage without data

n.items = 50
quantiles = (1:n.items)/(n.items+1)
perm = sample(n.items)

# if (esize=="large") {d1 = .15} else if (esize=="small") {d1 = .1} else if (esize=="no") {d1=0}
# if (esize=="large") {d2 = .1} else if (esize=="small") {d2 = .03} else if (esize=="no") {d2=0}

altpars = list(
  a = qlnorm(quantiles,sdlog = .07),
  d = qnorm(quantiles)[perm]
)


hyp <- setup_hypothesis(type = "1PLvs2PL", altpars = altpars)

ncps <- mmlpwrpackage:::ncp.sim.wald(hyp,n=1,10000)
power(hyp=hyp,ncp=ncps,alpha=.05,ssize=c(500,1000,3000))


## Usage without data

n.items = 50
quantiles = (1:n.items)/(n.items+1)
perm = sample(n.items)

# if (esize=="large") {d1 = .15} else if (esize=="small") {d1 = .1} else if (esize=="no") {d1=0}
# if (esize=="large") {d2 = .1} else if (esize=="small") {d2 = .03} else if (esize=="no") {d2=0}

altpars = list(
  a = qlnorm(quantiles,sdlog = .04),
  d = qnorm(quantiles)[perm]
)


hyp <- setup_hypothesis(type = "1PLvs2PL", altpars = altpars)

ncps <- mmlpwrpackage:::ncp.sim.wald(hyp,n=1,10000)
power(hyp=hyp,ncp=ncps,alpha=.05,ssize=c(500,1000,3000))


