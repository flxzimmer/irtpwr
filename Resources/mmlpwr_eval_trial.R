
# prep --------------------------------------------------------------------

source('mmlpwr_func5.R')


# Single Runs -------------------------------------------------------------

null.hypothesis = setup.null.hypothesis(n.items=5,preset="1PLvs2PL")

# null.hypothesis = setup.null.hypothesis(n.items=5,preset="PCMvsGPCM")

# null.hypothesis = setup.null.hypothesis(n.items=5,preset="3PLspec")

# null.hypothesis = setup.null.hypothesis(n.items=5,preset="DIF2PL")

alternative.hypothesis = setup.alternative.hypothesis(null.hypothesis = null.hypothesis,esize="small")

dataset = setup.data(null.hypothesis = null.hypothesis,alternative.hypothesis = alternative.hypothesis,n = 1000)
fitted = mirt.fit(null.hypothesis = null.hypothesis,alternative.hypothesis = alternative.hypothesis, dataset = dataset)

wald(fitted,null.hypothesis)
lr(fitted)
score(fitted,null.hypothesis)

ncp.sim(null.hypothesis,alternative.hypothesis,n.pers=100000,runs=10,n=5000)+4

ncp.wald(null.hypothesis,alternative.hypothesis,n=5000,method="mirtOakessim")+4

ncp.wald(null.hypothesis,alternative.hypothesis,n=5000)+4
ncp.lr(null.hypothesis,alternative.hypothesis,n=5000)+4
ncp.score(null.hypothesis,alternative.hypothesis,n=5000)+4


# Simulation Runs ---------------------------------------------------------

null.hypothesis = setup.null.hypothesis(n.items=5,preset="1PLvs2PL")

# null.hypothesis = setup.null.hypothesis(n.items=5,preset="PCMvsGPCM")

# null.hypothesis = setup.null.hypothesis(n.items=5,preset="3PLspec")

# null.hypothesis = setup.null.hypothesis(n.items=5,preset="DIF2PL")

alternative.hypothesis = setup.alternative.hypothesis(null.hypothesis = null.hypothesis,esize="large")


sim.stats(null.hypothesis, alternative.hypothesis,runs=10, n.pers=1000)
  

# choosing ss according to expected power ------------------------------------------

null.hypothesis = setup.null.hypothesis(n.items=15,preset="1PLvs2PL")

# null.hypothesis = setup.null.hypothesis(n.items=5,preset="DIF2PL")

alternative.hypothesis = setup.alternative.hypothesis(null.hypothesis = null.hypothesis,esize="small")

mean.ncp = c(ncp.wald(null.hypothesis,alternative.hypothesis),
ncp.lr(null.hypothesis,alternative.hypothesis),
ncp.score(null.hypothesis,alternative.hypothesis)) %>% mean()

powerset = c(.2,.4,.6,.8)
sample.sizes = sapply(powerset, function(x) ssize(df=null.hypothesis$df,ncp=mean.ncp,alpha=.05,power=x)) 


# Misspecification -------------------------------------------------------

# Different distributions for the person parameter -> what happens?
# Parameters neither follow null or alternative hypothesis -> what happens?


null.hypothesis = setup.null.hypothesis(n.items=5,preset="1PLvs2PL")
# null.hypothesis = setup.null.hypothesis(n.items=5,preset="DIF2PL")
alternative.hypothesis = setup.alternative.hypothesis(null.hypothesis = null.hypothesis,esize="large")



dataset = setup.data(null.hypothesis = null.hypothesis,alternative.hypothesis = alternative.hypothesis,n = 1000,dist.type="unif")
dataset = setup.data(null.hypothesis = null.hypothesis,alternative.hypothesis = alternative.hypothesis,n = 1000,dist.type="skewed")


fitted = mirt.fit(null.hypothesis = null.hypothesis,alternative.hypothesis = alternative.hypothesis, dataset = dataset)

wald(fitted,null.hypothesis)
lr(fitted)
score(fitted,null.hypothesis)

# ncp.wald(null.hypothesis,alternative.hypothesis,n=5000,method="mirtOakessim")+4

ncp.wald(null.hypothesis,alternative.hypothesis,n=1000)+4
ncp.lr(null.hypothesis,alternative.hypothesis,n=1000)+4
ncp.score(null.hypothesis,alternative.hypothesis,n=1000)+4





