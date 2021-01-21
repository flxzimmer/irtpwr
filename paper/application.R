
library(mmlpwrpackage)
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


p1 = curve.plot.multi(list(ncps6,ncps7),hyp = hyp6)
pdf("app1.pdf",height=6,width=8);p1;dev.off()


# Density plot for Introduction -------------------------------------------------------

p2=chisq.dens(ncps7,df=nrow(hyp6$resmod$Amat),n=nrow(dat7),lim=35);p2
pdf("density.pdf",height=6,width=8);p2;dev.off()


# Density plot for Introduction -------------------------------------------------------

p3=chisq.dens.null(ncps7,df=nrow(hyp6$resmod$Amat),n=nrow(dat7),lim=35);p3
pdf("density2.pdf",height=6,width=8);p3;dev.off()





# Other -------------------------------------------------------------------
#
# #Density Plot
# p1=chisq.dens(ncp6,df=df,n=nrow(dat6),lim=20)
# grid.arrange(p1,p2,ncol=2)
#
#
#
#
# #Results Table
# tab1 = data.frame(pars=rbind(pars6$a,pars6$d,pars7$a,pars7$d))
# tab1 = cbind(rep(c("LSAT6","LSAT7"),each=2),rep(c("a","d"),2),tab1)
# names(tab1) = c("Dataset","Parameter","Item 1","Item 2","Item 3","Item 4","Item 5")
# xtable(tab1,digits=c(3)) #latex export
#
#
# tab2 = cbind(rbind(stat6,stat7),
#              rbind(pvals6,pvals7),rbind(pw6,pw7),rbind(ssize6,ssize7))
# tab2 = as.data.frame(tab2,row.names = c("LSAT6","LSAT7"))
# names(tab2) = c("w","lr","sc","p.w","p.lr","p.sc","pw.w","pw.lr","pw.sc","n.w","n.lr","n.sc")
#
# t(tab2)
# xtable(t(tab2),digits=c(3)) #latex export




# # Test the original hypothesis --------------------------------------------
#
# dat6 <- expand.table(LSAT6)
# dat7 <- expand.table(LSAT7)
#
# hyp = "2PLvs1PL"
# df= ncol(dat6)-1
#
# # Estimated pars
# pars6 = coefx(mirt(dat6,1,itemtype = c('2PL'),verbose = FALSE),pars="2PL")
# pars7 = coefx(mirt(dat7,1,itemtype = c('2PL'),verbose = FALSE),pars="2PL")
#
# # Estimated statistics
# stat6 = mirt.stats(dat6,hyp) %>% do.call(c,.)
# stat7 = mirt.stats(dat7,hyp) %>% do.call(c,.)
#
# # Test results (p-values)
# pvals6 = stat6 %>% pchisq(.,df=ncol(dat6)-1,ncp=0,lower.tail=FALSE)
# pvals7 = stat7 %>% pchisq(.,df=ncol(dat6)-1,ncp=0,lower.tail=FALSE)
#
# #Power given parameters are true
# ncp6 = c(ncp.wald(pars6,hyp),ncp.lr(pars6,hyp),ncp.sc(pars6,hyp))
# ncp7 = c(ncp.wald(pars7,hyp),ncp.lr(pars7,hyp),ncp.sc(pars7,hyp))
#
# pw6 = power(df,ncp6,alpha=.05,ssize=nrow(dat6))
# pw7 = power(df,ncp7,alpha=.05,ssize=nrow(dat7))
#
# round(pw6*100,1)
# round(pw7*100,1)
#
# #Power Curves
# #p4=curve.plot(ncp7,df=df)
# #grid.arrange(p3,p4,ncol=2)
# p1 = curve.plot.multi(list(ncp6,ncp7),df=df)
# pdf("app1.pdf",height=6,width=8);p1;dev.off()
#
#
# #Sample Size needed to find this effect with a power of 80%
# ssize6 = ssize(df,ncp6 ,alpha=.05,power=.80)
# ssize7 = ssize(df,ncp7 ,alpha=.05,power=.80)
#
#
# # Density for Intro -------------------------------------------------------
#
# #Density Plot
# p2=chisq.dens(ncp7,df=df,n=nrow(dat7),lim=35)
# pdf("density.pdf",height=6,width=8);p2;dev.off()
#
#
# # Other -------------------------------------------------------------------
#
# #Density Plot
# p1=chisq.dens(ncp6,df=df,n=nrow(dat6),lim=20)
# grid.arrange(p1,p2,ncol=2)
#
#
#
#
# #Results Table
# tab1 = data.frame(pars=rbind(pars6$a,pars6$d,pars7$a,pars7$d))
# tab1 = cbind(rep(c("LSAT6","LSAT7"),each=2),rep(c("a","d"),2),tab1)
# names(tab1) = c("Dataset","Parameter","Item 1","Item 2","Item 3","Item 4","Item 5")
# xtable(tab1,digits=c(3)) #latex export
#
#
# tab2 = cbind(rbind(stat6,stat7),
# rbind(pvals6,pvals7),rbind(pw6,pw7),rbind(ssize6,ssize7))
# tab2 = as.data.frame(tab2,row.names = c("LSAT6","LSAT7"))
# names(tab2) = c("w","lr","sc","p.w","p.lr","p.sc","pw.w","pw.lr","pw.sc","n.w","n.lr","n.sc")
#
# t(tab2)
# xtable(t(tab2),digits=c(3)) #latex export
#
#
