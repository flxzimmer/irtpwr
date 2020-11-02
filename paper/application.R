
library(mmlpwrpackage)
library(mirt)

# Test the original hypothesis --------------------------------------------

dat6 <- expand.table(LSAT6)
dat7 <- expand.table(LSAT7)

# Estimated pars
pars6 = mirt(dat6,1,itemtype = c('2PL'),verbose = FALSE) %>% coef.short(.,itemtype="2PL")
pars7 = mirt(dat7,1,itemtype = c('2PL'),verbose = FALSE) %>% coef.short(.,itemtype="2PL")

# Setting up hypotheses
null.hyp = setup.null.hypothesis(n.items=5,preset="1PLvs2PL")
alt.hyp6 = setup.alternative.hypothesis(null.hypothesis = null.hyp,altpars=pars6)
alt.hyp7 = setup.alternative.hypothesis(null.hypothesis = null.hyp,altpars=pars7)

# Fitting mirt models
fitted6 = mirt.fit(null.hypothesis = null.hyp,alternative.hypothesis = alt.hyp6, dataset = dat6)
fitted7 = mirt.fit(null.hypothesis = null.hyp,alternative.hypothesis = alt.hyp7, dataset = dat7)

# Estimated statistics
stat6 = c(wald(fitted6,null.hyp),lr(fitted6),score(fitted6,null.hyp))
stat7 = c(wald(fitted7,null.hyp),lr(fitted7),score(fitted7,null.hyp))

# Test results (p-values)
pvals6 = stat6 %>% pchisq(.,df=null.hyp$df,ncp=0,lower.tail=FALSE)
pvals7 = stat7 %>% pchisq(.,df=null.hyp$df,ncp=0,lower.tail=FALSE)

#Power given parameters are true
ncp6 = c(ncp.wald(null.hyp,alt.hyp6),ncp.lr(null.hyp,alt.hyp6),ncp.score(null.hyp,alt.hyp6))
ncp7 = c(ncp.wald(null.hyp,alt.hyp7),ncp.lr(null.hyp,alt.hyp7),ncp.score(null.hyp,alt.hyp7))

pw6 = power(df=null.hyp$df,ncp6,alpha=.05,ssize=nrow(dat6))
pw7 = power(df=null.hyp$df,ncp7,alpha=.05,ssize=nrow(dat7))

round(pw6,3)
round(pw7,3)

#Sample Size needed to find this effect with a power of 80%
ssize6 = ssize(null.hyp$df,ncp6 ,alpha=.05,power=.80)
ssize7 = ssize(null.hyp$df,ncp7 ,alpha=.05,power=.80)


#Power Curves
#p4=curve.plot(ncp7,df=df)
#grid.arrange(p3,p4,ncol=2)
p1 = curve.plot.multi(list(ncp6,ncp7),df=null.hyp$df)
pdf("app1.pdf",height=6,width=8);p1;dev.off()



# Density for Intro -------------------------------------------------------

#Density Plot
p2=chisq.dens(ncp7,df=null.hyp$df,n=nrow(dat7),lim=35)
pdf("density.pdf",height=6,width=8);p2;dev.off()


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
