
###########################################################################
# Export simulation results
###########################################################################

source('C:/Users/admin/Google Drive/4 irt/mmlpwrpackage/paper/plot_functions.R')

library(mmlpwrpackage)
library(car)
library(dplyr)
library(digest)
library(ggplot2)
library(grid)
library(xtable)


load(file= "C:/Users/admin/Google Drive/4 irt/mmlpwrpackage/paper/Rdata/results1.Rdata")
load(file= "C:/Users/admin/Google Drive/4 irt/mmlpwrpackage/paper/Rdata/results1ncp.Rdata")
load(file= "C:/Users/admin/Google Drive/4 irt/mmlpwrpackage/paper/Rdata/results2.Rdata")


# merge stats and ncps ----------------------------------------------------

# add the ncp results res1ncp to the other results res
# go over res1, fetch the correct ncps from the ncps results

res.merged = lapply(res1, function(x) {

 index = which(sapply(res1ncp, function(y) digest(y$condition[c("type","n.items","esize")])==digest(x$condition[c("type","n.items","esize")]))) %>% as.numeric()

 if (length(index)>0) {
   x$ncp.analytical = res1ncp[[index]]$analytical * x$condition$n.pers
   x$ncp.simbased = res1ncp[[index]]$simbased * x$condition$n.pers
 } else {
   x$ncp.analytical = c(0,0,0)
   x$ncp.simbased = c(0,0,0)
 }
 x$obs = do.call(rbind,x$obs)
 return(x)
})

res = res.merged


res.merged2 = lapply(res2, function(x) {

   index = which(sapply(res1ncp, function(y) digest(y$condition[c("type","n.items","esize")])==digest(x$condition[c("type","n.items","esize")]))) %>% as.numeric()

   if (length(index)>0) {
      x$ncp.analytical = res1ncp[[index]]$analytical * x$condition$n.pers
      x$ncp.simbased = res1ncp[[index]]$simbased * x$condition$n.pers
   } else {
      x$ncp.analytical = c(0,0,0)
      x$ncp.simbased = c(0,0,0)
   }
   x$obs = do.call(rbind,x$obs)
   return(x)
})



# qqplots -------------------------------------------------------------------

p1 = qq.plot(res = res.merged,type="1PLvs2PL",n.items=10,analytical = TRUE)
p2 = qq.plot(res = res.merged,type="1PLvs2PL",n.items=10,analytical = FALSE)
p3= qq.plot(res = res.merged,type="1PLvs2PL",n.items=50,analytical = FALSE)
p4 = qq.plot(res = res.merged,type="DIF2PL",n.items=10,analytical = TRUE)
p5 = qq.plot(res = res.merged,type="DIF2PL",n.items=10,analytical = FALSE)
p6 = qq.plot(res = res.merged,type="DIF2PL",n.items=50,analytical = FALSE)

p1;p2;p3;p4;p5;p6


pdf("paper/figures/1PLvs2PL_qq_10_analytical.pdf");p1;dev.off()
pdf("paper/figures/1PLvs2PL_qq_10_simbased.pdf");p2;dev.off()
pdf("paper/figures/1PLvs2PL_qq_50_simbased.pdf");p3;dev.off()
pdf("paper/figures/DIF2PL_qq_10_analytical.pdf");p4;dev.off()
pdf("paper/figures/DIF2PL_qq_10_simbased.pdf");p5;dev.off()
pdf("paper/figures/DIF2PL_qq_50_simbased.pdf");p6;dev.off()


# powerplots --------------------------------------------------------------

p7 = power.plot(res = res.merged,analytical = TRUE)
p8 = power.plot(res = res.merged,analytical = FALSE)

pdf("paper/figures/pw_10.pdf",height=6,width=8);grid.draw(p7);dev.off()
pdf("paper/figures/pw_50.pdf",height=6,width=8);grid.draw(p8);dev.off()


# # misspec -----------------------------------------------------------------

# qqplot

p9 = qq.plot.mis(res = res.merged2,type="1PLvs2PL",n.items=10,analytical = TRUE,dist= "unif")
p10 = qq.plot.mis(res = res.merged2,type="1PLvs2PL",n.items=10,analytical = TRUE,dist= "skewed")

pdf("paper/figures/qq_unif.pdf");p9;dev.off()
pdf("paper/figures/qq_skewed.pdf");p10;dev.off()


# power plot
p11 = power.plot.mis(res = res.merged2,analytical = TRUE,dist= "unif")
p12 = power.plot.mis(res = res.merged2,analytical = TRUE,dist= "skewed")

pdf("paper/figures/pw_unif.pdf",height=6,width=8);grid.draw(p11);dev.off()
pdf("paper/figures/pw_skewed.pdf",height=6,width=8);grid.draw(p12);dev.off()



# kstable -----------------------------------------------------------------

# ks = rbind(ks.table(res.stat = res.stats1,res.ncp = res.ncps1),
#   ks.table(res.stat = res.stats2,res.ncp = res.ncps2))

ks1 = ks.table(res = res.merged,type="1PLvs2PL",n.items=10,analytical = TRUE)
ks2 = ks.table(res = res.merged,type="1PLvs2PL",n.items=10,analytical = FALSE)
ks3= ks.table(res = res.merged,type="1PLvs2PL",n.items=50,analytical = FALSE)
ks4 = ks.table(res = res.merged,type="DIF2PL",n.items=10,analytical = TRUE)
ks5 = ks.table(res = res.merged,type="DIF2PL",n.items=10,analytical = FALSE)
ks6 = ks.table(res = res.merged,type="DIF2PL",n.items=50,analytical = FALSE)

ks = rbind(ks1,ks2,ks3,ks4,ks5,ks6)

# ks= ks[order(ks$Hypothesis,ks$Items),]

ks %>% xtable() %>% print(.,include.rownames=FALSE)



# Power Table -------------------------------------------------------------

#Condition (Hypothesis, Number of Items, Sample Size, Method, Statistic), Observed Hit Rate, Expected Hit Rate, Confidence Envelope,

pt1 = power.table(res = res.merged,analytical = TRUE)
pt2 = power.table(res = res.merged,analytical = FALSE)

pt = rbind(pt1,pt2)

pt %>% xtable() %>% print(.,include.rownames=FALSE)



# Power Table Misspec -------------------------------------------------------------

#Condition (Hypothesis, Number of Items, Sample Size, Method, Statistic), Observed Hit Rate, Expected Hit Rate, Confidence Envelope,

pt1 = power.table.mis(res = res.merged2,analytical = TRUE,dist= "unif")
pt2 = power.table.mis(res = res.merged2,analytical = TRUE,dist= "skewed")

pt = rbind(pt1,pt2)

pt %>% xtable() %>% print(.,include.rownames=FALSE)
