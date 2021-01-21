###########################################################################
# Export simulation results
###########################################################################

source('C:/Users/admin/Google Drive/4 irt/mmlpwrpackage/paper/additional_functions.R')
source('C:/Users/felix/Google Drive/4 irt/mmlpwrpackage/paper/additional_functions.R')

library(mmlpwrpackage)
library(car)
library(dplyr)
library(digest)
library(ggplot2)
library(grid)
library(xtable)

load(file= "C:/Users/admin/Google Drive/4 irt/mmlpwr/results1.Rdata")
load(file= "C:/Users/admin/Google Drive/4 irt/mmlpwr/results1ncp.Rdata")
load(file= "C:/Users/admin/Google Drive/4 irt/mmlpwr/results2.Rdata")
load(file= "C:/Users/felix/Google Drive/4 irt/mmlpwr/results2temp.Rdata")
load(file= "C:/Users/felix/Google Drive/4 irt/mmlpwr/results1ncptemp.Rdata")

# merge stats and ncps ----------------------------------------------------

# add the ncp results res1ncp to the other results res
# go over res1, fetch the correct ncps from the ncps results


res.merged = lapply(res2, function(x) {

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


pdf("1PLvs2PL_qq_10_analytical_temp.pdf");p1;dev.off()
pdf("DIF2PL_qq_10_analytical_temp.pdf");p4;dev.off()


# powerplots --------------------------------------------------------------

p7 = power.plot(res = res.merged,analytical = TRUE)
p8 = power.plot(res = res.merged,analytical = FALSE)

pdf("pw_10_temp.pdf",height=6,width=8);grid.draw(p7);dev.off()
pdf("pw_50.pdf",height=6,width=8);grid.draw(p8);dev.off()

