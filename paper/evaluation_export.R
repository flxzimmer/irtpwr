###########################################################################
# Export simulation results
###########################################################################


source('C:/Users/felix/Google Drive/4 irt/mmlpwrpackage/paper/additional_functions.R')

library(mmlpwrpackage)

load(file= "C:/Users/felix/Google Drive/4 irt/mmlpwr/results1.Rdata")
load(file= "C:/Users/felix/Google Drive/4 irt/mmlpwr/results1ncp.Rdata")


# qqplots -------------------------------------------------------------------

p1 = qq.plot(res.stat = res.stats1,res.ncp = res.ncps1,preset="1PLvs2PL")

pdf("1PLvs2PL_qq_10.pdf");p1;dev.off()






# manage server -----------------------------------------------------------

# file.copy("mmlpwr_eval_setup_temp.R", "Y:/mmlpwr_eval_setup_temp.R", overwrite = T)
# file.copy("mmlpwr_eval.R", "Y:/mmlpwr_eval.R", overwrite = T)
# file.copy("mmlpwr_eval_misspec.R", "Y:/mmlpwr_eval_misspec.R", overwrite = T)
# file.copy("mmlpwr_func5.R", "Y:/mmlpwr_func5.R", overwrite = T)
# file.copy("mmlpwr_func_2PL.R", "Y:/mmlpwr_func_2PL.R", overwrite = T)
# file.copy("mmlpwr_func_3PL.R", "Y:/mmlpwr_func_3PL.R", overwrite = T)
# file.copy("mmlpwr_func_GPCM.R", "Y:/mmlpwr_func_GPCM.R", overwrite = T)

# nohup nice Rscript mmlpwr_eval_setup_temp.R > evalsetup.out &
# nohup nice Rscript mmlpwr_eval.R > eval.out &
# nohup nice Rscript mmlpwr_eval_misspec.R > misspec.out &
# pkill -u fzimmer



# prep --------------------------------------------------------------------

# source('mmlpwr_func5.R')

# load(file= "results - Copy (24).Rdata")
#
# load(file="results1.Rdata")
# load(file="results2.Rdata")

# mis1= res1
# mis2= res2
#
# mis1 = lapply(mis1,function(x) x$stats)
# mis2 = lapply(mis2,function(x) x$stats)
#
# res.stats = lapply(res,function(x) x$stats)
# res.ncps = lapply(res,function(x) x$ncps)
#
# res.stats1 = res.stats[1:18]
# res.stats2 = res.stats[(1:18)+18]
# res.stats3 = res.stats[37:(37+7)]
# res.stats4 = res.stats[45:(45+7)]
#
# res.ncps1 = res.ncps[1:18]
# res.ncps2 = res.ncps[(1:18)+18]
# res.ncps3 = res.ncps[37:(37+7)]
# res.ncps4 = res.ncps[45:(45+7)]


# qqplots -------------------------------------------------------------------

p1 = qq.plot(res.stat = res.stats1,res.ncp = res.ncps1,preset="1PLvs2PL")
p2 = qq.plot(res.stat = res.stats1,res.ncp = res.ncps1,preset="DIF2PL")
p3 = qq.plot(res.stat = res.stats2,res.ncp = res.ncps2,preset="1PLvs2PL")
p4 = qq.plot(res.stat = res.stats2,res.ncp = res.ncps2,preset="DIF2PL")

pdf("1PLvs2PL_qq_10.pdf");p1;dev.off()
pdf("DIF2PL_qq_10.pdf");p2;dev.off()
pdf("1PLvs2PL_qq_50.pdf");p3;dev.off()
pdf("DIF2PL_qq_50.pdf");p4;dev.off()


# powerplots --------------------------------------------------------------

# what were the used sample sizes?
sapply(res.stats3,function(x) paste(x$sim$preset,x$sim$n.pers))
sapply(res.stats4,function(x) paste(x$sim$preset,x$sim$n.pers))

# load(file= "sim4 - Copy (2).Rdata")
# sapply(sim,function(x) x$n.pers)


p5 = power.plot(res.stat = res.stats3,res.ncp = res.ncps3,preset="1PLvs2PL")
p6 = power.plot(res.stat = res.stats3,res.ncp = res.ncps3,preset="DIF2PL")
p7 = power.plot(res.stat = res.stats4,res.ncp = res.ncps4,preset="1PLvs2PL")
p8 = power.plot(res.stat = res.stats4,res.ncp = res.ncps4,preset="DIF2PL")
pdf("1PLvs2PL_pw_10.pdf",height=6,width=8);grid.draw(p5);dev.off()
pdf("DIF2PL_pw_10.pdf",height=6,width=8);grid.draw(p6);dev.off()
pdf("1PLvs2PL_pw_50.pdf",height=6,width=8);grid.draw(p7);dev.off()
pdf("DIF2PL_pw_50.pdf",height=6,width=8);grid.draw(p8);dev.off()


# misspec -----------------------------------------------------------------


p9 = power.plot(res.stat = mis1,res.ncp = res.ncps3,preset="1PLvs2PL",analytical.only = TRUE)
p10 = power.plot(res.stat = mis2,res.ncp = res.ncps3,preset="1PLvs2PL",analytical.only = TRUE)
pdf("missspec.pdf",height=6,width=8);grid.draw(grid.arrange(p9,p10,ncol=2));dev.off()

sapply(mis1,function(x) paste(x$sim$preset,x$sim$n.pers))

  # pdf("miss_unif.pdf",height=6,width=8);grid.draw(p9);dev.off()
# pdf("miss_skewed.pdf",height=6,width=8);grid.draw(p10);dev.off()




# kstable -----------------------------------------------------------------

ks = rbind(ks.table(res.stat = res.stats1,res.ncp = res.ncps1),
  ks.table(res.stat = res.stats2,res.ncp = res.ncps2))

ks= ks[order(ks$Hypothesis,ks$Items),]

ks %>% xtable() %>% print(.,include.rownames=FALSE)


