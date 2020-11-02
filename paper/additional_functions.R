
load.libraries = function() {

  library(mirt)
  library(dplyr)
  library(spatstat)
  library(Deriv)
  library(digest)
  library(Matrix)
  library(MASS)
  library(sn)
}

# Plots -------------------------------------------------------------------


qq.plot = function(res.stat,res.ncp,preset) {

  setup.data = function(x) {

    df = x$sim$null.hypothesis$df

    n.pers = x$sim$n.pers
    if(!is.null(x$ncp$analytical)) {
      ncp.w = (x$ncp$analytical*n.pers)[1]
      ncp.l = (x$ncp$analytical*n.pers)[2]
      ncp.s = (x$ncp$analytical*n.pers)[3]
    } else {
      # ncp.w = (x$ncp$simbased*n.pers)[1]
      ncp.w = (x$ncp$simbased.wald*n.pers)
      ncp.l = (x$ncp$simbased*n.pers)[2]
      ncp.s = (x$ncp$simbased*n.pers)[3]
    }

    exp.w = qchisq((1:length(x$stats$Wald))/(length(x$stats$Wald)+1) ,df = df, ncp = ncp.w) %>% as.numeric()
    exp.l = qchisq((1:length(x$stats$LR))/(length(x$stats$LR)+1) ,df = df, ncp = ncp.l) %>% as.numeric()
    exp.s = qchisq((1:length(x$stats$Score))/(length(x$stats$Score)+1) ,df = df, ncp = ncp.s) %>% as.numeric()

    obs.w = as.numeric(quantile(x$stats$Wald,(1:length(x$stats$Wald))/(length(x$stats$Wald)+1)))
    obs.l = as.numeric(quantile(x$stats$LR,(1:length(x$stats$LR))/(length(x$stats$LR)+1)))
    obs.s = as.numeric(quantile(x$stats$Score,(1:length(x$stats$Score))/(length(x$stats$Score)+1)))

    ks.w = ks.test(x$stats$Wald, "pchisq",df=df,ncp=ncp.w)$p.value
    ks.l = ks.test(x$stats$LR, "pchisq",df=df,ncp=ncp.l)$p.value
    ks.s = ks.test(x$stats$Score, "pchisq",df=df,ncp=ncp.s)$p.value

    df = data.frame(obs=c(obs.w,obs.l,obs.s),exp=c(exp.w,exp.l,exp.s),kind = factor(c(rep("Wald",length(obs.w)),rep("LR",length(obs.l)),rep("Score",length(obs.s))),levels=c("Wald","LR","Score")),n = n.pers,esize = x$sim$esize)
    df2 = data.frame(ks=c(ks.w,ks.l,ks.s),kind = factor(c("Wald","LR","Score")),n = n.pers,esize = x$sim$esize)
    return(list(df,df2))
  }

  for (i in 1:length(res.stat)) {
    res.stat[[i]]$ncp = res.ncp[[i]]
  }

  indi = sapply(res.stat,function(x) isTRUE(x$sim$preset==preset))
  resx = res.stat[indi]


  datax = lapply(resx,FUN=setup.data)
  data = do.call(rbind,lapply(datax,function(x) x[[1]]))
  dataks = do.call(rbind,lapply(datax,function(x) x[[2]]))
  dataks$ks=dataks$ks*nrow(dataks)
  dataks$sig =""

  dataks$sig[dataks$ks<.05] = paste(dataks$kind[dataks$ks<.05],"*",sep="")


  sequence = seq(1,nrow(dataks),3)
  dataks$sig[sequence] = do.call(rbind,lapply(sequence,function(x) paste(dataks$sig[x:(x+2)],collapse=" ")))
  dataks$sig[c(sequence+1,sequence+2)] = ""

  # New facet label names for dose variable
  n.labs <- paste("n =",unique(data$n))
  names(n.labs) <- as.character(unique(data$n))

  # New facet label names for supp variable
  esize.labs <- c("No effect", "Small effect","Large effect")
  names(esize.labs) <- c("no", "small","large")

  re = ggplot(data, aes(x=exp, y=obs,group=kind)) +
    geom_point(aes(color=kind),alpha = 1) +
    xlab("Expected Quantiles")+  ylab("Observed Quantiles") + geom_abline(linetype="dashed", color = "black") + theme(legend.position = "bottom") + facet_grid(vars(esize),vars(n),scale="free", labeller = labeller(n = n.labs, esize = esize.labs)) + labs(color = "Statistic") + scale_color_grey(start=.0,end=.6) + theme_bw() +
    geom_text(dataks,mapping=aes(x=Inf,y=-Inf,label=sig),hjust=1.2,vjust=-1 )


  return(re)
}



ks.table = function(res.stat,res.ncp,preset) {
  #COLS: Hypothesis, Effect Size, Sample Size, Wald(D,p), LR(D,p), Score(D,p)
  # *significant after bonferroni correction

  setup.data = function(x) {

    df = x$sim$null.hypothesis$df

    n.pers = x$sim$n.pers

    if(!is.null(x$ncp$analytical)) {
      ncp.w = (x$ncp$analytical*n.pers)[1]
      ncp.l = (x$ncp$analytical*n.pers)[2]
      ncp.s = (x$ncp$analytical*n.pers)[3]
    } else {
      # ncp.w = (x$ncp$simbased*n.pers)[1]
      ncp.w = (x$ncp$simbased.wald*n.pers)
      ncp.l = (x$ncp$simbased*n.pers)[2]
      ncp.s = (x$ncp$simbased*n.pers)[3]

    }

    ks.wx = ks.test(x$stats$Wald, "pchisq",df=df,ncp=ncp.w)
    ks.lx = ks.test(x$stats$LR, "pchisq",df=df,ncp=ncp.l)
    ks.sx = ks.test(x$stats$Score, "pchisq",df=df,ncp=ncp.s)


    ds = format(round(c(ks.wx$statistic,ks.lx$statistic,ks.sx$statistic),3),nsmall=3)
    ps = paste(format(round(c(ks.wx$p.value,ks.lx$p.value,ks.sx$p.value),3),nsmall=3),c(if(ks.wx$p.value*27<.05) {"*"}else{" "},if(ks.lx$p.value*27<.05) {"*"}else{" "},if(ks.sx$p.value*27<.05) {"*"}else{" "}),sep="")

    df = c(
      hyp = as.character(x$sim$null.hypothesis$preset),
      items = x$sim$n.items,
      esize = as.character(x$sim$esize),
      n.pers = x$sim$n.pers,
      d1 = ds[1],p1 = ps[1],d1 = ds[2],p1 = ps[2],d1 = ds[3],p1 = ps[3]
    )

    names(df) = c("Hypothesis", "Items","Effect size", "Sample size","D","p","D","p","D","p")
    return(df)
  }

  for (i in 1:length(res.stat)) {
    res.stat[[i]]$ncp = res.ncp[[i]]
  }

  # indi = sapply(res.stat,function(x) isTRUE(x$sim$preset==preset))
  # resx = res.stat[indi]
  resx = res.stat

  re = as.data.frame(do.call(rbind,lapply(resx,FUN=setup.data)))

  return(re)
}



power.plot = function(res.stat,res.ncp,preset,alpha=.05,analytical.only=FALSE) {

  pw.sd = function(ncp,df,n.runs,ci=.95,alpha=.05) {

    crit = qchisq(1-alpha,df = df, ncp = 0) %>% as.numeric()

    res0 = c()

    for (i in 1:1000) {
      res0[i]= mean(rchisq(n.runs,df=df,ncp=ncp)>crit)
    }

    re = quantile(res0,c((1-ci)/2,1-(1-ci)/2))

    return(re)
  }

  pw.sd2 = function(pw,n.runs,ci=.95) {

    res = sapply(1:10000,function(x) {(runif(n.runs) < pw) %>% mean() })

    re = quantile(res,c((1-ci)/2,1-(1-ci)/2))

    return(re)
  }

  setup.data = function(x,simbased.only) {
    df = x$sim$null.hypothesis$df

    runs = length(x$stats$Wald)
    n.pers = x$sim$n.pers

    crit = qchisq(1-alpha,df = df, ncp = 0) %>% as.numeric()
    pwro.w = mean(x$stats$Wald>crit)
    pwro.l = mean(x$stats$LR>crit)
    pwro.s = mean(x$stats$Score>crit)

    # ncp.w1 = (x$ncp$simbased*n.pers)[1]
    ncp.w1 = (x$ncp$simbased.wald*n.pers)
    ncp.l1 = (x$ncp$simbased*n.pers)[2]
    ncp.s1 = (x$ncp$simbased*n.pers)[3]

    pwre.w1 = 1 - pchisq(q = crit, df = df, ncp = ncp.w1)
    pwre.l1 = 1 - pchisq(q = crit, df = df, ncp = ncp.l1)
    pwre.s1 = 1 - pchisq(q = crit, df = df, ncp = ncp.s1)

    # pw.w.sd1 = pw.sd(ncp=ncp.w1,df=df,n.runs=runs)
    # pw.l.sd1 = pw.sd(ncp=ncp.l1,df=df,n.runs=runs)
    # pw.s.sd1 = pw.sd(ncp=ncp.s1,df=df,n.runs=runs)
    pw.w.sd1 = pw.sd2(pwro.w,n.runs=runs)
    pw.l.sd1 = pw.sd2(pwro.s,n.runs=runs)
    pw.s.sd1 = pw.sd2(pwro.l,n.runs=runs)

    if (!simbased.only) {
      ncp.w = (x$ncp$analytical*n.pers)[1]
      ncp.l = (x$ncp$analytical*n.pers)[2]
      ncp.s = (x$ncp$analytical*n.pers)[3]

      pwre.w = 1 - pchisq(q = crit, df = df, ncp = ncp.w)
      pwre.l = 1 - pchisq(q = crit, df = df, ncp = ncp.l)
      pwre.s = 1 - pchisq(q = crit, df = df, ncp = ncp.s)

      # pw.w.sd = pw.sd(ncp=ncp.w,df=df,n.runs=runs)
      # pw.l.sd = pw.sd(ncp=ncp.l,df=df,n.runs=runs)
      # pw.s.sd = pw.sd(ncp=ncp.s,df=df,n.runs=runs)

      pw.w.sd = pw.w.sd1
      pw.l.sd = pw.l.sd1
      pw.s.sd = pw.s.sd1

      re = data.frame(
        stat = factor(c("Wald","LR","Score"),levels = c("Wald","LR","Score")),
        pw.exp= c(pwre.w ,pwre.l,pwre.s),
        pw.exp1= c(pwre.w1 ,pwre.l1,pwre.s1),
        pw.obs= c(pwro.w ,pwro.l,pwro.s),
        # pw.lower= c(pw.sdx[1],pw.sdx[1],pw.sdx[1]),
        # pw.upper= c(pw.sdx[2],pw.sdx[2],pw.sdx[2]),
        pw.lower= c(pw.w.sd[1],pw.l.sd[1],pw.s.sd[1]),
        pw.upper= c(pw.w.sd[2],pw.l.sd[2],pw.s.sd[2]),
        # pw.lower1= c(pw.sdx1[1],pw.sdx1[1],pw.sdx1[1]),
        # pw.upper1= c(pw.sdx1[2],pw.sdx1[2],pw.sdx1[2]))
        pw.lower1= c(pw.w.sd1[1],pw.l.sd1[1],pw.s.sd1[1]),
        pw.upper1= c(pw.w.sd1[2],pw.l.sd1[2],pw.s.sd1[2]))
    }


    if(simbased.only) {
      re = data.frame(
        stat = factor(c("Wald","LR","Score"),levels = c("Wald","LR","Score")),
        pw.exp1= c(pwre.w1 ,pwre.l1,pwre.s1),
        pw.obs= c(pwro.w ,pwro.l,pwro.s),
        pw.lower1= c(pw.w.sd1[1],pw.l.sd1[1],pw.s.sd1[1]),
        pw.upper1= c(pw.w.sd1[2],pw.l.sd1[2],pw.s.sd1[2]))
    }
    return(re)
  }

  for (i in 1:length(res.stat)) {
    res.stat[[i]]$ncp = res.ncp[[i]]
  }
  simbased.only = is.null(res.ncp[[1]]$analytical)

  indi = sapply(res.stat,function(x) isTRUE(x$sim$preset==preset))
  resx = res.stat[indi]

  dat = lapply(resx,function(x) setup.data(x,simbased.only)) %>% do.call(rbind,.) %>% as.data.frame()

  re = ggplot(dat, aes(pw.exp1, pw.obs, colour = stat)) +
    geom_point(size=1) +
    geom_abline(linetype="dashed", color = "black") +
    geom_errorbar(data=dat, mapping= aes(ymin = pw.lower1, ymax = pw.upper1),width=.05) +
    scale_color_grey(start=.0,end=.6) + labs(colour = "Statistic") + theme_bw() + scale_y_continuous("Observed Power",breaks=seq(0,1,.2),limits = c(0, 1)) + scale_x_continuous("Expected Power",breaks=seq(0,1,.2),limits = c(0, 1))

  if(!simbased.only) {

    re2 = ggplot(dat, aes(pw.exp, pw.obs, colour = stat)) +
      geom_point(size=1) +
      geom_abline(linetype="dashed", color = "black") +
      geom_errorbar(data=dat, mapping= aes(ymin = pw.lower, ymax = pw.upper),width=.05) +
      scale_color_grey(start=.0,end=.6)  + xlab("Expected Power") + labs(colour = "Statistic") + theme_bw() + scale_y_continuous("Observed Power",breaks=seq(0,1,.2),limits = c(0, 1)) + scale_x_continuous("Expected Power",breaks=seq(0,1,.2),limits = c(0, 1))

    re = grid.arrange(re2,re,ncol=2)

    if(isTRUE(analytical.only)) {re = re2}
  }
  return(re)
}


curve.plot.multi = function(ncpslist,df,alpha=.05,nolegend=FALSE,scales="free") {
  #ncps is a list here


  ncps = ncpslist[[1]]

  upper =  uniroot(f = function(x) .999 - power(df=df,ncp=ncps[1],alpha=alpha,ssize=x),interval = c(0, 1000000),tol = .Machine$double.eps ^ 0.5)$root

  ssizes = seq(0,upper,10)
  res = list()
  for (i in 1:length(ssizes)) {
    res[[i]] = sapply(ncps,function (x) power(df=df,ncp=x,alpha=alpha,ssize=ssizes[i]))
  }

  dat = data.frame(do.call(rbind,res))

  names(dat) = c("Wald","LR","Score")
  dat$ssize = ssizes
  dat1 = as.data.frame(pivot_longer(dat,cols=names(dat)[1:3]))


  ncps = ncpslist[[2]]

  upper =  uniroot(f = function(x) .999 - power(df=df,ncp=ncps[1],alpha=alpha,ssize=x),interval = c(0, 1000000),tol = .Machine$double.eps ^ 0.5)$root

  ssizes = seq(0,upper,10)
  res = list()
  for (i in 1:length(ssizes)) {
    res[[i]] = sapply(ncps,function (x) power(df=df,ncp=x,alpha=alpha,ssize=ssizes[i]))
  }

  dat = data.frame(do.call(rbind,res))
  names(dat) = c("Wald","LR","Score")
  dat$ssize = ssizes
  dat2 = as.data.frame(pivot_longer(dat,cols=names(dat)[1:3]))

  dat1$set = "LSAT6"
  dat2$set = "LSAT7"
  dat = rbind(dat1,dat2)
  dat$name = factor(dat$name,levels=c("Wald","LR","Score"))

  re = ggplot(dat, aes(x=ssize, y=value, group=name)) +
    geom_line(aes(color=name),size=1) + scale_color_grey(start=.0,end=.6) + labs(title = "")+  xlab("Sample size")+ ylab("Power") +labs(color = "Statistic")  + facet_wrap(~set,scales = scales) + scale_y_continuous("Power",breaks=seq(0,1,.2)) + theme_bw()


  return(re)
}



chisq.dens = function(ncp,df,n,alpha=.05,lim=30) {

  crit= qchisq(1-alpha,df)

  p = ggplot(data.frame(x = c(0, lim)), aes(x)) +

    stat_function(
      fun = dchisq,
      geom = "area",
      # alpha = .5,
      args = list(
        df=df
      ),n = 5000,aes(fill = "Null")
    ) +
    stat_function(
      fun = dchisq,
      geom = "area",
      # alpha = .5,
      xlim = c(crit,lim),
      args = list(
        df=df,
        ncp=ncp[1]*n
      ),aes(fill = "Wald")
    ) +
    stat_function(
      fun = dchisq,
      geom = "area",
      # alpha = .1,
      xlim = c(crit,lim),
      args = list(
        df=df,
        ncp=ncp[2]*n
      ),aes(fill = "LR")
    ) +
    stat_function(
      fun = dchisq,
      geom = "area",
      # alpha = .1,
      xlim = c(crit,lim),
      args = list(
        df=df,
        ncp=ncp[3]*n
      ),aes(fill = "Score")
    ) +
    stat_function(
      fun = dchisq,
      geom = "line",
      args = list(
        df=df
      ),n = 5000
    )+
    stat_function(
      fun = dchisq,
      geom = "line",
      args = list(
        df=df,
        ncp=ncp[1]*n
      )
    ) +
    stat_function(
      fun = dchisq,
      geom = "line",
      args = list(
        df=df,
        ncp=ncp[2]*n
      )
    ) +
    stat_function(
      fun = dchisq,
      geom = "line",
      args = list(
        df=df,
        ncp=ncp[3]*n
      )
    ) +
    labs(
      title = "",
      x = "Statistic",
      y = "Density"
    ) + scale_fill_manual("",breaks=c("Null","Wald","LR","Score"),values = c("grey90","black","grey30","grey60")) + theme_bw()

  return(p)
}


