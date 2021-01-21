
# Plots -------------------------------------------------------------------

qq.plot = function(res,type="DIF2PL",n.items=10,analytical = TRUE) {

  setup.data = function(x) {

    df = x$hyp$resmod$Amat %>% nrow()
    n.obs = nrow(x$obs)

    obs.w = quantile(x$obs[,1],(1:n.obs)/(n.obs+1)) %>% as.numeric()
    obs.l = quantile(x$obs[,2],(1:n.obs)/(n.obs+1)) %>% as.numeric()
    obs.s = quantile(x$obs[,3],(1:n.obs)/(n.obs+1)) %>% as.numeric()

    if(isTRUE(analytical)){
      ncp.w = x$ncp.analytical[1]
      ncp.l = x$ncp.analytical[2]
      ncp.s = x$ncp.analytical[3]
    } else {
      ncp.w = x$ncp.simbased[1]
      ncp.l = x$ncp.simbased[2]
      ncp.s = x$ncp.simbased[3]
    }
    exp.w = qchisq((1:n.obs)/(n.obs+1),df = df, ncp = ncp.w) %>% as.numeric()
    exp.l = qchisq((1:n.obs)/(n.obs+1),df = df, ncp = ncp.l) %>% as.numeric()
    exp.s = qchisq((1:n.obs)/(n.obs+1),df = df, ncp = ncp.s) %>% as.numeric()

    ks.w = ks.test(x$obs[,1], "pchisq",df=df,ncp=ncp.w)$p.value
    ks.l = ks.test(x$obs[,2], "pchisq",df=df,ncp=ncp.l)$p.value
    ks.s = ks.test(x$obs[,3], "pchisq",df=df,ncp=ncp.s)$p.value

    if(is.na(ks.w)){ks.w =0.10}
    if(is.na(ks.l)){ks.l =0.10}
    if(is.na(ks.s)){ks.s =0.10}

    df = data.frame(obs=c(obs.w,obs.l,obs.s),exp=c(exp.w,exp.l,exp.s),kind = factor(c(rep("Wald",length(obs.w)),rep("LR",length(obs.l)),rep("Score",length(obs.s))),levels=c("Wald","LR","Score")),n = x$n.pers,esize = x$esize)
    df2 = data.frame(ks=c(ks.w,ks.l,ks.s),kind = factor(c("Wald","LR","Score")),n = x$n.pers,esize = x$esize)
    return(list(df,df2))
  }

  # Get the indices of the desired condition
  indices = which(sapply(res,function(x) x$type == type & x$n.items == n.items ))
  resx = res[indices]

  #Setup the data and ks tests
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


  esize.labs <- c("H0: No effect","H1: Small effect","H1: Large effect")
  names(esize.labs) <- c("no", "small","large")

  # Generate Title
  if(type=="1PLvs2PL") {t.type = "Rasch vs 2PL"} else {"t.type" = "DIF"}
  if(isTRUE(analytical)){t.sim="analytical power"}else{t.sim="simulation-based power"}
  title = paste(t.type,", ",n.items," items, ",t.sim,sep="")

  re = ggplot(data, aes(x=exp, y=obs,group=kind)) +
    geom_point(aes(color=kind,shape=kind),alpha = 1) +
    xlab("Expected Quantiles")+  ylab("Observed Quantiles") + geom_abline(linetype="dashed", color = "black") + theme(legend.position = "bottom") + facet_grid(vars(esize),vars(n),scale="free", labeller = labeller(n = n.labs, esize = esize.labs)) +
    labs(colour = "Statistic",shape = "Statistic") +
    theme_bw() + ggtitle(title)+
    geom_text(dataks,mapping=aes(x=Inf,y=-Inf,label=sig),hjust=1.2,vjust=-1 )


  return(re)
}

qq.plot.mis = function(res,type="DIF2PL",n.items=10,analytical = TRUE,dist = "unif") {

  setup.data = function(x) {

    df = x$hyp$resmod$Amat %>% nrow()
    n.obs = nrow(x$obs)

    obs.w = quantile(x$obs[,1],(1:n.obs)/(n.obs+1)) %>% as.numeric()
    obs.l = quantile(x$obs[,2],(1:n.obs)/(n.obs+1)) %>% as.numeric()
    obs.s = quantile(x$obs[,3],(1:n.obs)/(n.obs+1)) %>% as.numeric()

    if(isTRUE(analytical)){
      ncp.w = x$ncp.analytical[1]
      ncp.l = x$ncp.analytical[2]
      ncp.s = x$ncp.analytical[3]
    } else {
      ncp.w = x$ncp.simbased[1]
      ncp.l = x$ncp.simbased[2]
      ncp.s = x$ncp.simbased[3]
    }
    exp.w = qchisq((1:n.obs)/(n.obs+1),df = df, ncp = ncp.w) %>% as.numeric()
    exp.l = qchisq((1:n.obs)/(n.obs+1),df = df, ncp = ncp.l) %>% as.numeric()
    exp.s = qchisq((1:n.obs)/(n.obs+1),df = df, ncp = ncp.s) %>% as.numeric()

    ks.w = ks.test(obs.w, "pchisq",df=df,ncp=ncp.w)$p.value
    ks.l = ks.test(obs.l, "pchisq",df=df,ncp=ncp.l)$p.value
    ks.s = ks.test(obs.s, "pchisq",df=df,ncp=ncp.s)$p.value

    if(is.na(ks.w)){ks.w =0.10}
    if(is.na(ks.l)){ks.l =0.10}
    if(is.na(ks.s)){ks.s =0.10}

    df = data.frame(obs=c(obs.w,obs.l,obs.s),exp=c(exp.w,exp.l,exp.s),kind = factor(c(rep("Wald",length(obs.w)),rep("LR",length(obs.l)),rep("Score",length(obs.s))),levels=c("Wald","LR","Score")),n = x$n.pers,esize = x$esize)
    df2 = data.frame(ks=c(ks.w,ks.l,ks.s),kind = factor(c("Wald","LR","Score")),n = x$n.pers,esize = x$esize)
    return(list(df,df2))
  }

  # Get the indices of the desired condition
  indices = which(sapply(res,function(x) x$type == type & x$n.items == n.items & x$dist == dist ))
  resx = res[indices]

  #Setup the data and ks tests
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
  esize.labs <- c("H0: No effect","H1: Small effect","H1: Large effect")
  names(esize.labs) <- c("no", "small","large")

  # Generate Title
  if(type=="1PLvs2PL") {t.type = "Rasch vs 2PL"} else {"t.type" = "DIF"}
  if(isTRUE(analytical)){t.sim="analytical power"}else{t.sim="simulation-based power"}
  if(dist=="unif"){t.dist="uniform distribution"}else{t.dist="skewed normal distribution"}
  title = paste(t.type,", ",n.items," items, ",t.sim,", ",t.dist,sep="")

  re = ggplot(data, aes(x=exp, y=obs,group=kind)) +
    geom_point(aes(color=kind,shape=kind),alpha = 1) +
    xlab("Expected Quantiles")+  ylab("Observed Quantiles") + geom_abline(linetype="dashed", color = "black") + theme(legend.position = "bottom") + facet_grid(vars(esize),vars(n),scale="free", labeller = labeller(n = n.labs, esize = esize.labs)) +
    labs(colour = "Statistic",shape = "Statistic") +
    theme_bw() + ggtitle(title)
    # geom_text(dataks,mapping=aes(x=Inf,y=-Inf,label=sig),hjust=1.2,vjust=-1 )


  return(re)
}



ks.table = function(res,type="DIF2PL",n.items=10,analytical = TRUE) {
  #COLS: Hypothesis, Effect Size, Sample Size, Wald(D,p), LR(D,p), Score(D,p)
  # *significant after bonferroni correction

  setup.data = function(x) {

    df = x$hyp$resmod$Amat %>% nrow()

    n.obs = nrow(x$obs)

    obs.w = quantile(x$obs[,1],(1:n.obs)/(n.obs+1)) %>% as.numeric()
    obs.l = quantile(x$obs[,2],(1:n.obs)/(n.obs+1)) %>% as.numeric()
    obs.s = quantile(x$obs[,3],(1:n.obs)/(n.obs+1)) %>% as.numeric()

    if(isTRUE(analytical)){
      ncp.w = x$ncp.analytical[1]
      ncp.l = x$ncp.analytical[2]
      ncp.s = x$ncp.analytical[3]
    } else {
      ncp.w = x$ncp.simbased[1]
      ncp.l = x$ncp.simbased[2]
      ncp.s = x$ncp.simbased[3]
    }

    ks.wx = ks.test(x$obs[,1], "pchisq",df=df,ncp=ncp.w)
    ks.lx = ks.test(x$obs[,2], "pchisq",df=df,ncp=ncp.l)
    ks.sx = ks.test(x$obs[,3], "pchisq",df=df,ncp=ncp.s)

# if(is.na(ks.wx$p.value)) {browser()}
    ds = format(round(c(ks.wx$statistic,ks.lx$statistic,ks.sx$statistic),3),nsmall=3)
    ps = paste(format(round(c(ks.wx$p.value,ks.lx$p.value,ks.sx$p.value),3),nsmall=3),c(if(ks.wx$p.value*27<.05) {"*"}else{" "},if(ks.lx$p.value*27<.05) {"*"}else{" "},if(ks.sx$p.value*27<.05) {"*"}else{" "}),sep="")

    hyp = as.character(x$type)
    if(hyp=="DIF2PL") {hyp="DIF"}
    if(hyp=="1PLvs2PL") {hyp="Rasch vs 2PL"}

    if(isTRUE(analytical)){method="analytical"}else{method="simbased"}

    df = c(
      hyp = hyp,
      items = x$n.items,
      esize = as.character(x$esize),
      n.pers = x$n.pers,
      method = method,
      d1 = ds[1],p1 = ps[1],d1 = ds[2],p1 = ps[2],d1 = ds[3],p1 = ps[3]
    )

    names(df) = c("Hypothesis", "Items","Effect size", "Sample size","Method","D","p","D","p","D","p")
    return(df)
  }

  indices = which(sapply(res,function(x) x$type %in% type & x$n.items %in% n.items ))
  resx = res[indices]

  # for (i in 1:length(res.stat)) {
  #   res.stat[[i]]$ncp = res.ncp[[i]]
  # }

  # resx = res.stat

  re = as.data.frame(do.call(rbind,lapply(resx,FUN=setup.data)))

  return(re)
}


power.plot = function(res,analytical = TRUE,alpha=.05) {

  setup.data = function(x) {

    df = x$hyp$resmod$Amat %>% nrow()
    crit = qchisq(1-alpha,df = df, ncp = 0) %>% as.numeric()

    pwro.w = mean(x$obs[,1]>crit)
    pwro.l = mean(x$obs[,2]>crit)
    pwro.s = mean(x$obs[,3]>crit)

    if(isTRUE(analytical)){
      ncp.w = x$ncp.analytical[1]
      ncp.l = x$ncp.analytical[2]
      ncp.s = x$ncp.analytical[3]
    } else {
      ncp.w = x$ncp.simbased[1]
      ncp.l = x$ncp.simbased[2]
      ncp.s = x$ncp.simbased[3]
    }

    pwre.w = 1 - pchisq(q = crit, df = df, ncp = ncp.w)
    pwre.l = 1 - pchisq(q = crit, df = df, ncp = ncp.l)
    pwre.s = 1 - pchisq(q = crit, df = df, ncp = ncp.s)

    re = data.frame(
      stat = factor(c("Wald","LR","Score"),levels = c("Wald","LR","Score")),
      n.items = x$n.items,
      type = x$type,
      pw.exp= c(pwre.w ,pwre.l,pwre.s),
      pw.obs= c(pwro.w ,pwro.l,pwro.s))

    return(re)
  }

  # Get the indices of the desired condition and subset results accordingly
  if (isTRUE(analytical)) {
    indices = which(sapply(res,function(x) x$n.items == 10 ))
  } else {
    indices = 1:length(res)
  }
  resx = res[indices]

  # Calculate Observed and expected power
  dat = lapply(resx,function(x) setup.data(x)) %>% do.call(rbind,.) %>% as.data.frame()

  #Prepare Confidence Intervals
  p = (1:10000)/(10000+1)
  q = 1-p
  nruns = res[[1]]$simruns
  ci = sqrt(nruns*p*q)/nruns * qnorm(.975)
  cidat = data.frame(p=p,pl=p-ci,ph=p+ci)
  cidat[cidat<0] = 0
  cidat[cidat>1] = 1


  # dat$pw.exp[is.na(dat$pw.exp)] =.1 # remove for real data!


  # Facet Labels
  item.labs <- paste(unique(sapply(resx,function(x) x$n.items)),"items")
  names(item.labs) <- as.character(unique(sapply(resx,function(x) x$n.items)))
  type.labs <- c("Rasch vs 2PL", "DIF")
  names(type.labs) <- c("1PLvs2PL", "DIF2PL")

  # Generate Title
  if(isTRUE(analytical)){title="Analytical power"}else{title="Simulation-based power"}

  re = ggplot() +
    geom_ribbon(data = cidat, aes(p,p,ymin = pl, ymax = ph), fill = "grey90") +
    geom_abline(linetype="solid", color = "black") +
    geom_point(data = dat, aes(pw.exp, pw.obs, colour = stat,shape=stat)) +
    # geom_errorbar(data=dat, mapping= aes(ymin = pw.lower1, ymax = pw.upper1),width=.05) +
    # scale_color_grey(start=.0,end=.6) +
    labs(colour = "Statistic",shape = "Statistic") +
    theme_bw() + scale_y_continuous("Observed Hit Rate",breaks=seq(0,1,.2),limits = c(0, 1)) + scale_x_continuous("Expected Hit Rate",breaks=seq(0,1,.2),limits = c(0, 1)) +
    facet_grid(vars(n.items),vars(type),scale="free", labeller = labeller(n.items = item.labs, type = type.labs)) + ggtitle(title)

  return(re)
}

power.plot.mis = function(res,analytical = TRUE,alpha=.05,dist ="unif") {

  setup.data = function(x) {

    df = x$hyp$resmod$Amat %>% nrow()
    crit = qchisq(1-alpha,df = df, ncp = 0) %>% as.numeric()

    pwro.w = mean(x$obs[,1]>crit)
    pwro.l = mean(x$obs[,2]>crit)
    pwro.s = mean(x$obs[,3]>crit)

    if(isTRUE(analytical)){
      ncp.w = x$ncp.analytical[1]
      ncp.l = x$ncp.analytical[2]
      ncp.s = x$ncp.analytical[3]
    } else {
      ncp.w = x$ncp.simbased[1]
      ncp.l = x$ncp.simbased[2]
      ncp.s = x$ncp.simbased[3]
    }

    pwre.w = 1 - pchisq(q = crit, df = df, ncp = ncp.w)
    pwre.l = 1 - pchisq(q = crit, df = df, ncp = ncp.l)
    pwre.s = 1 - pchisq(q = crit, df = df, ncp = ncp.s)

    re = data.frame(
      stat = factor(c("Wald","LR","Score"),levels = c("Wald","LR","Score")),
      n.items = x$n.items,
      type = x$type,
      pw.exp= c(pwre.w ,pwre.l,pwre.s),
      pw.obs= c(pwro.w ,pwro.l,pwro.s))

    return(re)
  }

  # Get the indices of the desired condition and subset results accordingly
  if (isTRUE(analytical)) {
    indices = which(sapply(res,function(x) x$n.items == 10 & x$dist==dist))
  } else {
    indices = 1:length(res)
  }
  resx = res[indices]

  # Calculate Observed and expected power
  dat = lapply(resx,function(x) setup.data(x)) %>% do.call(rbind,.) %>% as.data.frame()

  #Prepare Confidence Intervals
  p = (1:10000)/(10000+1)
  q = 1-p
  nruns = res[[1]]$simruns
  ci = sqrt(nruns*p*q)/nruns * qnorm(.975)
  cidat = data.frame(p=p,pl=p-ci,ph=p+ci)
  cidat[cidat<0] = 0
  cidat[cidat>1] = 1


  # dat$pw.exp[is.na(dat$pw.exp)] =.1 # remove for real data!


  # Facet Labels
  item.labs <- paste(unique(sapply(resx,function(x) x$n.items)),"items")
  names(item.labs) <- as.character(unique(sapply(resx,function(x) x$n.items)))
  type.labs <- c("Rasch vs 2PL", "DIF")
  names(type.labs) <- c("1PLvs2PL", "DIF2PL")

  # Generate Title
  if(isTRUE(analytical)){title="Analytical power"}else{title="Simulation-based power"}
  if(dist=="unif"){t.dist="uniform distribution"}else{t.dist="skewed normal distribution"}
  title = paste(title,", ",t.dist,sep="")

  re = ggplot() +
    geom_ribbon(data = cidat, aes(p,p,ymin = pl, ymax = ph), fill = "grey90") +
    geom_abline(linetype="solid", color = "black") +
    geom_point(data = dat, aes(pw.exp, pw.obs, colour = stat,shape=stat)) +
    # geom_errorbar(data=dat, mapping= aes(ymin = pw.lower1, ymax = pw.upper1),width=.05) +
    # scale_color_grey(start=.0,end=.6) +
    labs(colour = "Statistic",shape = "Statistic") +
    theme_bw() + scale_y_continuous("Observed Hit Rate",breaks=seq(0,1,.2),limits = c(0, 1)) + scale_x_continuous("Expected Hit Rate",breaks=seq(0,1,.2),limits = c(0, 1)) +
    facet_grid(vars(n.items),vars(type),scale="free", labeller = labeller(n.items = item.labs, type = type.labs)) + ggtitle(title)

  return(re)
}



curve.plot.multi = function(ncpslist,hyp,alpha=.05,nolegend=FALSE,scales="free") {
# Plot the expected power by sample size for 2 sets of ncps

  ncps = ncpslist[[1]]

  # upper =  uniroot(f = function(x) .999 - power(df=df,ncp=ncps[1],alpha=alpha,ssize=x),interval = c(0, 1000000),tol = .Machine$double.eps ^ 0.5)$root

  upper =  uniroot(f = function(x) .999 - power(hyp=hyp,ncp=ncps[1],alpha=alpha,ssize=x),interval = c(0, 1000000),tol = .Machine$double.eps ^ 0.5)$root

  ssizes = seq(0,upper,10)
  res = list()
  for (i in 1:length(ssizes)) {
    res[[i]] = sapply(ncps,function (x) power(hyp=hyp,ncp=x,alpha=alpha,ssize=ssizes[i]))
  }

  dat = data.frame(do.call(rbind,res))

  names(dat) = c("Wald","LR","Score")
  dat$ssize = ssizes
  dat1 = as.data.frame(pivot_longer(dat,cols=names(dat)[1:3]))

  ncps = ncpslist[[2]]

  upper =  uniroot(f = function(x) .999 - power(hyp=hyp,ncp=ncps[1],alpha=alpha,ssize=x),interval = c(0, 1000000),tol = .Machine$double.eps ^ 0.5)$root

  ssizes = seq(0,upper,10)
  res = list()
  for (i in 1:length(ssizes)) {
    res[[i]] = sapply(ncps,function (x) power(hyp=hyp,ncp=x,alpha=alpha,ssize=ssizes[i]))
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
    geom_line(aes(color=name),size=1) +
    # scale_color_grey(start=.0,end=.6) +
    labs(title = "")+  xlab("Sample size")+ ylab("Power") +labs(color = "Statistic")  + facet_wrap(~set,scales = scales) + scale_y_continuous("Power",breaks=seq(0,1,.2)) + theme_bw()


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
    ) +
    # scale_fill_manual("",breaks=c("Null","Wald","LR","Score"),values = c("grey90","black","grey30","grey60")) +
    scale_fill_manual("",breaks=c("Null","Wald","LR","Score"),values = c("grey90","#F8766D","#7CAE00","#00BFC4")) +
    theme_bw()

  return(p)
}


chisq.dens.null = function(ncp,df,n,alpha=.05,lim=30) {
  # Density of null hypothesis only
  crit= qchisq(1-alpha,df)

  p = ggplot(data.frame(x = c(0, lim)), aes(x)) +
    stat_function(
      fun = dchisq,
      geom = "area",
      # alpha = .5,
      args = list(
        df=df
      ),n = 5000,aes(fill = "Null")
    )  +
    stat_function(
      fun = dchisq,
      geom = "area",
      xlim = c(crit,lim),
      # alpha = .5,
      args = list(
        df=df
      ),n = 5000,aes(fill = "One")
    )  +
    stat_function(
      fun = dchisq,
      geom = "line",
      args = list(
        df=df
      ),n = 5000
    )+

    labs(
      title = "",
      x = "Statistic",
      y = "Density"
    ) +
    # scale_fill_manual("",breaks=c("Null","Wald","LR","Score"),values = c("grey90","black","grey30","grey60")) +
    scale_fill_manual("",breaks=c("Null","One"),values = c("grey90","black")) +
    theme_bw() + theme(legend.position = "none")

  return(p)
}



power.table = function(res,analytical = TRUE,alpha=.05) {

  setup.data = function(x) {

    df = x$hyp$resmod$Amat %>% nrow()
    crit = qchisq(1-alpha,df = df, ncp = 0) %>% as.numeric()

    pwro.w = mean(x$obs[,1]>crit)
    pwro.l = mean(x$obs[,2]>crit)
    pwro.s = mean(x$obs[,3]>crit)

    if(isTRUE(analytical)){
      ncp.w = x$ncp.analytical[1]
      ncp.l = x$ncp.analytical[2]
      ncp.s = x$ncp.analytical[3]
    } else {
      ncp.w = x$ncp.simbased[1]
      ncp.l = x$ncp.simbased[2]
      ncp.s = x$ncp.simbased[3]
    }

    pwre.w = 1 - pchisq(q = crit, df = df, ncp = ncp.w)
    pwre.l = 1 - pchisq(q = crit, df = df, ncp = ncp.l)
    pwre.s = 1 - pchisq(q = crit, df = df, ncp = ncp.s)

    envelopes = sapply(c(pwre.w ,pwre.l,pwre.s),function(x) {
      p = x
      q = 1-p
      nruns = res[[1]]$simruns
      ci = sqrt(nruns*p*q)/nruns * qnorm(.975)
      if(p-ci<0){lower = 0} else {lower = p-ci}
      if(p+ci>1){higher = 1} else {higher = p+ci}
      string = paste("[",round(lower,3),", ",round(higher,3),"]",sep="")
      return(string)
    })

    hyp = as.character(x$type)
    if(hyp=="DIF2PL") {hyp="DIF"}
    if(hyp=="1PLvs2PL") {hyp="Rasch vs 2PL"}

    if(isTRUE(analytical)){method="analytical"}else{method="simbased"}

    df = data.frame(
      hyp = hyp,
      items = as.character(x$n.items),
      esize = as.character(x$esize),
      n.pers = as.character(x$n.pers),
      method = method,
      pw.obs1= pwro.w %>% round(.,3),
      pw.exp1= pwre.w %>% round(.,3),
      env1 = envelopes[1],
      pw.obs2= pwro.l %>% round(.,3),
      pw.exp2= pwre.l %>% round(.,3),
      env2 = envelopes[2],
      pw.obs3= pwro.s %>% round(.,3),
      pw.exp3= pwre.s %>% round(.,3),
      env3 = envelopes[3])

    names(df) = c("Hypothesis", "Items","Effect size", "Sample size","Method","OHR","EHR","Env","OHR","EHR","Env","OHR","EHR","Env")

    return(df)
  }

  # Get the indices of the desired condition and subset results accordingly
  if (isTRUE(analytical)) {
    indices = which(sapply(res,function(x) x$n.items == 10 ))
  } else {
    indices = 1:length(res)
  }
  resx = res[indices]

  # Calculate Observed and expected power
  dat = lapply(resx,function(x) setup.data(x)) %>% do.call(rbind,.) %>% as.data.frame()


  return(dat)
}



power.table.mis = function(res,analytical = TRUE,alpha=.05,dist=NULL) {

  setup.data = function(x) {

    df = x$hyp$resmod$Amat %>% nrow()
    crit = qchisq(1-alpha,df = df, ncp = 0) %>% as.numeric()

    pwro.w = mean(x$obs[,1]>crit)
    pwro.l = mean(x$obs[,2]>crit)
    pwro.s = mean(x$obs[,3]>crit)

    if(isTRUE(analytical)){
      ncp.w = x$ncp.analytical[1]
      ncp.l = x$ncp.analytical[2]
      ncp.s = x$ncp.analytical[3]
    } else {
      ncp.w = x$ncp.simbased[1]
      ncp.l = x$ncp.simbased[2]
      ncp.s = x$ncp.simbased[3]
    }

    pwre.w = 1 - pchisq(q = crit, df = df, ncp = ncp.w)
    pwre.l = 1 - pchisq(q = crit, df = df, ncp = ncp.l)
    pwre.s = 1 - pchisq(q = crit, df = df, ncp = ncp.s)

    envelopes = sapply(c(pwre.w ,pwre.l,pwre.s),function(x) {
      p = x
      q = 1-p
      nruns = res[[1]]$simruns
      ci = sqrt(nruns*p*q)/nruns * qnorm(.975)
      if(p-ci<0){lower = 0} else {lower = p-ci}
      if(p+ci>1){higher = 1} else {higher = p+ci}
      string = paste("[",round(lower,3),", ",round(higher,3),"]",sep="")
      return(string)
    })

    hyp = as.character(x$type)
    if(hyp=="DIF2PL") {hyp="DIF"}
    if(hyp=="1PLvs2PL") {hyp="Rasch vs 2PL"}

    if(isTRUE(analytical)){method="analytical"}else{method="simbased"}

    if(dist=="unif"){t.dist="uniform"}else{t.dist="skewed normal"}

    df = data.frame(
      esize = as.character(x$esize),
      n.pers = as.character(x$n.pers),
      dist = t.dist,
      pw.obs1= pwro.w %>% round(.,3),
      pw.exp1= pwre.w %>% round(.,3),
      env1 = envelopes[1],
      pw.obs2= pwro.l %>% round(.,3),
      pw.exp2= pwre.l %>% round(.,3),
      env2 = envelopes[2],
      pw.obs3= pwro.s %>% round(.,3),
      pw.exp3= pwre.s %>% round(.,3),
      env3 = envelopes[3])

    names(df) = c("Effect size", "Sample size","Distribution","OHR","EHR","Env","OHR","EHR","Env","OHR","EHR","Env")

    return(df)
  }

  # Get the indices of the desired condition and subset results accordingly
  if (isTRUE(analytical)) {
    indices = which(sapply(res,function(x) x$n.items == 10 & x$dist==dist))
  } else {
    indices = 1:length(res)
  }
  resx = res[indices]


  # Calculate Observed and expected power
  dat = lapply(resx,function(x) setup.data(x)) %>% do.call(rbind,.) %>% as.data.frame()


  return(dat)
}


# power.table2 = function(res,analytical = TRUE,alpha=.05) {
#
#   setup.data = function(x) {
#
#     df = x$hyp$resmod$Amat %>% nrow()
#     crit = qchisq(1-alpha,df = df, ncp = 0) %>% as.numeric()
#
#     pwro.w = mean(x$obs[,1]>crit)
#     pwro.l = mean(x$obs[,2]>crit)
#     pwro.s = mean(x$obs[,3]>crit)
#
#     if(isTRUE(analytical)){
#       ncp.w = x$ncp.analytical[1]
#       ncp.l = x$ncp.analytical[2]
#       ncp.s = x$ncp.analytical[3]
#     } else {
#       ncp.w = x$ncp.simbased[1]
#       ncp.l = x$ncp.simbased[2]
#       ncp.s = x$ncp.simbased[3]
#     }
#
#     pwre.w = 1 - pchisq(q = crit, df = df, ncp = ncp.w)
#     pwre.l = 1 - pchisq(q = crit, df = df, ncp = ncp.l)
#     pwre.s = 1 - pchisq(q = crit, df = df, ncp = ncp.s)
#
#     envelopes = sapply(c(pwre.w ,pwre.l,pwre.s),function(x) {
#       p = x
#       q = 1-p
#       nruns = res[[1]]$simruns
#       ci = sqrt(nruns*p*q)/nruns * qnorm(.975)
#       if(p-ci<0){lower = 0} else {lower = p-ci}
#       if(p+ci>1){higher = 1} else {higher = p+ci}
#       string = paste("[",round(lower,3),", ",round(higher,3),"]",sep="")
#       return(string)
#     })
#
#     hyp = as.character(x$type)
#     if(hyp=="DIF2PL") {hyp="DIF"}
#     if(hyp=="1PLvs2PL") {hyp="Rasch vs 2PL"}
#
#     if(isTRUE(analytical)){method="analytical"}else{method="simbased"}
#
#     df = data.frame(
#       hyp = hyp,
#       items = x$n.items,
#       esize = as.character(x$esize),
#       n.pers = x$n.pers,
#       method = method,
#       stat = factor(c("Wald","LR","Score"),levels = c("Wald","LR","Score")),
#       pw.obs= c(pwro.w ,pwro.l,pwro.s) %>% sapply(.,function(x) round(x,3)),
#       pw.exp= c(pwre.w ,pwre.l,pwre.s) %>% sapply(.,function(x) round(x,3)) ,
#       envelopes = envelopes)
#
#     names(df) = c("Hypothesis", "Items","Effect size", "Sample size","Method","Statistic","Observed Hit Rate","Expected Hit Rate","Envelope")
#
#     return(df)
#   }
#
#   # Get the indices of the desired condition and subset results accordingly
#   if (isTRUE(analytical)) {
#     indices = which(sapply(res,function(x) x$n.items == 10 ))
#   } else {
#     indices = 1:length(res)
#   }
#   resx = res[indices]
#
#   # Calculate Observed and expected power
#   dat = lapply(resx,function(x) setup.data(x)) %>% do.call(rbind,.) %>% as.data.frame()
#
#
#   return(dat)
# }
#
