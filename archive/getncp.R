


re= c()
re2 = c()
for (i in 1:10000) {x = rchisq(1,df=4,ncp=1000)
;re[i] = get_ncp(x,df=4)$ncp;re2[i] = mean(x)-4}
mean(abs(re-1000));mean(abs(re2-1000))
mean(re);mean(re2)


