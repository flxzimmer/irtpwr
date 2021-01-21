
# 2PL, 1PL, DIF Figure


library(ggplot2)
library(gridExtra)
library(grid)

pl2 = function(x,a,d) {
  1/(1+exp(-(a*x+d)))
}


# 2PL  --------------------------------------------------------------------

colors <- c("1" = "green", "2" = "red", "3" = "orange","4"="blue","5" = "black")

p1 = ggplot(data.frame(x = c(-4, 4)), aes(x)) +
  stat_function(
    fun = pl2,
    geom = "line",
    args = list(
      a = 1,
      d = 0
    ),aes(color="1")) +
  stat_function(
    fun = pl2,
    geom = "line",
    args = list(
      a = .7,
      d = .5
    ),aes(color="2")) +
  stat_function(
    fun = pl2,
    geom = "line",
    args = list(
      a = 1.3,
      d = -.5
    ),aes(color="3")) +
  stat_function(
    fun = pl2,
    geom = "line",
    args = list(
      a = 1.1,
      d = -1.5
    ),aes(color="4")) +
  stat_function(
    fun = pl2,
    geom = "line",
    args = list(
      a = 1.1,
      d = .2
    ),aes(color="5")) +
  labs(
    title = "",
    x = expression(theta),
    y = "Density",
    color = "Items"
  ) +
  scale_color_manual(values = colors) +
  # ggtitle("2PL")+
  theme_bw()


# 1PL ---------------------------------------------------------------------

colors <- c("1" = "green", "2" = "red", "3" = "orange","4"="blue","5" = "black")

p2 = ggplot(data.frame(x = c(-4, 4)), aes(x)) +
  stat_function(
    fun = pl2,
    geom = "line",
    args = list(
      a = 1,
      d = 0
    ),aes(color="1")) +
  stat_function(
    fun = pl2,
    geom = "line",
    args = list(
      a = 1,
      d = .5
    ),aes(color="2")) +
  stat_function(
    fun = pl2,
    geom = "line",
    args = list(
      a = 1,
      d = -.5
    ),aes(color="3")) +
  stat_function(
    fun = pl2,
    geom = "line",
    args = list(
      a = 1,
      d = -1.5
    ),aes(color="4")) +
  stat_function(
    fun = pl2,
    geom = "line",
    args = list(
      a = 1,
      d = .2
    ),aes(color="5")) +
  labs(
    title = "",
    x = expression(theta),
    y = "Density",
    color = "Items"
  ) +
  scale_color_manual(values = colors) +
  # ggtitle("1PL")+
  theme_bw()


# DIF ---------------------------------------------------------------------


colors <- c("1" = "gray", "2" = "gray", "3" = "gray","4"="blue","5" = "gray")

p3 = ggplot(data.frame(x = c(-4, 4)), aes(x)) +
  stat_function(
    fun = pl2,
    geom = "line",
    args = list(
      a = 1,
      d = 0
    ),aes(color="1")) +
  stat_function(
    fun = pl2,
    geom = "line",
    args = list(
      a = .7,
      d = .5
    ),aes(color="2")) +
  stat_function(
    fun = pl2,
    geom = "line",
    args = list(
      a = 1.3,
      d = -.5
    ),aes(color="3")) +
  stat_function(
    fun = pl2,
    geom = "line",
    args = list(
      a = 1.1,
      d = -1.5
    ),aes(color="4")) +
  stat_function(
    fun = pl2,
    geom = "line",
    args = list(
      a = 1.1,
      d = .2
    ),aes(color="5")) +
  labs(
    title = "",
    x = expression(theta),
    y = "Density"
  ) +
  scale_color_manual(values = colors) +
  theme_bw()+
  ggtitle("Group A") +
  theme(legend.position = "none")

p4 = ggplot(data.frame(x = c(-4, 4)), aes(x)) +
  stat_function(
    fun = pl2,
    geom = "line",
    args = list(
      a = 1,
      d = 0
    ),aes(color="1")) +
  stat_function(
    fun = pl2,
    geom = "line",
    args = list(
      a = .7,
      d = .5
    ),aes(color="2")) +
  stat_function(
    fun = pl2,
    geom = "line",
    args = list(
      a = 1.3,
      d = -.5
    ),aes(color="3")) +
  stat_function(
    fun = pl2,
    geom = "line",
    args = list(
      a = 1.1,
      d = 1.5
    ),aes(color="4")) +
  stat_function(
    fun = pl2,
    geom = "line",
    args = list(
      a = 1.1,
      d = .2
    ),aes(color="5")) +
  labs(
    title = "",
    x = expression(theta),
    y = "",
    color = "Item") +
  scale_color_manual(values = colors) +
  theme_bw() +
  ggtitle("Group B")

p5 = grid.arrange(p3,p4,ncol=2,widths = c(.46,.54))


# export ------------------------------------------------------------------


pdf("fig1.pdf",height=6,width=8);p1;dev.off()
pdf("fig2.pdf",height=6,width=8);p2;dev.off()
pdf("fig3.pdf",height=6,width=10);grid.draw(p5);dev.off()



#
#
# set.seed(1234)
# altpars <- list(
#   a = c(.5,1.4,.8,1,1.5),
#   d = c(-2,1,0,1,2)
# )
#
# hyp <- setup_hypothesis(type = "1PLvs2PL", altpars = altpars)
#
# data <- setup.data(hyp=hyp,n=5000)$data
#
#
# mod1 = mirt(data,1,itemtype = c('Rasch'),verbose = FALSE)
#
# mod2 = mirt(data,1,itemtype = c('2PL'),verbose = FALSE)
#
# plot(mod1,drop2 = TRUE,type="trace",facet_items=FALSE,main="")
#
# plot(mod2,drop2 = TRUE,type="trace",facet_items=FALSE,main="")
#
# # DIF Figure
#
#
#
#
#
#
#
#
# # modx = mirt(Science,1,itemtype = c('Rasch'),verbose = FALSE)
# # plot(modx,type="trace",facet_items=FALSE,main="")
#
#
# # dat6 <- expand.table(LSAT6)
# # # Estimate pars
# # pars6 = mirt(dat6,1,verbose = FALSE)
# # plot(pars6 ,type="trace",facet_items=FALSE,main="")
# set.seed(1234)
# group <- sample(c('g1','g2'), nrow(Science), TRUE)
# x2 <- multipleGroup(Science, 1, group)
# plot(x2)
# plot(x2, type = 'trace')
# plot(x2, type = 'trace', which.items = 1:2)
# plot(x2, type = 'itemscore', which.items = 1:2)
# plot(x2, type = 'trace', which.items = 1, facet_items = FALSE) #facet by group
# plot(x2, type = 'info')
#
#
#
#
#
#
#
#
# # DIF
# group1 = group2 <- list(
#   a = rlnorm(5,sdlog = .2),
#   d = rnorm(5)
# )
#
# group2$a[1] = (group2$a[1])^2
# group2$d[1] = group2$d[1] + .5
#
# altpars <- list(group1,group2)
#
# hyp <- setup_hypothesis(type = "DIF2PL", altpars = altpars)
#
# data <- setup.data(hyp=hyp,n=500)
#
#
#
# mod1 = mirt(dat7,1,itemtype = c('Rasch'),verbose = FALSE)
#
# mod2 = mirt(dat7,1,itemtype = c('2PL'),verbose = FALSE)
#
# p1 = plot(mod1,type="trace",facet_items=FALSE,main="")
# p2 = plot(mod2,type="trace",facet_items=FALSE,main="")
#
#
# pdf("fig1.pdf",height=6,width=8);p1;dev.off()
#
# pdf("fig2.pdf",height=6,width=8);p2;dev.off()
