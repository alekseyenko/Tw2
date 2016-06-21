sr = read.table("../results/simulation_results2.txt", sep="\t", header=T)

params = expand.grid(n1 = c(5,10,15,20), n2 = c(5,10,15,20), fsd = c(1, .8, .6, .4, .2), eff = 0:5)

res = data.frame()
for(i in 1:(dim(params)[1])){
  adonis.pos = sum(subset(sr, 
                          n1 == params[i,]$n1 & 
                            n2 == params[i,]$n2 &
                            fracSd == params[i,]$fsd &
                            effect == params[i,]$eff)$adonis.p < 0.05)
  welchT.pos = sum(subset(sr, 
                          n1 == params[i,]$n1 & 
                            n2 == params[i,]$n2 &
                            fracSd == params[i,]$fsd &
                            effect == params[i,]$eff)$welchT.p < 0.05)
  welchTS.pos = sum(subset(sr, 
                          n1 == params[i,]$n1 & 
                            n2 == params[i,]$n2 &
                            fracSd == params[i,]$fsd &
                            effect == params[i,]$eff)$welchTS.p < 0.05)
  res = rbind(res,
              c(params[i,], adonis.pos=adonis.pos/1000, welchT.pos=welchT.pos/1000, welchTS.pos=welchTS.pos/1000))
}

res = data.frame(res)
colnames(res) = c("n1", "n2", "Frac.S.D.", "Effect", "PERMANOVA", "Tw2", "welchTS.pos")
#subset(res, Frac.S.D.==1 & Effect==0)

library(ggplot2)
library(reshape2)
library(Hmisc)

pdf('../results/simulation_heatmaps.pdf', width=12, height=10)
ggplot(data = subset(res), aes(x = n1, y = n2)) + 
  geom_tile(aes(fill=adonis.pos), colour = "white") + scale_fill_gradient(low = "green", high = "red") + facet_grid(Frac.S.D.~Effect)

ggplot(data = subset(res), aes(x = n1, y = n2)) + 
  geom_tile(aes(fill=welchT.pos), colour = "white") + scale_fill_gradient(low = "green", high = "red") + facet_grid(Frac.S.D.~Effect)

ggplot(data = subset(res), aes(x = n1, y = n2)) + 
  geom_tile(aes(fill=welchTS.pos), colour = "white") + scale_fill_gradient(low = "green", high = "red") + facet_grid(Frac.S.D.~Effect)
dev.off()

#The following will reproduce results of Figure 1
pdf('../results/adonis_0_8.pdf', width=7, height=5)
ppall = melt(subset(res, Frac.S.D. == .8 & Effect %in% c(0, 2, 4, 5)), id.vars = c('n1', 'n2', 'Frac.S.D.', 'Effect'))
ppall = subset(ppall, variable == "PERMANOVA")
ppall = cbind(ppall, binconf(ppall$value*1000, 1000, method='exact'))
ppall$balanced = rep("No", nrow(ppall))
ppall$balanced[ppall$n1 == ppall$n2] = "Yes"
ppall$balanced = factor(ppall$balanced)
ggplot(data = ppall, aes(x = n1, y = value))+ 
  geom_point(aes(size=n2, shape = balanced)) + geom_line(aes(group=n2))+ 
  facet_grid(. ~ Effect, labeller = label_both) + geom_hline(yintercept=0.05, lty="dotted") +
  ylab("Fraction of rejected null hypotheses") + xlab("Observations in most dispersed sample") + 
  scale_size_continuous(name="Observations\nin least dispersed\nsample") + 
  scale_shape(name = "Balanced sample size")+
  ggtitle(expression(paste("Empirical type I error and power (", alpha, " = 0.05)"))) + theme_bw()
dev.off()


# The following provides the data for Table 1

ppall = melt(subset(res, Frac.S.D. == 1 & Effect == 0), id.vars = c('n1', 'n2', 'Frac.S.D.', 'Effect'))
ppall = subset(ppall, variable == "PERMANOVA")
ppall = cbind(ppall, binconf(ppall$value*1000, 1000, method='exact'))
tapply(ppall$PointEst, list(ppall$n1, ppall$n2), function(x) x[1])

ppall = melt(subset(res, Frac.S.D. == 1 & Effect == 2), id.vars = c('n1', 'n2', 'Frac.S.D.', 'Effect'))
ppall = subset(ppall, variable == "PERMANOVA")
ppall = cbind(ppall, binconf(ppall$value*1000, 1000, method='exact'))
tapply(ppall$PointEst, list(ppall$n1, ppall$n2), function(x) x[1])


# The following reproduces Figure 2
pdf('../results/adonis_vs_Tw.pdf', width=7, height=7)
ppall = melt(subset(res, Frac.S.D. %in% c(.2, .8, 1) & Effect %in% c(0, 4, 5)), id.vars = c('n1', 'n2', 'Frac.S.D.', 'Effect'))
ppall = subset(ppall, variable %in% c("PERMANOVA", "Tw2"))
ppall = cbind(ppall, binconf(ppall$value*1000, 1000, method='exact'))
ppall$balanced = rep("No", nrow(ppall))
ppall$balanced[ppall$n1 == ppall$n2] = "Yes"
ppall$balanced = factor(ppall$balanced)
ggplot(data = ppall, aes(x = n1, y = value))+ 
  geom_point(aes(color=variable, size=n2, shape=balanced)) + geom_hline(yintercept=0.05, lty="dotted") +
  geom_line(aes(color=variable, group=interaction(n2, variable)))+ facet_grid(Frac.S.D. ~ Effect, labeller = label_both) + 
  ylab("Fraction of rejected null hypotheses") +
  xlab("Observations in most dispersed sample") + scale_size_continuous(name="Observations\nin least dispersed\nsample") + 
  scale_color_manual(name="Statistical\ntest", values = c("black", "red")) + 
  scale_shape(name = "Balanced sample size")+
  ggtitle(expression(paste("Empirical type I error and power (", alpha, " = 0.05)"))) + theme_bw()+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()



