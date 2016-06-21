sr = read.table("../results/simulation10_results.txt", sep="\t", header=T)

params = expand.grid(n1 = unique(sr$n1), n2 = unique(sr$n2), fsd = unique(sr$fracSd), eff = unique(sr$effect))

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
library(ggplot2)
library(reshape2)
library(Hmisc)


ppall = melt(res, id.vars = c('n1', 'n2', 'Frac.S.D.', 'Effect'))
ppall = subset(ppall, variable %in% c("PERMANOVA", "Tw2") & n1 %in% 8:12 & n2 %in% 8:12 & Effect %in% c(0, 4, 5))
ppall = cbind(ppall, binconf(ppall$value*1000, 1000, method='exact'))
ppall$balanced = rep("No", nrow(ppall))
ppall$balanced[ppall$n1 == ppall$n2] = "Yes"
ppall$balanced = factor(ppall$balanced)

# The following reproduced Figure 3(a)
pdf('../results/adonis_Tw_10.pdf', width=7, height=7)
ggplot(data = ppall, aes(x = n1, y = value))+ 
  geom_point(aes(color=variable, size=n2, shape=balanced)) +
  geom_line(aes(color=variable, group=interaction(n2, variable)))+ facet_grid(Frac.S.D. ~ Effect, labeller = label_both) + 
  ylab("Fraction of rejected null hypotheses") + xlab("Observations in most dispersed sample") + 
  scale_size_continuous(name="Observations\nin least dispersed\nsample") + 
  scale_color_manual(name="Statistical\ntest", values = c("black", "red")) + 
  scale_shape(name = "Balanced sample size")+
  ggtitle(expression(paste("Empirical type I error and power (", alpha, " = 0.05)"))) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


