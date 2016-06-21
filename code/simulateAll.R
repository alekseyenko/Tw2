source('Tw2.R')
source('simulateDists.R')
library(ade4)
library(vegan)
library(doParallel)
registerDoParallel(cores=8)


test3 = function(group.size1 = 10, group.size2=group.size1, fsd=0.9, effect = 1, make.plot=F){
  z = simulate.normal2(group.size1, group.size2, fsd, effect)
  dm = dist(t(z))
  f = factor(c(rep(1, group.size1), rep(2, group.size2)))
  if(make.plot){
    pco = dudi.pco(dm, scannf=F, full=T)
    s.class(pco$li, f)
  }
  a.p.value = adonis(dm ~ f)$aov.tab[1,6]
  w.p.value = WT.test(dm, f)$p.value
  w.p.value2 = Wstar.test(dm, f)$p.value
  c(n1=group.size1, n2=group.size2, fsd=fsd, effect=effect, a.p.value=a.p.value, w.p.value=w.p.value, w.p.value2=w.p.value2)
}

#test3(group.size1 = 5, group.size2=10, fsd=.2, effect=4.5, make.plot=T)
#test3(group.size1 = 20, group.size2=5, fsd=.2, effect=0, make.plot=T)
#date(); pe = replicate(100, test3(group.size1 = 5, group.size2 = 10, fsd=.2, effect=0, make.plot=F)); date()
#pe = t(pe)
#sum(pe[,"a.p.value"] < 0.05)
#sum(pe[,"w.p.value"] < 0.05)
#sum(pe[,"w.p.value2"] < 0.05)


params = expand.grid(replicate = 1:1000, n1 = c(5,10,15,20), n2 = c(5,10,15,20), fsd = c(1, .8, .6, .4, .2), eff = 0:5)
res = foreach(i=1:(dim(params)[1]), .combine=rbind) %dopar%{
  c(params[i,]$replicate, 
    test3(group.size1 = params[i,]$n1, 
          group.size2=params[i,]$n2, 
          fsd=params[i,]$fsd, 
          effect=params[i,]$eff, make.plot=F))
}
colnames(res) = c("Replicate", "n1", "n2", "fracSd", "effect", "adonis.p", "welchT.p", "welchTS.p")

write.table(res, '../results/simulation_results2.txt', sep='\t')


