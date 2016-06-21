#(с) 2015 Alexander V. Alekseyenko

library(phyloseq)
library(ape)
library(doParallel)
library(vegan)
library(ade4)
registerDoParallel(cores=3)

## Load the data using phyloseq
file = "../data/seqs_otu_table.txt"
map = "../data/obesity_mapping.txt"
all.phy = import_qiime(otufilename = file, mapfilename = map)

## function computes omega square
omega2 = function(adonis.res){
  f = adonis.res$aov.tab[1,4]
  dferr = adonis.res$aov.tab["Total",1]
  dfeff = adonis.res$aov.tab[1,1]
  (f - 1)/(f+(dferr+1)/dfeff)
}

## compute bray curtis distances
bray = distance(all.phy, method="bray")
CvsAbx = rep('Abx', nsamples(all.phy))
CvsAbx[sample_data(all.phy)$Treatment=='C'] = 'C'
CvsAbx = factor(CvsAbx)

bray.pco = dudi.pco(cailliez(bray), scannf=F, nf=2)

# Get % inertia for PC1 and PC2
(bray.pco$eig/sum(bray.pco$eig))[1:2]

# The following reproduces figure 4.
pdf('../results/stat_pcoa.pdf', width=3.5, height=3.5)
s.class(bray.pco$li, interaction(sample_data(all.phy)$Location, CvsAbx), 
        cellipse = 0, cstar=0, col=c('red','orange', 'black', 'gray75'), 
        grid=F, sub="PC1 - 28.6%; PC2 - 8.7%", cpoint=0)
s.class(bray.pco$li, interaction(sample_data(all.phy)$Location, CvsAbx), 
        cellipse = 0, cstar=0, col=c('red','orange', 'black', 'gray75'), add.plot=T,
        grid=F, clabel=0)
dev.off()

source('Tw2.R')

compare2 = function(dm, f, nperm=99999){
  ar = adonis(dm ~ f, permutations = nperm)
  wt = WT.test(dm, f, nrep = nperm)
  c(table(f), omega2(ar), dist.cohen.d(dm, f), ar$aov.tab[1,6], wt$p.value)
}
compare2(bray, sample_data(all.phy)$Location)

cecal = subset_samples(all.phy, Location == 'cecal')
cecal.CvsAbx = rep('Abx', nsamples(cecal))
cecal.CvsAbx[sample_data(cecal)$Treatment=='C'] = 'C'
cecal.CvsAbx = factor(cecal.CvsAbx)

bray.cecal = as.dist(as.matrix(bray)[sample_data(all.phy)$Location =='cecal', 
                                     sample_data(all.phy)$Location =='cecal'])

## The following computes the data for Table 2
res= c('cecal', 'C', 'Abx', compare2(bray.cecal, cecal.CvsAbx))
for(lev in levels(sample_data(cecal)$Treatment)[-1]){
  Tr=sample_data(cecal)$Treatment
  subs = Tr %in% c("C", lev)
  sdm = as.dist(as.matrix(bray.cecal)[subs, subs])
  Tr = as.factor(as.character(Tr[subs]))
  res = rbind(res,
              c('cecal','C', lev, compare2(sdm, Tr)))
}
res

alld = diag(meandist(bray.cecal, cecal.CvsAbx))
alld/alld[2]

alld = diag(meandist(bray.cecal, sample_data(cecal)$Treatment))
alld/alld[1]


fecal = subset_samples(all.phy, Location == 'fecal')
fecal.CvsAbx = rep('Abx', nsamples(fecal))
fecal.CvsAbx[sample_data(fecal)$Treatment=='C'] = 'C'
fecal.CvsAbx = factor(fecal.CvsAbx)

bray.fecal = as.dist(as.matrix(bray)[sample_data(all.phy)$Location =='fecal', 
                                     sample_data(all.phy)$Location =='fecal'])
res= rbind(res,
           c('fecal', 'C', 'Abx', compare2(bray.fecal, fecal.CvsAbx)))

for(lev in levels(sample_data(fecal)$Treatment)[-1]){
  Tr=sample_data(fecal)$Treatment
  subs = Tr %in% c("C", lev)
  sdm = as.dist(as.matrix(bray.fecal)[subs, subs])
  Tr = as.factor(as.character(Tr[subs]))
  res = rbind(res,
              c('fecal', 'C', lev, compare2(sdm, Tr)))
}
res

alld = diag(meandist(bray.fecal, fecal.CvsAbx))
alld/alld[2]

alld = diag(meandist(bray.fecal, sample_data(fecal)$Treatment))
alld/alld[1]

colnames(res)=c("location", "grp1", "grp2", "n1", "n2", "omega2", "d", "permanova.p", "Tw2.p")

write.table(res, '../results/stat_significance.txt', sep="\t" )




