#(с) 2015 Alexander V. Alekseyenko

library(phyloseq)
library(ape)
library(doParallel)
registerDoParallel(cores=3)

file = "../data/otu_table.txt"
map = "../data/PsoriasisMetadata.txt"
all.phy = import_qiime(otufilename = file, mapfilename = map)

omega2 = function(adonis.res){
  f = adonis.res$aov.tab[1,4]
  dferr = adonis.res$aov.tab["Total",1]
  dfeff = adonis.res$aov.tab[1,1]
  (f - 1)/(f+(dferr+1)/dfeff)
}
source('Tw2.R')
compare2 = function(dm, f, nperm=99999){
  ar = adonis(dm ~ f, permutations = nperm)
  wt = WT.test(dm, f, nrep = nperm)
  c(table(f), omega2(ar), dist.cohen.d(dm, f), ar$aov.tab[1,6], wt$p.value)
}

css = subset_samples(all.phy, !is.na(TripletC))
table(sample_data(css)$Status)
length(unique(sample_data(css)$TripletC))

## The following computes the data for Table 3

lesion = subset_samples(css, Status %in% c("Control", "Lesion"))
bray.lesion = distance(lesion, method="bray")
compare2(bray.lesion, sample_data(lesion)$Status)
md = diag(meandist(bray.lesion, sample_data(lesion)$Status))
md/md[1]

normal = subset_samples(css, Status %in% c("Control", "Normal"))
bray.normal = distance(normal, method="bray")
compare2(bray.normal, sample_data(normal)$Status)
md = diag(meandist(bray.normal, sample_data(normal)$Status))
md/md[1]

psoriasis = subset_samples(css, Status %in% c("Lesion", "Normal"))
bray.psoriasis = distance(psoriasis, method="bray")
compare2(bray.psoriasis, sample_data(psoriasis)$Status)
md = diag(meandist(bray.psoriasis, sample_data(psoriasis)$Status))
md/md[1]

sqrt(WT(bray.psoriasis, sample_data(psoriasis)$Status)/(length(sample_data(psoriasis)$Status)))

