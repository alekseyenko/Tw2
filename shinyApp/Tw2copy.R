is.dist = function(x) any(class(x)=='dist')

WT.mean = function(dm, f){
  if(nlevels(f) != 2) return(NULL)
  lev = levels(f)
  ns = summary(f)
  N = sum(ns)
  dd = as.matrix(dm)
  dd[upper.tri(dd)]=0 ##
  SST = sum(dd^2)/N
  SSW1 = sum(dd[f==lev[1],f==lev[1]]^2)/ns[1] 
  SSW2 = sum(dd[f==lev[2],f==lev[2]]^2)/ns[2]
  SSW = SSW1 + SSW2
  
  (sqrt((ns[1]+ns[2])/(ns[1]*ns[2])))*sqrt(SST-SSW)
}


dist.sigma2 = function(dm){
  dd = as.matrix(dm)
  dd[upper.tri(dd)]=0 ##
  sum(dd^2)/nrow(dm)/(nrow(dm)-1) 
}

dist.cohen.d = function(dm, f){
  mean.diff = sqrt(WT.mean(dm, f))
  lev = levels(f)
  s1 = dist.sigma2(as.matrix(dm)[f==lev[1], f == lev[1]])
  s2 = dist.sigma2(as.matrix(dm)[f==lev[1], f == lev[1]])
  ns=table(f)
  mean.diff/sqrt(((ns[1]-1)*s1 + (ns[2]-1)*s2)/(sum(ns)-2))
}

# dist.cohen.d.strata = function(dm, f, strata){
#   p = unlist(tapply(f,strata,sample))
#   p2 = unlist(tapply(1:length(f), strata, identity))
#   # mean.diff = sqrt(WT.mean(dm, f))
#   mean.diff.strata = sqrt(WT.mean(dm, f[p[p2]]))
#   lev = levels(f)
#   s1 = dist.sigma2(as.matrix(dm)[f==lev[1], f == lev[1]])
#   s2 = dist.sigma2(as.matrix(dm)[f==lev[1], f == lev[1]])
#   ns=table(f)
#   mean.diff.strata/sqrt(((ns[1]-1)*s1 + (ns[2]-1)*s2)/(sum(ns)-2))
# }

WT = function(dm, f){
  if(nlevels(f) != 2) return(NULL)
  lev = levels(f)
  ns = summary(f)
  N = sum(ns)
  dd = as.matrix(dm)
  dd[upper.tri(dd)]=0 ##
  SST = sum(dd^2)/N
  SSW1 = sum(dd[f==lev[1],f==lev[1]]^2)/ns[1] 
  SSW2 = sum(dd[f==lev[2],f==lev[2]]^2)/ns[2]
  SSW = SSW1 + SSW2
  
  s1 = SSW1/(ns[1]-1)
  s2 = SSW2/(ns[2]-1)
  
  t.stat = (sqrt((ns[1]+ns[2])/(ns[1]*ns[2])))*sqrt(SST-SSW)/sqrt(s1/ns[1] + s2/ns[2])
  if(is.na(t.stat)){ ## SST - SSW must be negative
    t.stat = 0
  }
  t.stat
}

# WT.test = function(dm, f, nrep=999){
#   stats = c(WT(dm, f), replicate(nrep, WT(dm, f[sample(length(f))])))
#   p.value = sum(stats>=stats[1])/(nrep+1)
#   t.stat = stats[1]
#   list(p.value = p.value, t.stat = t.stat, nrep=nrep)  
# }

strata.permute = function(f, strata){
  p = unlist(tapply(f,strata,sample)) 
  p2 = unlist(tapply(1:length(f), strata, identity))
  f[p[p2]]
}

WT.test = function(dm, f, strata, nrep=999){
  stats = c(WT(dm, f), 
            replicate(nrep,
                      WT(dm, strata.permute(f, strata))
                      )
            )
  p.value = sum(stats>=stats[1])/(nrep+1)
  t.stat = stats[1]
  list(p.value = p.value, t.stat = t.stat, nrep=nrep)  
}

#Wstar.test=function(dm, f, nrep=999){ 
  ## This method is not currently implemented. 
  ## The function is provided as a placeholder.
#  list(p.value = 1, Wstar = "Unimplemented", nrep=nrep)
#}