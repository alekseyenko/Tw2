get.community = function(rd, otus=1000){
  xx=rep(0,otus)
  xx[sample(otus, round(otus*rd), replace=F)]=1
  xx
}


simulate.sample = function(group.size, rd1=.25, rd2=.25, effect=0.0){
  a=cbind(replicate(group.size, get.community(rd1)),
          replicate(group.size, get.community(rd2)))
  notus = dim(a)[1]
  if(round(notus*effect)>0){
    unique = sample(notus, round(notus*effect), replace=F)
    zero = matrix(0,nrow=length(unique), ncol=group.size)
    a= cbind(rbind(a[unique,1:group.size],
                   a[-unique, 1:group.size],
                   zero),
             rbind(zero,
                   a[unique, -(1:group.size)],
                   a[-unique, -(1:group.size)]))
  }
  a = a[rowSums(a)>0,] 
}


simulate.normal = function(group.size, fsd=1, effect=0.0){
  a=cbind(matrix(rnorm(group.size*1000, mean=0+effect/sqrt(1000), sd = 1), nrow=1000, ncol=group.size),
          matrix(rnorm(group.size*1000, mean=0, sd = 1*fsd), nrow=1000, ncol=group.size))
  a
}

simulate.normal2 = function(group.size1, group.size2, fsd=1, effect=0.0){
  a=cbind(matrix(rnorm(group.size1*1000, mean=0, sd = 1), nrow=1000, ncol=group.size1),
          matrix(rnorm(group.size2*1000, mean=0+effect/sqrt(1000), sd = 1*fsd), nrow=1000, ncol=group.size2))
  a
}