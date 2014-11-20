discret.equalwidth<-function(dataloc, n.class) {
# dataloc is a vector

  dataloc.discret=rep(NA, length(dataloc))

  dataloc.min=min(dataloc)
  dataloc.max=max(dataloc)
  dataloc.fra1=(dataloc.max-dataloc.min)/n.class

  dataloc.frac.lo=dataloc.min
  for (i in 1:n.class) {
  # i=1
    dataloc.frac.hi=dataloc.min+dataloc.fra1*i
    idx.frac=intersect(which(dataloc<dataloc.frac.hi), which(dataloc>=dataloc.frac.lo))
    dataloc.discret[idx.frac]=i
    dataloc.frac.lo=dataloc.frac.hi
  }

  dataloc.discret[which(dataloc==max(dataloc))]=n.class

  return(dataloc.discret)
}

entropy<-function(C) {
# calculate class entropy for the data S
# which has class label C
 class.data=unique(C); n.data=length(C)
 entropy.data=0
 for (class1.data in class.data) {
 # class1.data=class.data[1]
   if (any(which(C==class1.data))) {
     n.class1.data=length(which(C==class1.data))
     pi=n.class1.data/n.data
     entropy.data=entropy.data-(pi*log2(pi))
   }
 }
 return(entropy.data)
}

infogain<-function(s,a) {
# this function calculates 
# infogain(s|a)=H(s)-H(s|a)
# H is entropy
 class.a=unique(a)
 gain.sa=entropy(s)
 for (class1.a in class.a) {
   if (any(which(a==class1.a))) {
     idx.class1.a=which(a==class1.a)
     n.class1.a=length(idx.class1.a)

     s.class1.a=s[idx.class1.a]    
     gain.sa=gain.sa-n.class1.a/length(s)*entropy(s.class1.a)
   }
 }
 return(gain.sa)
}

infogain.normalized<-function(x,y) {
# symmetrical uncertainty
  su.xy=2*(infogain(x,y)/(entropy(x)+entropy(y)))
  return(su.xy)
}
