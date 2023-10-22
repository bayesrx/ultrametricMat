library(ape)
library(phylobase)

#### The same tree from rtree with the same nSmp ####
setwd("~/Research/UltrametricsMat/SimData/")
# set.seed(1234)
set.seed(123456)
nLeaf<-10
trueTr<-phylo4(rtree(n = nLeaf))
sinEdge<-runif(1)
plot(trueTr)
trueTrCov<-matrix(0,nr=nLeaf,nc=nLeaf)
for(j in 1:nLeaf){
  for(k in j:nLeaf){
    trueTrCov[j,k]<-trueTrCov[k,j]<-
      phylobase::nodeHeight(trueTr, node = MRCA(trueTr,paste0("t",c(j,k))),from="root")
  }
}
trueTrCov<-trueTrCov+sinEdge


nSmp.vec<-nLeaf*c(3,5,10,25,50)
nRep<-50
res.ls<-list()
for(nIdx in 1:length(nSmp.vec)){
  nSmp<-nSmp.vec[nIdx]
  df.ls<-list()
  for(rIdx in 1:nRep){
    n.ObsDf<- mvtnorm::rmvnorm(nSmp,sigma=trueTrCov)
    t4.ObsDf<-mvtnorm::rmvt(nSmp,sigma=trueTrCov, df=4)
    t3.ObsDf<-mvtnorm::rmvt(nSmp,sigma=trueTrCov, df=3)
    # cov(ObsDf)
    tmp.ls<-list(n.ObsDf=n.ObsDf, t3.ObsDf=t3.ObsDf, t4.ObsDf=t4.ObsDf)
    df.ls[[rIdx]]<-tmp.ls
  }
  res.ls[[nIdx]]<-df.ls
  names(res.ls)[nIdx]<-paste0("nSmp_",nSmp)
}

# saveRDS(list(df=res.ls, trueTr=trueTr, trueTrCov=trueTrCov),paste0("./rtree_nLeaf",nLeaf,".RDS"))
saveRDS(list(df=res.ls, trueTr=trueTr, trueTrCov=trueTrCov,rtEdge=sinEdge),paste0("./rtEdge_rtree_nLeaf",nLeaf,".RDS"))
