library(partitions)
library(stringr)
library(ape)
library(phylobase)
library(mvtnorm)
library(truncnorm)

list_allSplits<-function(nLeaf){
  n<-nLeaf+1
  intPar<-parts(n)
  intSplt<-intPar[,colSums(intPar>0)==2]
  intSplt<-matrix(intSplt[,colSums(intSplt>1)==2],nr=n)
  res<-data.frame(matrix(0,nr=2^nLeaf-nLeaf-2,ncol=3))
  colnames(res)<-c("eNames","par1","par2")
  res$eNames<-paste0("e",1:nrow(res))
  idx<-0
  for(cIdx in 1:ncol(intSplt)){
    cur.intSplt<-intSplt[1:2,cIdx]
    cur.leafSplt<-setparts(cur.intSplt)
    for(sIdx in 1:ncol(cur.leafSplt)){
      cur.splt<-cur.leafSplt[,sIdx]
      idx<-idx+1
      splt1<-which(cur.splt==1)-1
      if(any(splt1==0)){
        splt2<-which(cur.splt==2)-1
      }else{
        splt2<-splt1
        splt1<-which(cur.splt==2)-1
      }
      par1<-paste(splt1, collapse=",")
      par2<-paste(splt2, collapse=",")
      res[idx,2:3]<-c(par1,par2)
    }
  }
  return(res)
}

chk_compatible<-function(splt1,splt2){
  # input splt: a vector of two sub-blocks. Each block is represented by a string. E.g.: splt1: c("0,1,2,5,6,7","3,4")
  e1.par1<-as.numeric(str_split(splt1[1],",",simplify = T))
  e1.par2<-as.numeric(str_split(splt1[2],",",simplify = T))
  e2.par1<-as.numeric(str_split(splt2[1],",",simplify = T))
  e2.par2<-as.numeric(str_split(splt2[2],",",simplify = T))
  return(any(length(intersect(e1.par1,e2.par1))==0,
             length(intersect(e1.par1,e2.par2))==0,
             length(intersect(e1.par2,e2.par1))==0,
             length(intersect(e1.par2,e2.par2))==0))
}

allSplt_compatible<-function(inp.dfSplt){
  resMat<-matrix(NA,ncol=nrow(inp.dfSplt),nrow=nrow(inp.dfSplt))
  diag(resMat)<-0
  colnames(resMat)<-rownames(resMat)<-inp.dfSplt$eNames
  for(sIdx in 1:(nrow(inp.dfSplt)-1)){
    splt1<-c(inp.dfSplt[sIdx,2],inp.dfSplt[sIdx,3])
    for(nIdx in (sIdx+1):nrow(inp.dfSplt)){
      splt2<-c(inp.dfSplt[nIdx,2],inp.dfSplt[nIdx,3])
      resMat[sIdx,nIdx]<-as.numeric(chk_compatible(splt1,splt2))
    }
  }
  resMat[lower.tri(resMat)]<-t(resMat)[lower.tri(resMat)]
  return(resMat)
}

gen_rndTrFromCompMat<-function(inp.dfSplt, inp.compMat, nLeaf_, nTr, seed_=1){
  set.seed(seed_)
  tipTr<-spltTr<-list()
  for(tIdx in 1:nTr){
    spltMat<-data.frame(matrix(NA,nc=3,nr=nLeaf-2))
    colnames(spltMat)<-c("splitIdx","par","length")
    spltMat$length<-rexp(nrow(spltMat))
    e1<-sample(1:nrow(inp.compMat),1)
    spltMat[1,1]<-e1
    spltMat[1,2]<-inp.dfSplt[e1,"par2"]
    e2<-sample(which(inp.compMat[e1,]==1),1)
    spltMat[2,1]<-e2
    spltMat[2,2]<-inp.dfSplt[e2,"par2"]
    for(eIdx in 3:nrow(spltMat)){
      spltIdx<-spltMat[1:(eIdx-1),1]
      spltMat[eIdx,1]<-sample(which(colSums(inp.compMat[spltIdx,])==length(spltIdx)),1)
      spltMat[eIdx,2]<-inp.dfSplt[spltMat[eIdx,1],"par2"]
    }
    spltMat<-spltMat[,2:3]
    spltTr[[tIdx]]<-spltMat[order(spltMat[,"par"]),]
    tipMat<-data.frame(matrix(NA,nc=2,nr=nLeaf))
    colnames(tipMat)<-c("par","length")
    tipMat[,"par"]<-1:nLeaf
    tipMat[,"length"]<-rexp(nLeaf)
    tipTr[[tIdx]]<-tipMat
  }
  return(list(split=spltTr,tip=tipTr))
}

lomega_<-function(n,beta_){
  res<-lgamma(n+1+beta_)
  return(res)
}

betaSplt.tbl.gen<-function(n,beta_){
  res<-list()
  for(idx in 2:n){
    tmp<-as.matrix(parts(idx))
    intPart.idx<-apply(tmp,2,function(x) sum(x!=0)==2)
    intPart<-as.matrix(tmp[1:2,intPart.idx])
    lProb<-rep(NA,sum(intPart.idx))
    for(cIdx in 1:ncol(intPart)){
      if(intPart[1,cIdx]==intPart[2,cIdx]){
        k=intPart[1,cIdx]
        m=idx
        c<-choose(m-1,k-1)
      }else{
        k1=intPart[1,cIdx]
        m=idx
        c1<-choose(m-1,k1-1)
        k2=intPart[2,cIdx]
        c2<-choose(m-1,k2-1)
        c<-c1+c2
      }
      lProb[cIdx]<-lomega_(intPart[1,cIdx],beta_) + lomega_(intPart[2,cIdx],beta_) + log(c)
    }
    prob<-exp(lProb - matrixStats::logSumExp(lProb))
    intPart<-data.frame(cbind(t(intPart),prob))
    res[[idx-1]]<-intPart
  }
  names(res)<-paste0("n=",2:n)
  return(res)
}

betaSplt.prior<-function(inp.trTbl, inp.betaSplt.tbl, nLeaf){
  ancetr.Edges<-sort(unique(inp.trTbl$ancestor))
  res.prior<-0
  for(eIdx in 1:length(ancetr.Edges)){
    cur.anstr.edge<-ancetr.Edges[eIdx]
    desc.edge<-inp.trTbl[inp.trTbl$ancestor==cur.anstr.edge,]
    cur.Splt.ls<-lapply(desc.edge$par,function(x) as.numeric(str_split_1(x,",")))
    par1<-cur.Splt.ls[[1]]
    par2<-cur.Splt.ls[[2]]
    n.par1<-length(par1)
    n.par2<-length(par2)
    n.par1.ord<-max(n.par1,n.par2)
    n.par2.ord<-min(n.par1,n.par2)
    cur.betaSplt.tbl<-inp.betaSplt.tbl[[paste0("n=",n.par1+n.par2)]]
    prior.prob<-cur.betaSplt.tbl
    res.prior<-res.prior+log(cur.betaSplt.tbl$prob[cur.betaSplt.tbl$V1==n.par1.ord & cur.betaSplt.tbl$V2==n.par2.ord])
  }
  return(res.prior)
}

find_cmpSplt_noCompTbl<-function(inp.intMdt, inp.shrnkg.splt, nLeaf){
  intMdt.Tr<-inp.intMdt
  intMdt.Tr$parC<-unlist(apply(intMdt.Tr,1,function(x) paste(sort(setdiff(0:nLeaf,as.numeric(str_split_1(x[1],",")))),collapse = ",")))
  shrnkg.splt<-inp.shrnkg.splt
  shrnkg.splt.par<-as.numeric(str_split_1(shrnkg.splt$par,","))
  anstr.idx<-unlist(lapply(intMdt.Tr$par,function(x) all(shrnkg.splt.par %in% str_split_1(x,","))))
  anstr<-intMdt.Tr[anstr.idx,]
  
  # when parent is the root
  if(sum(anstr.idx)==0){
    parent.par<-paste(1:nLeaf,collapse = ",")
    siblingID<-intMdt.Tr[intMdt.Tr$parC==paste(sort(c(0,shrnkg.splt.par)),collapse = ","),]
  }else{
    parent<-anstr[which.min(nchar(anstr$par)),]
    parent.par<-parent$par
    siblingID<-intMdt.Tr[intMdt.Tr$par==paste(sort(setdiff(
      as.numeric(str_split_1(parent$par,",")),shrnkg.splt.par)),collapse = ","),]
  }
  
  # when sibling is a tip
  if(nrow(siblingID)==0){
    parent.parVec<-as.numeric(str_split_1(parent.par,","))
    sibling.par<-setdiff(parent.parVec,shrnkg.splt.par)
  }else{
    sibling.par<-as.numeric(str_split_1(siblingID$par,","))
  }
  
  # when the descendant of the shrinkage split are tips
  if(length(shrnkg.splt.par)==2){
    cmpt.Splt1<-paste(sort(c(sibling.par,shrnkg.splt.par[1])),collapse = ",")
    cmpt.Splt2<-paste(sort(c(sibling.par,shrnkg.splt.par[2])),collapse = ",")
  }else{
    desc.idx<-unlist(lapply(intMdt.Tr$par,function(x) length(setdiff(as.numeric(str_split_1(x,",")),shrnkg.splt.par))==0))
    desc<-intMdt.Tr[desc.idx,]
    desc1<-desc[which.max(nchar(desc$par)),]
    desc1.par<-as.numeric(str_split_1(desc1$par,","))
    desc2.par<-sort(setdiff(shrnkg.splt.par,desc1.par))
    cmpt.Splt1<-paste(sort(c(sibling.par,desc1.par)),collapse = ",")
    cmpt.Splt2<-paste(sort(c(sibling.par,desc2.par)),collapse = ",")
  }
  return(c(cmpt.Splt1, cmpt.Splt2))
}

NNIpropTr<-function(inp.curSplt, inp.compMat=NULL, nLeaf, ctsSplt.sd){
  # set.seed(seed_)
  n.shrnkgEdge=1
  # inp.dfSplt[inp.curSpltID[,1],]
  shrkgEdge.idx<-sample(1:nrow(inp.curSplt),n.shrnkgEdge)
  intMdt.Tr<-inp.curSplt[-shrkgEdge.idx,]
  shrkg.Splt<-inp.curSplt[shrkgEdge.idx,]
  newEdge<-data.frame(matrix(NA,nrow=n.shrnkgEdge,ncol=2))
  colnames(newEdge)<-colnames(shrkg.Splt)
  l.tplg.prob<-0
  
  if(is.null(inp.compMat)){
    # print("NNI")
    #### find compatible split through NNI alg
    # candt.noMat<-find_cmpSplt_noCompTbl(inp.intMdt.ID = intMdt.Tr, inp.shrnkg.spltID = shrkg.SpltID, inp.dfSplt = inp.dfSplt)
    candt<-find_cmpSplt_noCompTbl(inp.intMdt = intMdt.Tr, inp.shrnkg.splt = shrkg.Splt, nLeaf = nLeaf)
  }else{
    # print("comp table")
    spltIdx<-inp.compMat$spltIdx
    compMat<-inp.compMat$comp.mat
    intMdt.Tr.idx<-spltIdx$par2 %in% intMdt.Tr$par
    #### find compatible split through full table
    tmpMat<-matrix(compMat[intMdt.Tr.idx,],ncol=ncol(compMat))
    candtID<-which(colSums(tmpMat)==nrow(intMdt.Tr))
    # candt.mat<-candt[candt!=shrkg.SpltID[,1]]
    candt<-spltIdx$par2[candtID]
    candt<-candt[candt!=shrkg.Splt$par]
  }
  candt<-sort(candt)
  
  # candt.noMat<-find_cmpSplt_noCompTbl(inp.intMdt.ID = intMdt.Tr, inp.shrnkg.spltID = shrkg.SpltID, inp.dfSplt = inp.dfSplt)
  # tmpMat<-matrix(inp.compMat[intMdt.Tr[,1],],ncol=ncol(inp.compMat))
  # candt<-which(colSums(tmpMat)==nrow(intMdt.Tr))
  # candt.mat<-candt[candt!=shrkg.SpltID[,1]]
  # if(all(candt.noMat%in% candt.mat)){
  #   candt=candt.mat
  # }else{
  #   stop("different candidate compatible splits")
  # }
  if(ctsSplt.sd==0){
    # print("discrete")
    newEdge[1,1]<-sample(candt,1)
    l.tplg.prob<-l.tplg.prob+log(length(candt))
    newEdge[1,2]<-inp.curSplt[shrkgEdge.idx,2]
    intMdt.Tr<-rbind(intMdt.Tr[,1:2],newEdge)
    intMdt.Tr<-intMdt.Tr[order(intMdt.Tr[,1]),]
    resTr<-intMdt.Tr
    colnames(resTr)<-c("par","length")
  }else{
    # print("partial cts")
    # tmpSpltBr<-rnorm(1,mean=0,sd=ctsSplt.sd)
    # if(tmpSpltBr>0){
    #   newEdge[1,1]<-candt[1]
    # }else{
    #   newEdge[1,1]<-candt[2]
    # }
    
    SpltBrDiff<-rnorm(1,mean=0,sd=ctsSplt.sd)
    newSpltBrLen<-inp.curSplt[shrkgEdge.idx,"length"]+SpltBrDiff
    if(newSpltBrLen>0){
      newEdge[1,1]<-inp.curSplt[shrkgEdge.idx,"par"]
      newEdge[1,2]<-newSpltBrLen
      l.tplg.prob<-dnorm(x=SpltBrDiff,mean=0,sd=ctsSplt.sd, log=T)
    }else{
      newSplt<-sample(candt,1)
      newEdge[1,1]<-newSplt
      newEdge[1,2]<-abs(newSpltBrLen)
      l.tplg.prob<-dnorm(x=SpltBrDiff,mean=0,sd=ctsSplt.sd, log=T)+log(1/2)
    }
    intMdt.Tr<-rbind(intMdt.Tr[,1:2],newEdge)
    intMdt.Tr<-intMdt.Tr[order(intMdt.Tr[,1]),]
    resTr<-intMdt.Tr
    colnames(resTr)<-c("par","length")
  }
  resTr<-resTr[order(resTr$par),]
  return(list(propTr=resTr,l.jpProb=l.tplg.prob))
}

NNIjumpPr<-function(inp.tgtSplt, inp.orgSplt, ctsSplt.sd){
  comnSpltIdx<-inp.tgtSplt[,1] %in% inp.orgSplt[,1]
  if(all(comnSpltIdx)){
    newLenIdx<-!(inp.tgtSplt[,2] %in% inp.orgSplt[,2])
    newLen<-inp.tgtSplt[newLenIdx,2]
    orgLenIdx<-!(inp.orgSplt[,2] %in% inp.tgtSplt[,2])
    orgLen<-inp.orgSplt[orgLenIdx,2]
    q.lprob<-dnorm(x=abs(newLen-orgLen),mean=0,sd=ctsSplt.sd,log=T)
  }else{
    newSplt<-inp.tgtSplt[!comnSpltIdx,]
    replacedSpltIdx<-!(inp.orgSplt[,1] %in% inp.tgtSplt[,1])
    replacedSplt<-inp.orgSplt[replacedSpltIdx,]
    q.lprob<-dnorm(x=abs(newSplt[,2])+abs(replacedSplt[,2]),mean=0,sd=ctsSplt.sd,log=T) + log(1/2)
  }
  return(q.lprob)
}

update_SpltID<-function(inp.curSplt, inp.tip, inp.rtEdge,
                        inp.compMat=NULL, inp.Obsdf, inp.betaSplt.tbl, inp.ctsSplt.sd, 
                        hypParam_mean,  nLeaf){
  oldTrCov<-splt2TrCov(inp.splt = inp.curSplt, inp.pedEdge = inp.tip,inp.rtEdge = inp.rtEdge)
  # propSpltID.ls<-propTr(inp.curSpltID = inp.curSpltID, inp.compMat = inp.compMat,n.shrnkgEdge = 1)
  propSplt.ls<-NNIpropTr(inp.curSplt = inp.curSplt, nLeaf=nLeaf, ctsSplt.sd = inp.ctsSplt.sd, inp.compMat = inp.compMat)
  propSplt<-propSplt.ls$propTr
  propTrCov<-splt2TrCov(inp.splt = propSplt, inp.pedEdge = inp.tip,inp.rtEdge = inp.rtEdge)
  
  # prior
  old.TrTbl<-splt2TrTbl(inp.splt = inp.curSplt, inp.pedEdge = inp.tip)
  prop.TrTbl<-splt2TrTbl(inp.splt = propSplt, inp.pedEdge = inp.tip)
  old.prior<-betaSplt.prior(inp.trTbl = old.TrTbl, inp.betaSplt.tbl = inp.betaSplt.tbl, nLeaf)
  prop.prior<-betaSplt.prior(inp.trTbl = prop.TrTbl, inp.betaSplt.tbl = inp.betaSplt.tbl, nLeaf)
  prior.diff<-prop.prior - old.prior
  if(inp.ctsSplt.sd>0){
    propLenIdx<- !(propSplt[,2] %in% inp.curSplt[,2])
    propLen<-propSplt[propLenIdx,2]
    oldLenIdx<-!(inp.curSplt[,2] %in% propSplt[,2])
    oldLen<-inp.curSplt[oldLenIdx,2]
    oldLen.prior<-dexp(oldLen,hypParam_mean)
    propLen.prior<-dexp(propLen,hypParam_mean)
    prior.diff<-prior.diff + (propLen.prior-oldLen.prior)
  }
  
  #llh
  old.llh<-sum(dmvnorm(inp.Obsdf,sigma = oldTrCov,log = T))
  prop.llh<-sum(dmvnorm(inp.Obsdf,sigma = propTrCov,log = T))
  llh.diff<-prop.llh-old.llh
  
  #jump Prob
  if(inp.ctsSplt.sd==0){
    q.diff<-0
  }else{
    old2prop<-propSplt.ls$l.jpProb
    prop2old<-NNIjumpPr(inp.tgtSplt = inp.curSplt, inp.orgSplt = propSplt, ctsSplt.sd = inp.ctsSplt.sd)
    q.diff<-prop2old-old2prop
    # print(q.diff)
  }
  
  diff<-prior.diff + llh.diff + q.diff
  u<-runif(1)
  if(diff>log(u)){
    ac<-1
    return(list(split= propSplt, llh=prop.llh, acpt=ac))
  }else{
    ac<-0
    return(list(split= inp.curSplt, llh=old.llh, acpt=ac))
  }
}

update_IntLenTr<-function(inp.curSplt, inp.tip, inp.rtEdge,
                          inp.Obsdf, hypParam_mean ,inp.intBr.sd, nLeaf){
  resTrCov<-splt2TrCov(inp.splt = inp.curSplt, inp.pedEdge = inp.tip, inp.rtEdge = inp.rtEdge)
  nIntBr<-nrow(inp.curSplt)
  resIntBr<-inp.curSplt
  # acpt.vec<-rep(0,nIntBr)
  for(sIdx in 1:nIntBr){
    intTrCov<-resTrCov

    cur.splt<-inp.curSplt[sIdx,]
    cur.par<-as.numeric(str_split_1(cur.splt[1,1],","))
    tmp.E<-matrix(0,nrow=nLeaf,ncol=nLeaf)
    tmp.E[cur.par,cur.par]<-1
    intTrCov <- intTrCov - tmp.E*cur.splt[1,2]
    
    # propIntBr<-rexp(1)
    propIntBr<-rtruncnorm(1, a=0, b=Inf, mean = cur.splt[1,2], sd=inp.intBr.sd)
    propTrCov<-intTrCov 
    propTrCov <- propTrCov + tmp.E*propIntBr
    
    # prior
    old.prior<-dexp(cur.splt[1,2],hypParam_mean)
    prop.prior<-dexp(propIntBr,hypParam_mean)
    prior.diff<-prop.prior - old.prior
    
    #llh
    res.llh<-sum(dmvnorm(inp.Obsdf,sigma = resTrCov,log = T))
    prop.llh<-sum(dmvnorm(inp.Obsdf,sigma = propTrCov,log = T))
    llh.diff<- prop.llh - res.llh
    
    #jump Prob
    old2prop<-dtruncnorm(propIntBr, a=0, b=Inf, mean = cur.splt[1,2], sd=inp.intBr.sd)
    prop2old<-dtruncnorm(cur.splt[1,2], a=0, b=Inf, mean = propIntBr, sd=inp.intBr.sd)
    q.diff<-prop2old-old2prop
    
    diff<-prior.diff + llh.diff + q.diff
    u<-runif(1)
    if(diff>log(u)){
      resTrCov<-propTrCov
      resIntBr[sIdx,2]<-propIntBr
    }
  }
  resIntBr<-resIntBr[order(resIntBr[,1]),]
  return(list(TrCov=resTrCov, intBr=resIntBr))
}

update_TipLenTr<-function(inp.curSplt, inp.tip, inp.Obsdf, inp.rtEdge,
                          hypParam_mean, inp.tipBr.sd, nLeaf){
  resTrCov<-splt2TrCov(inp.splt = inp.curSplt, inp.pedEdge = inp.tip, inp.rtEdge = inp.rtEdge)
  resTip<-inp.tip
  acpt.vec<-rep(0,nLeaf)
  
  for(tIdx in 1:nLeaf){
    intTrCov<-resTrCov
    intTrCov[tIdx,tIdx]<-intTrCov[tIdx,tIdx]-inp.tip[tIdx,2]
    
    # propTip<-rexp(1)
    propTip<-rtruncnorm(1, a=0, b=Inf, mean = inp.tip[tIdx,2], sd=inp.tipBr.sd)
    propTrCov<-intTrCov 
    propTrCov[tIdx,tIdx]<-propTrCov[tIdx,tIdx] + propTip
    
    # prior
    old.prior<-dexp(inp.tip[tIdx,2],hypParam_mean)
    prop.prior<-dexp(propTip,hypParam_mean)
    prior.diff<-prop.prior - old.prior
    
    res.llh<-sum(dmvnorm(inp.Obsdf,sigma = resTrCov,log = T))
    prop.llh<-sum(dmvnorm(inp.Obsdf,sigma = propTrCov,log = T))
    llh.diff<- prop.llh - res.llh
    
    #jump Prob
    old2prop<-dtruncnorm(propTip, a=0, b=Inf, mean = inp.tip[tIdx,2], sd=inp.tipBr.sd)
    prop2old<-dtruncnorm(inp.tip[tIdx,2], a=0, b=Inf, mean = propTip, sd=inp.tipBr.sd)
    q.diff<-prop2old-old2prop
    
    diff<-prior.diff + llh.diff + q.diff
    u<-runif(1)
    
    if(diff>log(u)){
      acpt.vec[tIdx]<-1
      resTrCov<-propTrCov
      resTip[tIdx,2]<-propTip
      res.llh<-prop.llh
    }
  }
  resTip<-resTip[order(as.numeric(resTip[,1])),]
  return(list(TrCov=resTrCov, tipBr=resTip ,acpt=acpt.vec, llh=res.llh))
}

update_rtEdgeLen<-function(inp.curSplt, inp.tip, inp.Obsdf, inp.rtEdge,
                           hypParam_mean, inp.tipBr.sd, nLeaf){
  resTrCov<-splt2TrCov(inp.splt = inp.curSplt, inp.pedEdge = inp.tip, inp.rtEdge = inp.rtEdge)
  resRtEdge<-inp.rtEdge
  acpt<-0
  
  intTrCov<-resTrCov
  intTrCov<-intTrCov-inp.rtEdge
  
  # propTip<-rexp(1)
  propRtEdge<-rtruncnorm(1, a=0, b=Inf, mean = inp.rtEdge, sd=inp.tipBr.sd)
  propTrCov<-intTrCov 
  propTrCov<-propTrCov + propRtEdge
  
  # prior
  old.prior<-dexp(inp.rtEdge,hypParam_mean)
  prop.prior<-dexp(propRtEdge,hypParam_mean)
  prior.diff<-prop.prior - old.prior
  
  res.llh<-sum(dmvnorm(inp.Obsdf,sigma = resTrCov,log = T))
  prop.llh<-sum(dmvnorm(inp.Obsdf,sigma = propTrCov,log = T))
  llh.diff<- prop.llh - res.llh
  
  #jump Prob
  old2prop<-dtruncnorm(propRtEdge, a=0, b=Inf, mean = inp.rtEdge, sd=inp.tipBr.sd)
  prop2old<-dtruncnorm(inp.rtEdge, a=0, b=Inf, mean = propRtEdge, sd=inp.tipBr.sd)
  q.diff<-prop2old-old2prop
  
  diff<-prior.diff + llh.diff + q.diff
  u<-runif(1)
  
  if(diff>log(u)){
    acpt<-1
    resTrCov<-propTrCov
    resRtEdge<-propRtEdge
    res.llh<-prop.llh
  }
  return(list(TrCov=resTrCov, rtEdgeBr=resRtEdge ,acpt=acpt, llh=res.llh))
}

splt2TrCov<-function(inp.splt,inp.pedEdge,inp.rtEdge=0){
  nLeaf<-nrow(inp.pedEdge)
  res.TrCov<-matrix(0,nrow=nLeaf,ncol=nLeaf)
  for(sIdx in 1:nrow(inp.splt)){
    cur.par<-inp.splt[sIdx,1]
    cur.par<-as.numeric(str_split_1(cur.par,","))
    cur.length<-inp.splt[sIdx,2]
    tmp.E<-matrix(0,nrow=nLeaf,ncol=nLeaf)
    tmp.E[cur.par,cur.par]<-1
    res.TrCov<-res.TrCov + tmp.E * cur.length
  }
  for(lIdx in 1:nrow(inp.pedEdge)){
    leaf<-as.numeric(inp.pedEdge[lIdx,1])
    res.TrCov[leaf,leaf]<-res.TrCov[leaf,leaf]+as.numeric(inp.pedEdge[lIdx,2])
  }
  res.TrCov<-res.TrCov+inp.rtEdge
  return(res.TrCov)
}

splt2TrTbl<-function(inp.splt, inp.pedEdge){ # TrTbl: Table from Phylo4 print output
  nLeaf<-nrow(inp.pedEdge)
  resTbl<-data.frame(inp.splt)
  tipTbl<-data.frame(inp.pedEdge)
  colnames(tipTbl)<-c("par","length")
  resTbl<-rbind(tipTbl,resTbl)
  resTbl$ancestor<-NA
  for(eIdx in 1:nrow(resTbl)){
    cur.edge<-resTbl[eIdx,]
    cur.leaf.vec<-as.numeric(str_split_1(cur.edge$par,","))
    rest.edge<-resTbl[-eIdx,]
    anstr.idx<-unlist(lapply(rest.edge$par, function(x) all(cur.leaf.vec %in% str_split_1(x,","))))
    anstr.edge<-rest.edge[anstr.idx,]
    parent.idx<-which.min(unlist(lapply(anstr.edge$par,function(x) length(str_split_1(x,",")))))
    cur.parent<-anstr.edge[parent.idx,]
    if(length(parent.idx)==0){
      resTbl$ancestor[eIdx]<-0
    }else{
      resTbl$ancestor[eIdx]<-cur.parent[,"par"]
    }
  }
  return(resTbl)
}

trTbl2Tr<-function(inp.trTbl){
  wrkTbl<-inp.trTbl
  wrkTbl$edge<-NA
  wrkTbl$edge[str_detect(wrkTbl$par,",")]<-"i"
  wrkTbl$edge[!str_detect(wrkTbl$par,",")]<-"t"
  wrkTbl$edge[wrkTbl$edge=="i"]<-paste0("i",1:sum(wrkTbl$edge=="i"))
  wrkTbl$edge[wrkTbl$edge=="t"]<-paste0("t",as.numeric(wrkTbl$par[wrkTbl$edge=="t"]))
  anstr.idx<-wrkTbl$ancestor=="0" 
  descd.edge<-wrkTbl[anstr.idx,]
  final.txt<-paste0("(",paste(paste(descd.edge$edge,descd.edge$length,sep=":"),collapse = ","),"):0;")
  wrkTbl<-wrkTbl[!anstr.idx,]
  cur.anstr.set<-descd.edge
  cur.anstr.set<-cur.anstr.set[!str_detect(cur.anstr.set$edge,"t"),]
  while(nrow(wrkTbl)!=0){
    nxt.anstr.set<-data.frame(matrix(NA,nr=0,nc=4))
    colnames(nxt.anstr.set)<-colnames(wrkTbl)
    for(aIdx in 1:nrow(cur.anstr.set)){
      cur.anstr<-cur.anstr.set[aIdx,]
      anstr.idx<-wrkTbl$ancestor==cur.anstr$par
      descd.edge<-wrkTbl[anstr.idx,]
      subTr.txt<-paste0("(",paste(paste(descd.edge$edge,descd.edge$length,sep=":"),collapse = ","),")")
      final.txt<-str_replace(final.txt,cur.anstr$edge,subTr.txt)
      wrkTbl<-wrkTbl[!anstr.idx,]
      nxt.anstr.set<-rbind(nxt.anstr.set,descd.edge)
      nxt.anstr.set<-nxt.anstr.set[!str_detect(nxt.anstr.set$edge,"t"),]
    }
    cur.anstr.set<-nxt.anstr.set
  }
  tr<-read.tree(text=final.txt)
  return(list(phylo4=phylo4(tr),txt=final.txt))
}

# trTbl2Tr<-function(inp.trTbl){
#   wrkTbl<-inp.trTbl
#   wrkTbl$edge<-NA
#   wrkTbl$edge[str_detect(wrkTbl$par,",")]<-"i"
#   wrkTbl$edge[!str_detect(wrkTbl$par,",")]<-"t"
#   wrkTbl$edge[wrkTbl$edge=="i"]<-paste0("i",1:sum(wrkTbl$edge=="i"))
#   wrkTbl$edge[wrkTbl$edge=="t"]<-paste0("t",as.numeric(wrkTbl$par[wrkTbl$edge=="t"]))
#   anstr.idx<-wrkTbl$ancestor=="0" 
#   descd.edge<-wrkTbl[anstr.idx,]
#   final.txt<-paste0("(",descd.edge[1,"edge"],":",descd.edge[1,"length"],",",
#                     descd.edge[2,"edge"],":",descd.edge[2,"length"],");")
#   wrkTbl<-wrkTbl[!anstr.idx,]
#   cur.anstr.set<-descd.edge
#   cur.anstr.set<-cur.anstr.set[!str_detect(cur.anstr.set$edge,"t"),]
#   while(nrow(wrkTbl)!=0){
#     nxt.anstr.set<-data.frame(matrix(NA,nr=0,nc=4))
#     colnames(nxt.anstr.set)<-colnames(wrkTbl)
#     for(aIdx in 1:nrow(cur.anstr.set)){
#       cur.anstr<-cur.anstr.set[aIdx,]
#       anstr.idx<-wrkTbl$ancestor==cur.anstr$par
#       descd.edge<-wrkTbl[anstr.idx,]
#       subTr.txt<-paste0("(",descd.edge[1,"edge"],":",descd.edge[1,"length"],",",
#                         descd.edge[2,"edge"],":",descd.edge[2,"length"],")")
#       final.txt<-str_replace(final.txt,cur.anstr$edge,subTr.txt)
#       wrkTbl<-wrkTbl[!anstr.idx,]
#       nxt.anstr.set<-rbind(nxt.anstr.set,descd.edge)
#       nxt.anstr.set<-nxt.anstr.set[!str_detect(nxt.anstr.set$edge,"t"),]
#     }
#     cur.anstr.set<-nxt.anstr.set
#   }
#   tr<-read.tree(text=final.txt)
#   return(list(phylo4=phylo4(tr),txt=final.txt))
# }

splt2Phylo4<-function(inp.splt,inp.pedEdge){
  tmp.TrTbl<-splt2TrTbl(inp.splt = inp.splt, inp.pedEdge = inp.pedEdge)
  res.ls<-trTbl2Tr(inp.trTbl = tmp.TrTbl)
  return(res.ls$phylo4)
}

phylo42Splt<-function(phylo4_){
  nLeaf<-nTips(phylo4_)
  spltMat<-data.frame(matrix(NA,nc=2,nr=nLeaf-2))
  tipMat<-data.frame(matrix(NA,nc=2,nr=nLeaf))
  colnames(tipMat)<-colnames(spltMat)<-c("par","length")
  for(sIdx in 1:nrow(spltMat)){
    intNode<-nLeaf+1+sIdx
    spltMat[sIdx,1]<-paste(sort(as.numeric(str_sub(names(
      descendants(phy=phylo4_,node=intNode,type="tips")),2))),collapse = ",")
    spltMat[sIdx,2]<-edgeLength(phylo4_)[getEdge(phylo4_,intNode,"descendant")]
  }
  for(tIdx in 1:nLeaf){
    tipNode<-paste0("t",tIdx)
    tipMat[tIdx,1]<-as.character(tIdx)
    tipMat[tIdx,2]<-edgeLength(phylo4_)[getEdge(phylo4_,tipNode,"descendant")]
  }
  spltMat<-spltMat[order(spltMat$par),]
  return(list(split=spltMat,tip=tipMat))
}

# get_BHV_fromFixedTr<-function(fixedTr4, trTxt.ls){
#   outFile<-"tmp.txt"
#   BHV.vec<-rep(NA,length(trTxt.ls))
#   for(lIdx in 1:length(trTxt.ls)){
#     write.tree(unroot(as(fixedTr4,"phylo")),file=outFile,append=F)
#     txt.tr<-trTxt.ls[[lIdx]]
#     write.tree(unroot(read.tree(text = txt.tr)),file=outFile,append=T)
#     java_codeTxt<-paste0("java -jar /Users/shakingyao/Documents/Research_Veera/PYDT/phTree_CI/GTP/gtp.jar -u ",outFile)
#     system(java_codeTxt, ignore.stdout = T)
#     GTPOutFileName <- "./output.txt"
#     GTPFile <- file(GTPOutFileName,open="r")
#     GTPLines <-readLines(GTPFile)
#     close(GTPFile)
#     ele<-as.numeric(unlist(strsplit(GTPLines[1],split="\t") ))
#     BHV.vec[lIdx]<-ele[3]
#   }
#   return(BHV.vec)
# }

# get_BHV<-function(trueTr, trueTr_rtBr=0,
#                   estTr.txt, estTr_rtBr=0, 
#                   outFile="/home/yaots/Research/UltrametricsMat/tmpTxt/bhvOutFile.txt",
#                   trFile="/home/yaots/Research/UltrametricsMat/tmpTxt/tmp.txt",
#                   jarFile="/home/yaots/Research/GTP/gtp.jar",
#                   javaExec="/home/yaots/jdk-20.0.1/bin/java"){
#   # add root to true_tree
#   trueTr.txt<-write.tree(as(trueTr,"phylo"))
#   end.idx<-str_locate(string = trueTr.txt, "\\);")
#   end2.idx<-str_locate(string = trueTr.txt, "\\)0\\);")
#   if(!any(is.na(end.idx))){
#     trueTr.txt<-paste0("(",str_replace(string = trueTr.txt, pattern = "\\);", paste0("):",trueTr_rtBr,");")))
#   }else if(any(is.na(end2.idx))){
#     trueTr.txt<-paste0(str_replace(string = trueTr.txt, pattern = "\\):0\\);", paste0("):",trueTr_rtBr,");")))
#   }
#   
#   end.idx<-str_locate(string = estTr.txt, "\\):0\\);")
#   if(any(is.na(end.idx))){
#     stop("wrong estTr.txt format: no root specified in text \\):0\\);")
#   }else{
#     if(estTr_rtBr>0){
#       # estTr.txt<-paste0("(",str_replace(string = estTr.txt, pattern = "\\):0;", paste0("):",estTr_rtBr,");")))
#       estTr.txt<-paste0(str_replace(string = estTr.txt, pattern = "\\):0\\);", paste0("):",estTr_rtBr,");")))
#     }
#   }
#   
#   trF<-file(trFile)
#   writeLines(text = c(trueTr.txt,estTr.txt),con = trF)
#   close(trF)
#   
#   # java_codeTxt<-paste(javaExec,"-jar",jarFile,"-v -o",outFile,trFile)
#   # system(java_codeTxt,ignore.stdout = F)
#   java_codeTxt<-paste(javaExec,"-jar",jarFile,"-o",outFile,trFile)
#   system(java_codeTxt,ignore.stdout = T)
#   GTPF <- file(outFile,open="r")
#   GTPLines <-readLines(GTPF)
#   close(GTPF)
#   ele<-as.numeric(unlist(strsplit(GTPLines[1],split="\t") ))
#   BHV<-ele[3]
#   return(BHV)
# }

# trCov2Splt<-function(inp.trCov, tol=10^{-10}){
#   wrkTrCov<-inp.trCov
#   nLeaf<-nrow(inp.trCov)
#   rownames(wrkTrCov)<-colnames(wrkTrCov)<-paste0("t",1:nLeaf)
#   int.count<-1
#   cur.trCov.list<-list()
#   res.splt<-data.frame(matrix(NA,nr=nLeaf-2,nc=2))
#   res.tip<-data.frame(matrix(NA,nr=nLeaf,nc=2))
#   colnames(res.splt)<-colnames(res.tip)<-c("par","length")
#   res.tip$par<-as.character(1:nLeaf)
#   
#   
#   singRootLn=min(wrkTrCov)
#   wrkTrCov<-wrkTrCov - singRootLn
#   l<-rownames(wrkTrCov)[1]
#   desc1<-names(which(wrkTrCov[rownames(wrkTrCov)==l,]==0))
#   idx1<-rownames(wrkTrCov)%in%desc1
#   subTrCov1<-as.matrix(wrkTrCov[idx1,idx1])
#   if(nrow(subTrCov1)==1){
#     res.tip[res.tip$par==str_sub(desc1,2),2]<-as.numeric(subTrCov1)
#   }else{
#     len1<-min(subTrCov1[lower.tri(subTrCov1)])
#     subTrCov1<-subTrCov1-len1
#     subTrCov1[subTrCov1<tol]<-0
#     res.splt[int.count,1]<-paste(sort(as.numeric(str_sub(desc1,2))),collapse = ",")
#     res.splt[int.count,2]<-as.numeric(len1)
#     int.count<-int.count+1
#     if(nrow(subTrCov1)==2){
#       l1.idx<-rownames(subTrCov1)==desc1[1]
#       l2.idx<-rownames(subTrCov1)==desc1[2]
#       res.tip[res.tip$par==str_sub(desc1[1],2),2]<-as.numeric(subTrCov1[l1.idx,l1.idx])
#       res.tip[res.tip$par==str_sub(desc1[2],2),2]<-as.numeric(subTrCov1[l2.idx,l2.idx])
#     }else{
#       cur.trCov.list[[length(cur.trCov.list)+1]]<-subTrCov1
#     }
#   } 
#   
#   desc2<-setdiff(rownames(wrkTrCov),desc1)
#   idx2<-rownames(wrkTrCov)%in%desc2
#   subTrCov2<-as.matrix(wrkTrCov[idx2,idx2])
#   if(nrow(subTrCov2)==1){
#     res.tip[res.tip$par==str_sub(desc2,2),2]<-as.numeric(subTrCov2)
#   }else{
#     len2<-min(subTrCov2[lower.tri(subTrCov2)])
#     subTrCov2<-subTrCov2-len2
#     subTrCov2[subTrCov2<tol]<-0
#     res.splt[int.count,1]<-paste(sort(as.numeric(str_sub(desc2,2))),collapse = ",")
#     res.splt[int.count,2]<-as.numeric(len2)
#     int.count<-int.count+1
#     if(nrow(subTrCov2)==2){
#       l1.idx<-rownames(subTrCov2)==desc2[1]
#       l2.idx<-rownames(subTrCov2)==desc2[2]
#       res.tip[res.tip$par==str_sub(desc2[1],2),2]<-as.numeric(subTrCov2[l1.idx,l1.idx])
#       res.tip[res.tip$par==str_sub(desc2[2],2),2]<-as.numeric(subTrCov2[l2.idx,l2.idx])
#     }else{
#       cur.trCov.list[[length(cur.trCov.list)+1]]<-subTrCov2
#     }
#   } 
#   
#   while(length(cur.trCov.list)>0){
#     wrk.trCov.ls<-list()
#     for(mIdx in 1:length(cur.trCov.list)){
#       cur.TrCov<-cur.trCov.list[[mIdx]]
#       
#       l<-rownames(cur.TrCov)[1]
#       desc1<-names(which(cur.TrCov[rownames(cur.TrCov)==l,]==0))
#       idx1<-rownames(cur.TrCov)%in%desc1
#       subTrCov1<-as.matrix(cur.TrCov[idx1,idx1])
#       if(nrow(subTrCov1)==1){
#         res.tip[res.tip$par==str_sub(desc1,2),2]<-as.numeric(subTrCov1)
#       }else{
#         len1<-min(subTrCov1[lower.tri(subTrCov1)])
#         subTrCov1<-subTrCov1-len1
#         subTrCov1[subTrCov1<tol]<-0
#         res.splt[int.count,1]<-paste(sort(as.numeric(str_sub(desc1,2))),collapse = ",")
#         res.splt[int.count,2]<-as.numeric(len1)
#         int.count<-int.count+1
#         if(nrow(subTrCov1)==2){
#           l1.idx<-rownames(subTrCov1)==desc1[1]
#           l2.idx<-rownames(subTrCov1)==desc1[2]
#           res.tip[res.tip$par==str_sub(desc1[1],2),2]<-as.numeric(subTrCov1[l1.idx,l1.idx])
#           res.tip[res.tip$par==str_sub(desc1[2],2),2]<-as.numeric(subTrCov1[l2.idx,l2.idx])
#         }else{
#           wrk.trCov.ls[[length(wrk.trCov.ls)+1]]<-subTrCov1
#         }
#       }
#       
#       desc2<-setdiff(rownames(cur.TrCov),desc1)
#       idx2<-rownames(cur.TrCov)%in%desc2
#       subTrCov2<-as.matrix(cur.TrCov[idx2,idx2])
#       if(nrow(subTrCov2)==1){
#         res.tip[res.tip$par==str_sub(desc2,2),2]<-as.numeric(subTrCov2)
#       }else{
#         len2<-min(subTrCov2[lower.tri(subTrCov2)])
#         subTrCov2<-subTrCov2-len2
#         subTrCov2[subTrCov2<tol]<-0
#         res.splt[int.count,1]<-paste(sort(as.numeric(str_sub(desc2,2))),collapse = ",")
#         res.splt[int.count,2]<-as.numeric(len2)
#         int.count<-int.count+1
#         if(nrow(subTrCov2)==2){
#           l1.idx<-rownames(subTrCov2)==desc2[1]
#           l2.idx<-rownames(subTrCov2)==desc2[2]
#           res.tip[res.tip$par==str_sub(desc2[1],2),2]<-as.numeric(subTrCov2[l1.idx,l1.idx])
#           res.tip[res.tip$par==str_sub(desc2[2],2),2]<-as.numeric(subTrCov2[l2.idx,l2.idx])
#         }else{
#           wrk.trCov.ls[[length(wrk.trCov.ls)+1]]<-subTrCov2
#         }
#       }
#     }
#     cur.trCov.list<-wrk.trCov.ls
#   }
#   return(list(splt=res.splt,tip=res.tip, singRootLn= singRootLn))
# }

get_BHV_intSplt<-function(inp.splt1, inp.splt2, nLeaf,
                          outFile="/home/yaots/Research/UltrametricsMat/tmpTxt/bhvOutFile.txt",
                          trFile="/home/yaots/Research/UltrametricsMat/tmpTxt/tmp.txt",
                          jarFile="/home/yaots/Research/GTP/gtp.jar",
                          javaExec="/home/yaots/jdk-20.0.1/bin/java"){
  tmp.leaf<-data.frame(par=1:nLeaf,length=1)
  tr1.txt<-trTbl2Tr(splt2TrTbl(inp.splt = inp.splt1, inp.pedEdge = tmp.leaf))$txt
  tr2.txt<-trTbl2Tr(splt2TrTbl(inp.splt = inp.splt2, inp.pedEdge = tmp.leaf))$txt
  
  trF<-file(trFile)
  writeLines(text = c(tr1.txt,tr2.txt),con = trF)
  close(trF)
  
  # java_codeTxt<-paste(javaExec,"-jar",jarFile,"-v -o",outFile,trFile)
  # system(java_codeTxt,ignore.stdout = F)
  java_codeTxt<-paste(javaExec,"-jar",jarFile,"-o",outFile,trFile)
  system(java_codeTxt,ignore.stdout = T)
  GTPF <- file(outFile,open="r")
  GTPLines <-readLines(GTPF)
  close(GTPF)
  ele<-as.numeric(unlist(strsplit(GTPLines[1],split="\t") ))
  BHV<-ele[3]
  return(BHV)
}

get_BHV_leaf<-function(inp.pedEdge1, rtEdge1, inp.pedEdge2, rtEdge2){
  inp.pedEdge1<-inp.pedEdge1[order(as.numeric(inp.pedEdge1$par)),]
  inp.pedEdge2<-inp.pedEdge2[order(as.numeric(inp.pedEdge2$par)),]
  leaf.vec1<-c(rtEdge1,inp.pedEdge1$length)
  leaf.vec2<-c(rtEdge2,inp.pedEdge2$length)
  diff.vec=abs(leaf.vec1-leaf.vec2)
  dist=sqrt(sum(diff.vec^2))
  return(dist)
}

trCov2Splt<-function(inp.trCov, method="Tree", tol=10^{-8}){
  if(method=="MIP"){
    tol=10^{-4}-10^{-8}
  }
  wrkTrCov<-inp.trCov
  nLeaf<-nrow(inp.trCov)
  rownames(wrkTrCov)<-colnames(wrkTrCov)<-paste0("t",1:nLeaf)
  int.count<-1
  cur.trCov.list<-list()
  res.splt<-data.frame(matrix(NA,nr=nLeaf-2,nc=2))
  res.tip<-data.frame(matrix(NA,nr=nLeaf,nc=2))
  colnames(res.splt)<-colnames(res.tip)<-c("par","length")
  res.tip$par<-as.character(1:nLeaf)
  
  desc.vec<-character(0)
  singRootLn=min(wrkTrCov)
  if(singRootLn < tol){
    singRootLn<-0
  }
  wrkTrCov<-wrkTrCov - singRootLn
  wrkTrCov[wrkTrCov<tol]<-0
  for(lIdx in 1:nLeaf){
    l<-rownames(wrkTrCov)[lIdx]
    desc.tmp<-paste(names(which(wrkTrCov[rownames(wrkTrCov)==l,]>tol)),collapse = ",")
    desc.vec<-c(desc.vec,desc.tmp)
  }
  desc.vec<-unique(desc.vec)
  
  for(dIdx in 1:length(desc.vec)){
    desc<-str_split_1(desc.vec[dIdx],",")
    idx<-rownames(wrkTrCov)%in%desc
    subTrCov<-as.matrix(wrkTrCov[idx,idx])
    if(nrow(subTrCov)==1){
      res.tip[res.tip$par==str_sub(desc,2),2]<-as.numeric(subTrCov)
    }else{
      len<-min(subTrCov[lower.tri(subTrCov)])
      subTrCov<-subTrCov-len
      subTrCov[subTrCov<tol]<-0
      res.splt[int.count,1]<-paste(sort(as.numeric(str_sub(desc,2))),collapse = ",")
      res.splt[int.count,2]<-as.numeric(len)
      int.count<-int.count+1
      cur.trCov.list[[length(cur.trCov.list)+1]]<-subTrCov
    }
  }
  
  while(length(cur.trCov.list)>0){
    wrk.trCov.ls<-list()
    for(mIdx in 1:length(cur.trCov.list)){
      cur.TrCov<-cur.trCov.list[[mIdx]]
      
      desc.vec<-character(0)
      for(lIdx in 1:nrow(cur.TrCov)){
        l<-rownames(cur.TrCov)[lIdx]
        desc.tmp<-paste(names(which(cur.TrCov[rownames(cur.TrCov)==l,]>tol)),collapse = ",")
        desc.vec<-c(desc.vec,desc.tmp)
      }
      desc.vec<-unique(desc.vec)
      if(any(duplicated(str_split_1(paste(desc.vec,collapse=","),",")))){
        return(NA)
      }
      
      for(dIdx in 1:length(desc.vec)){
        desc<-str_split_1(desc.vec[dIdx],",")
        idx<-rownames(cur.TrCov)%in%desc
        subTrCov<-as.matrix(cur.TrCov[idx,idx])
        if(nrow(subTrCov)==1){
          res.tip[res.tip$par==str_sub(desc,2),2]<-as.numeric(subTrCov)
        }else{
          len<-min(subTrCov[lower.tri(subTrCov)])
          subTrCov<-subTrCov-len
          subTrCov[subTrCov<tol]<-0
          res.splt[int.count,1]<-paste(sort(as.numeric(str_sub(desc,2))),collapse = ",")
          res.splt[int.count,2]<-as.numeric(len)
          int.count<-int.count+1
          wrk.trCov.ls[[length(wrk.trCov.ls)+1]]<-subTrCov
        }
      }
    }
    cur.trCov.list<-wrk.trCov.ls
  }
  res.splt<-res.splt[!is.na(res.splt$par),]
  return(list(splt=res.splt,tip=res.tip, singRootLn=singRootLn))
}

MCMC_Tr<-function(Obsdf_, tipBr.sd_=0.5, intBr.sd_=0.5, betaSplt.tbl_, rtEdge=T, 
                  init.Splt=NULL, init.Tip=NULL, init.RtEdge=NULL,
                  ctsSplt.sd_=0, hypParam.brLen_=1,hypParam.tplg_=-1.5,
                  compMat_=NULL, seed_=1, iteNum=3000, burnIn=500){
  set.seed(seed_)
  nLeaf<-ncol(Obsdf_)
  
  if(rtEdge){
    if(is.null(init.RtEdge)){
      init.RtEdge<-rexp(1)
    }
    init.rtEdge<-init.RtEdge
  }else{
    init.rtEdge<-0
  }
  
  # initiation from rtree
  if(is.null(init.Splt) | is.null(init.Tip)){
    init.Tr.phylo4<-phylo4(rtree(nLeaf))
    init.Tr<-phylo42Splt(init.Tr.phylo4)
    init.Splt<-init.Tr$split
    init.Tip<-init.Tr$tip
  }else{
    init.Tr<-list()
    init.Tr$split<-init.Splt
    init.Tr$tip<-init.Tip
  }
  init.Splt<-init.Splt[order(init.Splt$par),]
  init.Tr$split<-init.Splt
  init.Tr$rtEdge<-init.rtEdge
  
  llh.vec<-rep(NA,iteNum)
  acpt.vec<-rep(0,1+nLeaf+1)
  names(acpt.vec)<-c("Splt",paste0("tipBrLen",1:nLeaf),"rtEdgeLen")
  splt.ls<-tip.ls<-list()
  rtEdge.vec<-rep(NA,iteNum-burnIn)
  
  if(is.null(betaSplt.tbl_)){
    betaSplt.tbl_<-betaSplt.tbl.gen(beta_ = hypParam.tplg)
  }
  
  if(length(ctsSplt.sd_)==1){
    ctsSplt.sd.vec<-rep(ctsSplt.sd_,iteNum)
  }else{
    if(length(ctsSplt.sd_)!=iteNum){
      stop("Wrong length of ctsSplt.sd_")
    }else{
      ctsSplt.sd.vec<-ctsSplt.sd_
    }
  } 
  
  TrSplt<-init.Splt
  TrTip<-init.Tip
  TrRtEdge<-init.rtEdge
  init.TrCov<-splt2TrCov(inp.splt = TrSplt, inp.pedEdge = TrTip, inp.rtEdge=init.rtEdge)
  init.llh<-sum(dmvnorm(Obsdf_,sigma = init.TrCov,log = T))
  
  for(ite in 1:iteNum){
    if(ite %% 1000==0){print(ite)}
    ctsSplt.sd.num<-ctsSplt.sd.vec[ite]
    newSplt<-update_SpltID(inp.curSplt = TrSplt, inp.tip = TrTip, inp.rtEdge=TrRtEdge,
                           inp.Obsdf = Obsdf_, inp.betaSplt.tbl = betaSplt.tbl_,
                           inp.compMat = compMat_, hypParam_mean = hypParam.brLen_,
                           nLeaf=nLeaf, inp.ctsSplt.sd=ctsSplt.sd.num)
    
    # plot(Splt2Phylo4(inp.spltID = TrSplt, inp.dfSplt = df.splt, inp.pedEdge = TrTip))
    TrSplt<-newSplt$split
    acpt.vec["Splt"]<-acpt.vec["Splt"]+newSplt$acpt
    
    newIntSplt<-update_IntLenTr(inp.curSplt = TrSplt, inp.tip = TrTip, nLeaf = nLeaf, 
                                inp.rtEdge=TrRtEdge, hypParam_mean = hypParam.brLen_,
                                inp.Obsdf = Obsdf_, inp.intBr.sd = 0.5)
    TrSplt<-newIntSplt$intBr
    
    newTip<-update_TipLenTr(inp.curSplt = TrSplt, inp.tip = TrTip, nLeaf = nLeaf,
                            inp.rtEdge=TrRtEdge, hypParam_mean = hypParam.brLen_,
                            inp.Obsdf = Obsdf_, inp.tipBr.sd = 0.5)
    TrTip<-newTip$tipBr
    acpt.vec[paste0("tipBrLen",1:nLeaf)]<-acpt.vec[paste0("tipBrLen",1:nLeaf)]+newTip$acpt
    
    if(rtEdge){
      newRtEdge<-update_rtEdgeLen(inp.curSplt = TrSplt, inp.tip = TrTip, nLeaf = nLeaf,
                                  inp.rtEdge=TrRtEdge, hypParam_mean = hypParam.brLen_,
                                  inp.Obsdf = Obsdf_, inp.tipBr.sd = 0.5)
    }
    TrRtEdge<-newRtEdge$rtEdgeBr
    acpt.vec["rtEdgeLen"]<-acpt.vec["rtEdgeLen"]+newRtEdge$acpt
    
    llh.vec[ite]<-newRtEdge$llh
    if(ite>burnIn){
      splt.ls[[ite-burnIn]]<-TrSplt
      tip.ls[[ite-burnIn]]<-TrTip
      rtEdge.vec[ite-burnIn]<-TrRtEdge
    }
  }
  return(list(llh=c(init.llh,llh.vec), acpt=acpt.vec, rtEdge=rtEdge.vec,
              split=splt.ls, tip=tip.ls, init.Tr=init.Tr))
}

MCMC_NyeAlg<-function(Obsdf_, tipBr.sd_=0.5, ctsSplt.sd_=0.5,
                      betaSplt.tbl_, rtEdge=T, 
                      init.Splt=NULL, init.Tip=NULL, init.RtEdge=NULL,
                      hypParam.brLen_=1,hypParam.tplg_=-1.5,
                      compMat_=NULL, seed_=1, iteNum=3000, burnIn=500){
  set.seed(seed_)
  nLeaf<-ncol(Obsdf_)
  
  if(rtEdge){
    if(is.null(init.RtEdge)){
      init.RtEdge<-rexp(1)
    }
    init.rtEdge<-init.RtEdge
  }else{
    init.rtEdge<-0
  }
  
  # initiation from rtree
  if(is.null(init.Splt) | is.null(init.Tip)){
    init.Tr.phylo4<-phylo4(rtree(nLeaf))
    init.Tr<-phylo42Splt(init.Tr.phylo4)
    init.Splt<-init.Tr$split
    init.Tip<-init.Tr$tip
  }else{
    init.Tr<-list()
    init.Tr$split<-init.Splt
    init.Tr$tip<-init.Tip
  }
  init.Splt<-init.Splt[order(init.Splt$par),]
  init.Tr$split<-init.Splt
  init.Tr$rtEdge<-init.rtEdge
  
  llh.vec<-rep(NA,iteNum)
  acpt.vec<-rep(0,1+nLeaf+1)
  names(acpt.vec)<-c("Splt",paste0("tipBrLen",1:nLeaf),"rtEdgeLen")
  splt.ls<-tip.ls<-list()
  rtEdge.vec<-rep(NA,iteNum-burnIn)
  
  if(is.null(betaSplt.tbl_)){
    betaSplt.tbl_<-betaSplt.tbl.gen(beta_ = hypParam.tplg)
  }
  
  if(length(ctsSplt.sd_)==1){
    ctsSplt.sd.vec<-rep(ctsSplt.sd_,iteNum)
  }else{
    if(length(ctsSplt.sd_)!=iteNum){
      stop("Wrong length of ctsSplt.sd_")
    }else{
      ctsSplt.sd.vec<-ctsSplt.sd_
    }
  } 
  
  TrSplt<-init.Splt
  TrTip<-init.Tip
  TrRtEdge<-init.rtEdge
  init.TrCov<-splt2TrCov(inp.splt = TrSplt, inp.pedEdge = TrTip, inp.rtEdge=init.rtEdge)
  init.llh<-sum(dmvnorm(Obsdf_,sigma = init.TrCov,log = T))
  
  for(ite in 1:iteNum){
    if(ite %% 1000==0){print(ite)}
    ctsSplt.sd.num<-ctsSplt.sd.vec[ite]
    newSplt<-update_SpltID(inp.curSplt = TrSplt, inp.tip = TrTip, inp.rtEdge=TrRtEdge,
                           inp.Obsdf = Obsdf_, inp.betaSplt.tbl = betaSplt.tbl_,
                           inp.compMat = compMat_, hypParam_mean = hypParam.brLen_,
                           nLeaf=nLeaf, inp.ctsSplt.sd=ctsSplt.sd.num)
    
    # plot(Splt2Phylo4(inp.spltID = TrSplt, inp.dfSplt = df.splt, inp.pedEdge = TrTip))
    TrSplt<-newSplt$split
    acpt.vec["Splt"]<-acpt.vec["Splt"]+newSplt$acpt
    
    newTip<-update_TipLenTr(inp.curSplt = TrSplt, inp.tip = TrTip, nLeaf = nLeaf,
                            inp.rtEdge=TrRtEdge, hypParam_mean = hypParam.brLen_,
                            inp.Obsdf = Obsdf_, inp.tipBr.sd = 0.5)
    TrTip<-newTip$tipBr
    acpt.vec[paste0("tipBrLen",1:nLeaf)]<-acpt.vec[paste0("tipBrLen",1:nLeaf)]+newTip$acpt
    
    if(rtEdge){
      newRtEdge<-update_rtEdgeLen(inp.curSplt = TrSplt, inp.tip = TrTip, nLeaf = nLeaf,
                                  inp.rtEdge=TrRtEdge, hypParam_mean = hypParam.brLen_,
                                  inp.Obsdf = Obsdf_, inp.tipBr.sd = 0.5)
    }
    TrRtEdge<-newRtEdge$rtEdgeBr
    acpt.vec["rtEdgeLen"]<-acpt.vec["rtEdgeLen"]+newRtEdge$acpt
    
    llh.vec[ite]<-newRtEdge$llh
    if(ite>burnIn){
      splt.ls[[ite-burnIn]]<-TrSplt
      tip.ls[[ite-burnIn]]<-TrTip
      rtEdge.vec[ite-burnIn]<-TrRtEdge
    }
  }
  return(list(llh=c(init.llh,llh.vec), acpt=acpt.vec, rtEdge=rtEdge.vec,
              split=splt.ls, tip=tip.ls, init.Tr=init.Tr))
}


# postMeanTr<-function(splt.ls, tip.ls,thinIn=10, 
#                      trOutFile="/home/yaots/Research/UltrametricsMat/tmpTxt/tr.txt", 
#                      meanOutFile="/home/yaots/Research/UltrametricsMat/tmpTxt/meanOutFile.txt", 
#                      jarFile="/home/yaots/Research/SturmMean/SturmMean.jar",
#                      javaExec="/home/yaots/jdk-20.0.1/bin/java"){
#   postSmpSize<-length(splt.ls)
#   smpIdx<-seq(thinIn,postSmpSize,thinIn)
#   splt.ls.in<-splt.ls[smpIdx]
#   tip.ls.in<-tip.ls[smpIdx]
#   trTxt.vec<-character()
#   for(sIdx in 1:length(smpIdx)){
#     trTbl<-splt2TrTbl(inp.splt = splt.ls.in[[sIdx]], inp.pedEdge = tip.ls.in[[sIdx]])
#     tr.ls<-trTbl2Tr(inp.trTbl = trTbl)
#     trTxt<-tr.ls$txt
#     # rtTrTxt<-str_replace(trTxt,":0;",paste0(":",rtEdge.vec.in[sIdx],";"))
#     trTxt.vec<-c(trTxt.vec, trTxt)
#   }
#   trF<-file(trOutFile)
#   writeLines(text = trTxt.vec,con = trF)
#   close(trF)
#   
#   java_codeTxt<-paste(javaExec,"-jar",jarFile,"-o",meanOutFile,trOutFile)
#   system(java_codeTxt,ignore.stdout = T)
#   MeanTrF <- file(meanOutFile,open="r")
#   MeanTrLines <-readLines(MeanTrF)
#   MeanTrTxt<-MeanTrLines[5]
#   close(MeanTrF)
#   # add root node
#   end.idx<-str_locate(string = MeanTrTxt, "\\);")
#   if(!any(is.na(end.idx))){
#     # MeanTrTxt<-paste0("(",str_replace(string = MeanTrTxt, pattern = "\\);", paste0("):",0,");")))
#     MeanTrTxt<-paste0(str_replace(string = MeanTrTxt, pattern = "\\);", paste0("):",0,";")))
#   }
#   return(list(txt=MeanTrTxt, phylo4=phylo4(read.tree(text=MeanTrTxt))))
# }

meanRtEdgeTr<-function(splt.ls, tip.ls, rtEdge.vec, thinIn=10, nLeaf,
                       trOutFile="/home/yaots/Research/UltrametricsMat/tmpTxt/tr.txt", 
                       meanOutFile="/home/yaots/Research/UltrametricsMat/tmpTxt/meanOutFile.txt", 
                       jarFile="/home/yaots/Research/SturmMean/SturmMean.jar",
                       javaExec="/home/yaots/jdk-20.0.1/bin/java"){
  postSmpSize<-length(splt.ls)
  smpIdx<-seq(thinIn,postSmpSize,thinIn)
  splt.ls.in<-splt.ls[smpIdx]
  tip.ls.in<-tip.ls[smpIdx]
  rtEdge.vec.in<-rtEdge.vec[smpIdx]
  trTxt.vec<-character()
  tmp.tipEdge<-data.frame(par=1:nLeaf,length=1)
  tmp.tipMat<-matrix(NA,nrow=nLeaf,ncol=length(smpIdx))
  for(sIdx in 1:length(smpIdx)){
    trTbl<-splt2TrTbl(inp.splt = splt.ls.in[[sIdx]], inp.pedEdge = tmp.tipEdge)
    tr.ls<-trTbl2Tr(inp.trTbl = trTbl)
    trTxt<-tr.ls$txt
    # rtTrTxt<-str_replace(trTxt,":0;",paste0(":",rtEdge.vec.in[sIdx],";"))
    trTxt.vec<-c(trTxt.vec, trTxt)
    tmp.tip<-tip.ls.in[[sIdx]]
    tmp.tip<-tmp.tip[order(as.numeric(tmp.tip$par)),]
    tmp.tipMat[,sIdx]<-tmp.tip$length
  }
  meanRtEdge<-mean(rtEdge.vec.in)
  mean.tip<-data.frame(par=1:nLeaf,length=apply(tmp.tipMat,1,mean))
  
  trF<-file(trOutFile)
  writeLines(text = trTxt.vec,con = trF)
  close(trF)
  
  java_codeTxt<-paste(javaExec,"-jar",jarFile,"-o",meanOutFile,trOutFile)
  system(java_codeTxt,ignore.stdout = T)
  MeanTrF <- file(meanOutFile,open="r")
  MeanTrLines <-readLines(MeanTrF)
  MeanTrTxt_int<-MeanTrLines[5]
  close(MeanTrF)
  # add root node
  end.idx<-str_locate(string = MeanTrTxt_int, "\\);")
  if(!any(is.na(end.idx))){
    # MeanTrTxt<-paste0("(",str_replace(string = MeanTrTxt, pattern = "\\);", paste0("):",0,");")))
    MeanTrTxt_int<-paste0(str_replace(string = MeanTrTxt_int, pattern = "\\);", paste0("):",0,";")))
  }
  MeanTrTxt<-MeanTrTxt_int
  for(i in 1:nLeaf){
    leaf.idx<-str_locate(string = MeanTrTxt,pattern = paste0("t",i,":1"))
    if(is.na(all(leaf.idx))){
      stop(paste0("No target leaf: ", i))
    }
    MeanTrTxt<-str_replace(string = MeanTrTxt,pattern = paste0("t",i,":1"),
                           replacement = paste0("t",i,":",mean.tip$length[i]))
  }
  
  return(list(int_txt=MeanTrTxt_int, final.txt=MeanTrTxt,
              phylo4=phylo4(read.tree(text=MeanTrTxt)),rtEdge=meanRtEdge))
}

phylo42TrCov<-function(inp.phylo4){
  nLeaf<-nTips(inp.phylo4)
  resTrCov<-matrix(NA,nr=nLeaf,nc=nLeaf)
  for(j in 1:nLeaf){
    for(k in 1:nLeaf){
      resTrCov[j,k]<-resTrCov[k,j]<-nodeHeight(inp.phylo4,node = MRCA(inp.phylo4,paste0("t",c(j,k))),"root")
    }
  }
  return(resTrCov)
}