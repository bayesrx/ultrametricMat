library(partitions)
library(stringr)
library(ape)
library(phylobase)
library(mvtnorm)
library(truncnorm)
library(ggplot2)
library(here)
source(here("ultrametricMat","Code","MCMC_BHV_Func.R"))

# args <- commandArgs(TRUE)
args<-c("rtEdge_rtree","10","1","1","1")
print(args)
trueTrType=trimws(args[1])
nLeaf<-as.numeric(args[2])
nIdx<-as.numeric(args[3])
lIdx<-as.numeric(args[4])
rIdx<-as.numeric(args[5])

trueList<-readRDS(here("ultrametricMat","Data",paste0("",trueTrType,"_nLeaf",nLeaf,".RDS")))
df.ls<-trueList$df
nSmp.len<-length(df.ls)
nRep<-length(df.ls[[1]])
nSmp.vec<-unlist(lapply(df.ls,function(x) nrow(x[[1]]$n.ObsDf)))
llh.typeVec<-c("n","t4","t3")

nSmp<-nSmp.vec[nIdx]
llh.type<-llh.typeVec[lIdx]

# beta_=-1.5
# betaSplt.tbl<-betaSplt.tbl.gen(n=nLeaf,beta_=beta_)
# saveRDS(betaSplt.tbl, file=here("ultrametricMat","MCMC_Res",paste0("betaSpltTbl_unif_nLeaf",nLeaf,".RDS")))
betaSplt.tbl<-readRDS(here("ultrametricMat","MCMC_Res",paste0("betaSpltTbl_unif_nLeaf",nLeaf,".RDS")))


#### simulation start ####
iteNum<-10000
burnIn<-9000


obsDf<-df.ls[[paste0("nSmp_",nSmp)]][[rIdx]][[paste0(llh.type,".","ObsDf")]]
print(paste0("nSmp: ",nSmp,"; llhType: ", llh.type,"; rep: ", rIdx))
NNI.t<-system.time({NNI_dis<-MCMC_Tr(Obsdf_ = obsDf, betaSplt.tbl_ = betaSplt.tbl, rtEdge=T,
                                     compMat_ = NULL, ctsSplt.sd_ = 0, seed_=1+10*nIdx+100*lIdx+1000*rIdx, 
                                     iteNum = iteNum, burnIn=burnIn)
})
# print(NNI.t)
# ggplot(data.frame(ite=0:iteNum,llh=NNI_dis$llh),aes(x=ite,y=llh))+geom_line()
# View(phylo42Splt(trueList$trueTr)$split)
# View(NNI_dis$split[[iteNum-burnIn]])
# saveRDS(list(res=NNI_dis,t=NNI.t[1:3]),
#         paste0("./Res/",trueTrType,"/MCMC/nLeaf",nLeaf,"_nSmp",nSmp,"_llh",llh.type,"_rep",rIdx,".RDS"))
saveRDS(list(res=NNI_dis,t=NNI.t[1:3]),
        file=here("UltrametricMat","MCMC_Res",paste0(trueTrType,"_nTr_nLeaf",nLeaf,"_nSmp",nSmp,"_llh",llh.type,"_rep",rIdx,".RDS")))


