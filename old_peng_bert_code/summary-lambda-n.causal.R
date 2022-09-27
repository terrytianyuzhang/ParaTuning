rm(list=ls()); gc()
options(stringsAsFactors = F)

#### software needed
plink="/data3/Software/Plink/plink"
plink2="/data3/Software/Plink2/plink2"

#### load the functions that are needed
source("R-code/simulation-functions.R")

library(data.table)
library(pryr) # check memory useage
library(lassosum) #transform p value to correlation
library(doParallel) # foreach
library(pROC) # for AUC 
library(R.utils)
library(snpStats)
library(wordspace)

RE.auc=NULL
for(i.sim in c("NC2000","NC4000","NC8000")){
  work.dir=paste0("Work/Sim-C0.80-Y0.60-rep",i.sim,"/")
  ped=list()
  ped[["CEU"]]=list()
  for(N in seq(20000,20000,4000)){
    ped[["CEU"]][[paste0("N",N)]]=fread(paste0(work.dir,"CEU.TRN/CEU.TRN-",N,".psam"),header=T,data.table=F)
    rownames(ped[["CEU"]][[paste0("N",N)]])=ped[["CEU"]][[paste0("N",N)]]$IID
  }
  ped[["YRI"]]=list()
  for(N in seq(20000,20000,4000)){
    ped[["YRI"]][[paste0("N",N)]]=fread(paste0(work.dir,"YRI.TRN/YRI.TRN-",N,".psam"),header=T,data.table=F)
    rownames(ped[["YRI"]][[paste0("N",N)]])=ped[["YRI"]][[paste0("N",N)]]$IID
  }
  
  gamma=c("0.50")
  for(trn.set in c("CEU.TRN","YRI.TRN")){
    print(trn.set); flush.console()
    anc=substr(trn.set,1,3)
    pgs =fread(paste0(work.dir,"CombinedLassoSum/PGS-",trn.set,"-20000-lasso-C20000-Y20000-gamma-",gamma,".pgs"),header=T,data.table=F)
    rownames(pgs)=pgs$V1
    pgs=pgs[,-1]
    for(N in names(ped[[anc]])){
      re.auc=NULL
      for(lambda in colnames(pgs)){
        re.auc=c(re.auc,auc(ped[[anc]][[N]]$PHENO1,pgs[rownames(ped[[anc]][[N]]),lambda])[1])
      }
      names(re.auc)=colnames(pgs)
    }
    RE.auc=rbind.data.frame(RE.auc,data.frame(i.sim,gwas.CEU=20000,gwas.YRI=20000,gamma="0.50",anc,N,t(re.auc)))
    
    
  }
  
}

RE.auc
