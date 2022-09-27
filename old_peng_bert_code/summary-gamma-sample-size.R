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

#### load the parameter information
i.sim=1
load("Work/Sim-C0.80-Y0.60-rep1/simulation-params.RData")

# set directories for this simulation
main.dir=params$run.info$main.dir
work.dir=params$run.info$work.dir

### read the pedigree files for the validation sets
ped=list()
ped[["CEU"]]=list()
for(N in seq(4000,20000,4000)){
  ped[["CEU"]][[paste0("N",N)]]=fread(paste0(work.dir,"CEU.TRN/CEU.TRN-",N,".psam"),header=T,data.table=F)
  rownames(ped[["CEU"]][[paste0("N",N)]])=ped[["CEU"]][[paste0("N",N)]]$IID
}
ped[["YRI"]]=list()
for(N in seq(4000,20000,4000)){
  ped[["YRI"]][[paste0("N",N)]]=fread(paste0(work.dir,"YRI.TRN/YRI.TRN-",N,".psam"),header=T,data.table=F)
  rownames(ped[["YRI"]][[paste0("N",N)]])=ped[["YRI"]][[paste0("N",N)]]$IID
}

### list of gammas
gamma=seq(0.2,0.8,length.out=5)


# find the files
files=list.files(paste0(work.dir,"CombinedLassoSum"),pattern = ".pgs")
file.df=data.frame(gwas.CEU=unlist(lapply(strsplit(files,"-"),`[[`,5)),
           gwas.YRI=unlist(lapply(strsplit(files,"-"),`[[`,6)),
           trn.pop=substr(unlist(lapply(strsplit(files,"-"),`[[`,2)),1,3),
           trn.N=unlist(lapply(strsplit(files,"-"),`[[`,3)),
           gamma=unlist(lapply(strsplit(files,"-"),`[[`,8)),
           file=files)
file.df[,"lambda:0.0010"]=NA
file.df[,"lambda:0.0027"]=NA
file.df[,"lambda:0.0071"]=NA
file.df[,"lambda:0.0188"]=NA
file.df[,"lambda:0.0500"]=NA


RE.auc=NULL
for(i in 1:nrow(file.df)){
  pgs=fread(paste0(work.dir,"CombinedLassoSum/",file.df$file[i]),header=T,data.table=F)
  anc=file.df$trn.pop[i]
  rownames(pgs)=pgs$V1
  pgs=pgs[,-1]
  
  for(N in paste0("N",seq(4000,20000,4000))){
    # run AUC on the European data
    re.auc=NULL
    for(lambda in colnames(pgs)){
      re.auc=c(re.auc,auc(ped[[anc]][[N]]$PHENO1,pgs[rownames(ped[[anc]][[N]]),lambda])[1])
    }
    names(re.auc)=colnames(pgs)
    RE.auc=rbind.data.frame(RE.auc,data.frame(gwas.CEU=file.df$gwas.CEU[i],gwas.YRI=file.df$gwas.YRI[i],gamma=file.df$gamma[i],anc,N,t(re.auc)))
  }
}

fwrite(RE.auc,"results-sample-size-testing.txt",row.names=F,col.names=T,quote=F,sep="\t")


gamma=seq(0.2,0.8,length=5)
lambda=c(0.0010,0.0027,0.0071,0.0188,0.0500)
N.trn=20000
anc="CEU"
RE.CEU=NULL
for(n.ceu in seq(4000,20000,4000)){
  for(n.yri in seq(4000,20000,4000)){
    sub.df=RE.auc[RE.auc$gwas.CEU == paste0("C",n.ceu) & RE.auc$gwas.YRI == paste0("Y",n.yri) & RE.auc$N == paste0("N",N.trn) & RE.auc$anc == anc,]
    ij=which(sub.df[,6:10] == max(sub.df[,6:10]),arr.ind=T)
    RE.CEU=rbind.data.frame(RE.CEU,data.frame(anc,gwas.CEU = n.ceu,gwas.YRI = n.yri, N=N.trn,gamma=gamma[ij[,1]],lambda=lambda[ij[,2]],auc=round(sub.df[ij[,1],ij[,2]+5],3)))
  }
}
RE.CEU

anc="YRI"
RE.YRI=NULL
for(n.ceu in seq(4000,20000,4000)){
  for(n.yri in seq(4000,20000,4000)){
    sub.df=RE.auc[RE.auc$gwas.CEU == paste0("C",n.ceu) & RE.auc$gwas.YRI == paste0("Y",n.yri) & RE.auc$N == paste0("N",N.trn) & RE.auc$anc == anc,]
    ij=which(sub.df[,6:10] == max(sub.df[,6:10]),arr.ind=T)
    RE.YRI=rbind.data.frame(RE.YRI,data.frame(anc,gwas.CEU = n.ceu,gwas.YRI = n.yri, N=N.trn,gamma=gamma[ij[,1]],lambda=lambda[ij[,2]],auc=round(sub.df[ij[,1],ij[,2]+5],3)))
  }
}
RE.YRI
