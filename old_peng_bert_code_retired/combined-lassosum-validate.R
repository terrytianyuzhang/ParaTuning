rm(list=ls()); gc()
options(stringsAsFactors = F)


#### software needed
plink="/data3/Software/Plink/plink"
plink2="/data3/Software/Plink2/plink2"

#### load the functions that are needed
source("R-code/simulation-functions.R")

### Load functions ###
source(paste0("R-code/LassosumFunctions/parseselect.R"))
source(paste0("R-code/LassosumFunctions/parseblocks.R"))
source(paste0("R-code/LassosumFunctions/ncol.bfile.R"))
source(paste0("R-code/LassosumFunctions/nrow.bfile.R"))
source(paste0("R-code/LassosumFunctions/read.table2.R"))
source(paste0("R-code/LassosumFunctions/selectregion.R"))
source(paste0("R-code/LassosumFunctions/parse.pheno.covar.R"))
source(paste0("R-code/LassosumFunctions/myelnet.R"))
source(paste0("R-code/LassosumFunctions/mylassosum.R"))
source(paste0("R-code/LassosumFunctions/splitgenome.R"))
source(paste0("R-code/LassosumFunctions/validation.R"))
source(paste0("R-code/LassosumFunctions/merge.mylassosum.R"))
Rcpp::sourceCpp(paste0("R-code/LassosumFunctions/myfunctions.cpp"))

library(data.table)
library(pryr) # check memory useage
library(lassosum) #transform p value to correlation
library(doParallel) # foreach
library(pROC) # for AUC 
library(R.utils)
library(snpStats)
library(wordspace)

#### load the parameter information
i.sim="NC8000"
load(paste0("Work/Sim-C0.80-Y0.60-rep",i.sim,"/simulation-params.RData"))

# set directories for this simulation
main.dir=params$run.info$main.dir
work.dir=params$run.info$work.dir


### load one of the lassosum runs
#load("Work/Sim-C0.80-Y0.60-rep1/CombinedLassoSum/Tmp/GWAS-lasso-C12000-Y12000-gamma-0.20.Rdata")

## next step is validation this needs to be done using the validate command, due to the transformation performed on the genotypes within lasso
## this transformation involves setting the mean to zero and the makingeach genotype vector have unit length (2-norm, sum of squares is 1)
## the approach will be to do this by chromosome, this means that each of the training/validation genotypes needs to be split by chromosome
## ughhh.



#system.time(re.validate <- validate(re.lasso, "Work/Sim-C0.80-Y0.60-rep1/CEU.TRN/CEU.TRN-4000", extract=NULL, keep=NULL, distance=NULL, mem.limit=mem.limit))

#save(re.validate,file="Tmp/re.validate-CEU-4000.RData")




pgsCalculationByCHR <- function(chr,work.dir,trn.set,trn.n,shrink,CHR,beta){
  # find the beta of interest
  snp=rownames(beta)[CHR == chr & rowSums(beta != 0) > 0]
  system.time(gnt<-read.plink(bed=paste0(work.dir,trn.set,"/",trn.set,"-",trn.n,"-chr",chr,".bed"),
                              bim=paste0(work.dir,trn.set,"/",trn.set,"-",trn.n,"-chr",chr,".bim"),
                              fam=paste0(work.dir,trn.set,"/",trn.set,"-",trn.n,"-chr",chr,".fam"),select.snps = snp))
  # transform the genotypes to a matrix, and reverse the call count. This snpStats package counts the A2 allele, not the A1
  system.time(gnt<-2-as(gnt$genotypes,Class="numeric"))
  # center the columns
  system.time(gnt <-gnt - rep(1, nrow(gnt)) %*% t(colMeans(gnt)))
  # normalize the calls to the unit 1 norm
  system.time(gnt<-normalize.cols(gnt,method="euclidean",p=2))
  # apply the proper shrinkage
  system.time(gnt<-gnt*sqrt(1-shrink))
  # calculate the pgs
  system.time(re.pgs<-gnt%*%beta[snp,])
  
  # return the results
  return(re.pgs)
}




### some type of function  i.set,work.df,work.dir,SNP,CHR
pgsCalculation <- function(i.set,work.df=work.df,work.dir=work.dir,SNP=SNP,CHR=CHR){
  # set parameters
  lasso.file=work.df$lasso.file[i.set]
  trn.set=work.df$trn.set[i.set]
  trn.n=work.df$trn.n[i.set]

  # load the lasso file and collect important information
  load(lasso.file)
  beta=re.lasso$beta
  rownames(beta)=SNP
  colnames(beta)=sprintf("lambda:%.4f",re.lasso$lambda)
  shrink=re.lasso$shrink
  rm(re.lasso)

  # run the validation
  system.time(re<-mclapply(1:22,pgsCalculationByCHR,work.dir=work.dir,trn.set=trn.set,trn.n=trn.n,CHR=CHR,shrink=shrink,beta=beta,mc.cores=10,mc.preschedule = F))
  
  ### sum the results
  pgs=re[[1]]
  for(i in 2:length(re)){
    pgs=pgs+re[[i]]
  }

  # write the results
  gwasset=unlist(strsplit(lasso.file,"/"))
  gwasset=gsub(".Rdata","",gsub("GWAS","",gwasset[length(gwasset)]))
  fwrite(as.data.frame(pgs),paste0(work.dir,"CombinedLassoSum/PGS-",trn.set,"-",trn.n,gwasset,".pgs"),row.names=T,col.names=T,quote=F,sep="\t")

  return(i.set)
}






map=fread("/data3/DownLoadedData/GWAS-Populations-SimulationInput/chr1-22-qc-frq-ld.block.map",header=T,data.table=F)[,c("CHROM","ID")]
CHR=gsub("chr","",map$CHROM)
SNP=map$ID
rm(map)

lasso.files=list.files(path=paste0(work.dir,"CombinedLassoSum/Tmp/"),pattern = ".Rdata",full.names = T)
work.df=data.frame(trn.set=rep(c("CEU.TRN","YRI.TRN"),length(lasso.files)),trn.n=20000,lasso.file=rep(lasso.files,each=2))

system.time(re.out<-mclapply(1:nrow(work.df),pgsCalculation,work.df=work.df,work.dir=work.dir,SNP=SNP,CHR=CHR,mc.cores=2,mc.preschedule = F))

