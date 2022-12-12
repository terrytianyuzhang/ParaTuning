#rm(list=ls()); 
gc()
options(stringsAsFactors = F)

library(data.table)
library(doParallel)
library(snpStats)
library(R.utils)
library(lassosum)
library(pROC)

#### software needed
# plink="/data3/Software/Plink/plink"
# plink2="/data3/Software/Plink2/plink2"
if(!exists("plink") | !exists("plink2")){
  plink <- "/usr/local/bin/plink"
  plink2 <- "/usr/local/bin/plink2"
}

if(!exists("i.sim")){
  i.sim <- 800
}

#### load the functions that are needed
source("simulation-functions.R")

# load the parameters for this simulation
load(paste0("Work/Sim-",i.sim,"/simulation-params.RData"))
main.dir=params$run.info$main.dir
work.dir=params$run.info$work.dir

train.sets=names(params)[grep(".TRN",names(params))] #"CEU.TRN" "YRI.TRN"
tune.sets=names(params)[grep(".TUNE",names(params))] #"CEU.TUNE" "YRI.TUNE"
test.sets=names(params)[grep(".TST",names(params))] #"CEU.TST" "YRI.TST"


#Define the ld populations
ld.sets=c("CEU","YRI"); names(ld.sets)=train.sets
# #> ld.sets
# CEU.TRN YRI.TRN 
# "CEU"   "YRI" 
ld.pops=c("EUR.hg38","AFR.hg38"); names(ld.pops)=train.sets


for(train.set in train.sets){
  # select the datasets and directories to use
  train.dir=paste0(work.dir,train.set,"/") #"/raid6/Tianyu/PRS/SimulationPipeline/Work/Sim-800/CEU.TRN/"
  for(tune.set in tune.sets){
    test.set=gsub("TUNE","TST",tune.set) #"CEU.TST"
    
    # information for reference genotypes
    ld.set=ld.sets[train.set]
    # CEU.TRN 
    # "CEU" 
    ld.pop=ld.pops[train.set]
    
    # data to work with
    tune.dir=paste0(work.dir,tune.set,"/")
    test.dir=paste0(work.dir,test.set,"/")
  
    # set parameteres for the lassosum
    lambda=seq(0.0001,0.1, length.out=20)
    s=c(0.2,0.5,0.9,1.0)

    # read the pvalues for this ancestry training set
    gwas=fread(paste0(work.dir,train.set,"/","Assoc/",train.set,".PHENO1.glm.logistic.hybrid"),header=T,data.table=F)
  
    # transform the pvalues to correlations with the phenotype, this is called "r" in the algorithm
    r <- p2cor(p = gwas$P, n = params[[train.set]]$n.case+params[[train.set]]$n.control, sign=log(gwas$OR))

    timestamp()
    re.lasso=NULL
    for(chr in 1:22){
      print(chr); flush.console()
      ld.bfile=paste0(main.dir,"Data/Reference-LDblocks/",ld.set,"/CHR/",ld.set,"-chr",chr)
      tune.bfile=paste0(work.dir,tune.set,"/CHR/",tune.set,"-chr",chr)
      
      timestamp()
      cl <- makeCluster(24, type="FORK") # Parallel over 10 nodes
      re.chr <- lassosum.pipeline(cor=r[gwas$`#CHROM`== chr], chr=gwas$`#CHROM`[gwas$`#CHROM`== chr], pos=gwas$POS[gwas$`#CHROM`== chr], 
                                         A1=gwas$A1[gwas$`#CHROM`== chr], A2=gwas$REF[gwas$`#CHROM`== chr], # A2 is not required but advised
                                         ref.bfile=ld.bfile, test.bfile=tune.bfile,s=s,lambda=lambda,
                                         cluster = cl,
                                         LDblocks = ld.pop)
      
      stopCluster(cl)  
      gc()
      timestamp()
      # merge the results
      if(is.null(re.lasso)){
        re.lasso=re.chr
      }else{
        re.lasso=merge(re.lasso,re.chr)
      }
    }
    timestamp()
    
  # save the results
  dir.create(paste0(tune.dir,"/Lassosum"),showWarnings = F, recursive = T)
  save(re.lasso,file=paste0(tune.dir,"Lassosum/",train.set,"-lassosum-",tune.set,".RData"))
  #"/raid6/Tianyu/PRS/SimulationPipeline/Work/Sim-800/CEU.TUNE/Lassosum/CEU.TRN-lassosum-CEU.TUNE.RData"

  # now analyze the results
  re.validate.lasso=validate(re.lasso)
  
  #### find the combination of s and lambda that gives the best AUC
  df.AUC=NULL
  for(cs in as.character(s)){
    for(i.lambda in 1:length(lambda)){
      auc=roc(re.validate.lasso$results.table$pheno~re.validate.lasso$pgs[[cs]][,i.lambda])$auc[1]
      df.AUC=rbind.data.frame(df.AUC,data.frame(s=as.numeric(cs),lambda=lambda[i.lambda],auc=auc))
    }
  }
  # save(df.AUC, file = paste0(tune.dir,"Lassosum/",tune.set,"/",train.set,"-",tune.set,"-lassosum.AUC.RData"))
  
  re.subset.lasso <- subset(re.lasso, s=df.AUC[which.max(df.AUC$auc),"s"], lambda=df.AUC[which.max(df.AUC$auc),"lambda"])
  
  # tuning set
  re.tune=re.validate.lasso$results.table
  i.s=which(re.validate.lasso$s == df.AUC[which.max(df.AUC$auc),"s"]) 
  i.lambda=which(re.validate.lasso$lambda == df.AUC[which.max(df.AUC$auc),"lambda"]) 
  re.tune[,"best.pgs"]=re.validate.lasso$pgs[[i.s]][,i.lambda]
  roc(pheno~best.pgs,data=re.tune)
  
  dir.create(paste0(tune.dir,"Lassosum/",tune.set),showWarnings=F,recursive=T)
  fwrite(re.tune,paste0(tune.dir,"Lassosum/",tune.set,"/",train.set,"-",tune.set,"-lassosum.score"),row.names=F,col.names=T,quote=F,sep="\t")
  #"/raid6/Tianyu/PRS/SimulationPipeline/Work/Sim-800/CEU.TUNE/Lassosum/CEU.TUNE/CEU.TRN-CEU.TUNE-lassosum.score"
  
  # testing set
  cl <- makeCluster(24, type="FORK") # Parallel over 10 nodes
  re.test <- validate(re.subset.lasso, test.bfile=paste0(work.dir,test.set,"/",test.set),cluster=cl)
  stopCluster(cl)  
  gc()
  roc(pheno~best.pgs,data=re.test$results.table)
  
  dir.create(paste0(test.dir,"Lassosum/",test.set),showWarnings=F,recursive=T)
  fwrite(re.test$results.table,paste0(test.dir,"Lassosum/",test.set,"/",train.set,"-",test.set,"-lassosum.score"),row.names=F,col.names=T,quote=F,sep="\t")
  #"/raid6/Tianyu/PRS/SimulationPipeline/Work/Sim-800/CEU.TST/Lassosum/CEU.TST/CEU.TRN-CEU.TST-lassosum.score"
  }
}

print('finished PGS-lassosum')
rm.list=ls()
rm.list=rm.list[!(rm.list %in% c("i.sim","sims"))]
rm(list = rm.list); flush.console()



  
  
  



