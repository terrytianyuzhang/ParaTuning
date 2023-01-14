#rm(list=ls()); gc()
#options(stringsAsFactors = F)

library(data.table)
library(doParallel)
library(snpStats)
library(R.utils)

#### software needed
plink="/data3/Software/Plink/plink"
plink2="/data3/Software/Plink2/plink2"


#### load the functions that are needed
source("simulation-functions.R")

#### load the parameters for the simulation
load(paste0("Work/Sim-",i.sim,"/simulation-params.RData"))
main.dir=params$run.info$main.dir
work.dir=params$run.info$work.dir

# Determine the PRS for each of the training set
train.sets=names(params)[grep(".TRN",names(params))]
tune.sets=names(params)[grep(".TUNE",names(params))]

for(train.set in train.sets){
  for(tune.set in tune.sets){
    train.dir=paste0(work.dir,train.set,"/")
    tune.dir=paste0(work.dir,tune.set,"/")
    # create an output file
    dir.create(paste0(tune.dir,"PRS/",train.set,"/"),showWarnings = F,recursive = T)
    
    # testing genotype file
    tune.pfile=paste0(tune.dir,tune.set)
    # training GWAS data files
    gwas.file=paste(train.dir,"Assoc/",train.set,".PHENO1.glm.logistic.hybrid",sep="")
    
    ### read the training gwas results 
    gwas=fread(gwas.file,header=T,data.table=F)
    # determine the regression coefficient from the odds ratios
    gwas$B=log(gwas$OR)
    # write a short version with the essential information
    fwrite(gwas[,c(3,6,15)],paste0(tune.dir,"PRS/",train.set,"/",train.set,".SNP.values"),row.names=F,col.names=T,quote=F,sep="\t")
    
    
    # calculate PRS for each r2 cut-off and then use a range of p-values cut-offs
    for(r2 in c(0.1,0.2,0.5,0.8)){
      # read the clumping file
      clump.file=paste0(train.dir,"Clumping/",train.set,sprintf("-r2_%2.1f",r2),".clumped")
      clumps=fread(clump.file,header=T,data.table=F)
      # write a file with the inde SNP id for each of the clumps
      write.table(as.data.frame(clumps$SNP),paste0(tune.dir,"PRS/",train.set,"/",train.set,sprintf("-r2_%2.1f",r2),".clumped.valid.snp"),row.names=F,col.names=F,quote=F,sep="\t")
      
      # use plink to calcualte the PRS score
      plink.command=paste(plink2,"--allow-no-sex","--nonfounders",
                          "--pfile",tune.pfile,
                          "--score",paste0(tune.dir,"PRS/",train.set,"/",train.set,".SNP.values"),"header",
                          "--q-score-range","Data/pval-range-list.txt",gwas.file,3,13,"header",
                          "--extract",paste0(tune.dir,"PRS/",train.set,"/",train.set,sprintf("-r2_%2.1f",r2),".clumped.valid.snp"),
                          "--out",paste0(tune.dir,"PRS/",train.set,"/",train.set,"-",tune.set,sprintf("-r2_%2.1f",r2),"-p"),
                          sep=" ")
      system(plink.command)
    }
  }
}

# retain information that you want to carry over to another part of the simulation
rm.list=ls()
rm.list=rm.list[!(rm.list %in% c("i.sim","sims"))]
rm(list = rm.list); flush.console()

