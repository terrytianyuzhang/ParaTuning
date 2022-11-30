#rm(list=ls()); gc()
options(stringsAsFactors = F)

library(data.table)
library(doParallel)
library(snpStats)
library(R.utils)

#### software needed
# plink="/data3/Software/Plink/plink"
# plink2="/data3/Software/Plink2/plink2"
# directory changes
if(!exists("plink") | !exists("plink2")){
  plink <- "/usr/local/bin/plink"
  plink2 <- "/usr/local/bin/plink2"
}


#### load the functions that are needed
source("simulation-functions.R")

#### load the parameters for the simulation
load(paste0("Work/Sim-",i.sim,"/simulation-params.RData"))
main.dir=params$run.info$main.dir
work.dir=params$run.info$work.dir

train.sets=names(params)[grep(".TRN",names(params))]
for(set in train.sets){
  set.dir=paste0(work.dir,set,"/")
  
  # translate the pgen format file to a bed format file
  plink2.command=paste(plink2,"--nonfounders",
                       "--pfile",paste(set.dir,set,sep=""),
                       "--make-bed",
                       "--out",paste(set.dir,set,sep=""),
                       sep=" ")
  system(plink2.command)
  
  ### create a directory for the clumping
  dir.create(paste0(set.dir,"Clumping/"),recursive = T,showWarnings = F)
  
  #### split the GWAS results by chromosome, this avoids all kinds of warning message
  gwas=fread(paste(set.dir,"Assoc/",set,".PHENO1.glm.logistic.hybrid",sep=""),header=T,data.table=F)
  chrs=1:22
  for(chr in chrs){
    fwrite(gwas[gwas$`#CHROM` == chr,],paste(set.dir,"Assoc/",set,"-chr",chr,".PHENO1.glm.logistic.hybrid",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
  }
  
  
  r2.list=c(0.1,0.2,0.5,0.8)
  for(r2 in r2.list){
    print(r2); flush.console()
    
    mclapply(chrs,clumpingFunction,r2=r2,set=set,set.dir=set.dir,plink=plink,mc.cores=8,mc.preschedule=F,mc.silent=T)
    
    
    clumps=NULL
    for(chr in chrs){
      clumps=rbind.data.frame(clumps,fread(paste0(set.dir,"Clumping/",set,"-chr",chr,sprintf("-r2_%2.1f",r2),".clumped"),header=T,data.table=F))
    }
    clumps=clumps[order(clumps$P),]
    fwrite(clumps,paste0(set.dir,"Clumping/",set,sprintf("-r2_%2.1f",r2),".clumped"),row.names=F,col.names=T,quote=F,sep="\t")
  
  }
  
  # remove the intermediate files
  system.command=paste("rm",paste0(set.dir,"Clumping/","*chr*"))
  system(system.command)
  
  
}

# retain information that you want to carry over to another part of the simulation
rm.list=ls()
rm.list=rm.list[!(rm.list %in% c("i.sim","sims"))]
rm(list = rm.list); flush.console()

