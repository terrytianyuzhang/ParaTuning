rm(list=ls()); gc()
options(stringsAsFactors = F)

### This is the master pipeline for the analysis of simulated data.
### The pipeline uses three sets of genotype files for two distinct ancestries:
###  1) training sets for calculating the GWAS, which after the first step will be referred to as TRN.ANC
###  2) tuning sets which are used for selection hyper parameters (TUNE.ANC)
###  3) testing sets which are used to calculated genetic scores and AUC given the selected hyoer parameters (TST.ANC)

# genetic scores will be calculated using 
# 1) pruning and tresholding using a single ancestry GWAS
# 2) Weighted P&T using a weighted value from the PT in 1) for the two ancestries


### these are simulations that are currently available
# 600-609 h2 = 0.80 in both populations
# 700-709 h2=0.80 in CEU and 0.60 in YRI
# 800-809 H2 approx 0.8, betas the same in both populations


# load the general information
source("general-pipeline-information.R")

# load the pipeline functions
source("R-code/PnT-pipeline-functions.R")

######### START THE PIPELINE
# select the simulation to use
for(sim in c(600:609,700:709,800:809)){
sim=paste("Sim",sim,sep="-")

#### copy the data files to Data in the current directory
# start the clock
ptm <- proc.time()
from.dir=paste0("/data3/Bert/PengLiu/SimulationCode-Feb2022/Work/",sim)
to.dir=paste0(data.dir,"/",sim)
dir.create(to.dir,recursive = T,showWarnings = F)
for(anc in c("CEU","YRI")){
  for(role in c("TRN","TUNE","TST")){
    system.command=paste0("cp"," ",from.dir,"/",anc,".",role,"/",anc,".",role,".p*"," ",to.dir)
    system(system.command)
  }
}
# stop the clock
proc.time() - ptm # total process takes 32 seconds

# create a directory to store the results
results.dir=paste0(work.dir,"/",sim,"/Results")
dir.create(results.dir,recursive = T,showWarnings = F)

# step 1 - prepare the genotype datasets so that they can be used across the pipeline
# this involves creating a plink bed version of the complete data as well as creating separate files for each of the chromosomes.
# start the clock
ptm <- proc.time()
for(role in c("TRN","TUNE","TST")){
  for(anc in c("CEU","YRI")){
    prepareGenotypes(GNT.SET=gsub("XXX",sim,GNT.SETS[[role]][[anc]]),ANC = anc, ROLE = role,WORK.DIR = paste0(work.dir,"/",sim,"/",anc,".",role), plink2 = plink2)
  }
}
# stop the clock
proc.time() - ptm # total process takes ~ 2 minutes

# step 2 - calculate the GWAS using the TRN data. This should be available when analyzing real data.
# start the clock
ptm <- proc.time()
for(anc in c("CEU","YRI")){
  runGWAS(ANC = anc,WORK.DIR=paste0(work.dir,"/",sim,"/",anc,".","TRN"), plink2 = plink2)
}
# stop the clock
proc.time() - ptm # total process takes 1 minute

# step 3a - pruning step for P&T - this is the time consuming step (for CEU - CEU - 0.8 it takes ~ 3 minutes)
# needed are GWAS and tuning genotype datasets for each ancestry
# start the clock
ptm <- proc.time()
for(TRN in c("CEU","YRI")){
  for(TUNE in c("CEU","YRI")){
    for(r2 in c(0.2,0.5,0.8)){
      system.time(re.PNT.PRUNING <-pntPruning(TRN = TRN,TUNE = TUNE,WORK.DIR = paste0(work.dir,"/",sim),r2 = r2,plink = plink))
    }
  }
}
# stop the clock
proc.time() - ptm # total process takes 20 minutes

# step 3b - find the optimal p-value cut-off for P&T using he TUNEing set genotypes and the pruning information derived from the TRaiNing set
# and TUNEing set genotypes.  
# start the clock
ptm <- proc.time()
re.PNT.TUNING=NULL
for(TRN in c("CEU","YRI")){
  for(TUNE in c("CEU","YRI")){
    re.pntTuning=NULL
    for(r2 in c(0.2,0.5,0.8)){
      re=pntTuning(TRN=TRN,TUNE=TUNE,WORK.DIR=paste0(work.dir,"/",sim),pval.range.file=paste0(data.dir,"/pval-range-list.txt"),r2=r2,plink2=plink2)
      re.pntTuning=rbind.data.frame(re.pntTuning,re)
    }
    # save the collected results
    fwrite(re.pntTuning,paste0(paste0(work.dir,"/",sim),"/",TUNE,".TUNE/PNT-TUNING/",TRN,".TRN/",TRN,"-PNT-TUNING.txt"),row.names=F,col.names=T,quote=F,sep="\t")
    # select the one with the largest AUC
    re.PNT.TUNING=rbind(re.PNT.TUNING,re.pntTuning[which.max(re.pntTuning$auc),])
  }
}
fwrite(re.PNT.TUNING,paste0(results.dir,"/pnt-tuning-parameters.txt"),row.names=F,col.names=T,quote=F,sep="\t")
# stop the clock
proc.time() - ptm # total process takes now 2 minutes (was 16 minutes)

#### step 3c now the optimal hyper-parameters are known, apply them to the testing set
# start the clock
ptm <- proc.time()
# load the hyper parameters to use
pnt.TuningParameters=fread(paste0(results.dir,"/pnt-tuning-parameters.txt"),header=T,data.table=F)
re.PNT.TESTING=NULL
for(TRN in c("CEU","YRI")){
  for(TST in c("CEU","YRI")){
    re.pntTesting<-pntTesting(TRN,TST,r2=pnt.TuningParameters[pnt.TuningParameters$TRN == TRN & pnt.TuningParameters$TUNE == TST,"r2"],
               pvalue=pnt.TuningParameters[pnt.TuningParameters$TRN == TRN & pnt.TuningParameters$TUNE == TST,"pvalue"],WORK.DIR=paste0(work.dir,"/",sim),plink2 = plink2)
    re.PNT.TESTING=rbind.data.frame(re.PNT.TESTING,re.pntTesting)
  }
}
fwrite(re.PNT.TESTING,paste0(results.dir,"/pnt-testing-auc.txt"),row.names=F,col.names=T,quote=F,sep="\t")
# stop the clock
proc.time() - ptm # total process takes now 13 seconds

#### step 4a will determine the optimal weights to combine the two P&T PRS, I will refere to this as weighted P&T
#### I am using the selected cut-offs from step 3b
# start the clock
ptm <- proc.time()
for(sim in c(600:609,700:709,800:809)){
  sim=paste("Sim",sim,sep="-")
  # create a directory to store the results
  results.dir=paste0(work.dir,"/",sim,"/Results")
  dir.create(results.dir,recursive = T,showWarnings = F)
  
  pnt.TuningParameters=fread(paste0(results.dir,"/pnt-tuning-parameters.txt"),header=T,data.table=F)
  re.WPNT.TUNING=NULL
  for(TUNE in c("CEU","YRI")){
    re.wpntTuning=wpntTuning(TUNE=TUNE,WORK.DIR = paste0(work.dir,"/",sim), tuning.parameters = pnt.TuningParameters)
    re.WPNT.TUNING=rbind.data.frame(re.WPNT.TUNING,re.wpntTuning)
  }
  fwrite(re.WPNT.TUNING,paste0(results.dir,"/wpnt-tuning-parameters.txt"),row.names=F,col.names=T,quote=F,sep="\t")
  # stop the clock
}
proc.time() - ptm # total process takes 2 seconds
##### step 4b perform the weighted PRS on the TeSTing datasets
# start the clock
ptm <- proc.time()
for(sim in c(600:609,700:709,800:809)){
  sim=paste("Sim",sim,sep="-")
  print(sim); flush.console()
  # create a directory to store the results
  results.dir=paste0(work.dir,"/",sim,"/Results")
  dir.create(results.dir,recursive = T,showWarnings = F)
  
  wpnt.TuningParameters=fread(paste0(results.dir,"/wpnt-tuning-parameters.txt"),header=T,data.table=F)
  re.WPNT.TESTING=NULL
  for(TST in c("CEU","YRI")){
    re.wpntTesting=wpntTesting(TST, WORK.DIR = paste0(work.dir,"/",sim),pnt.weights=wpnt.TuningParameters)
    re.WPNT.TESTING=rbind.data.frame(re.WPNT.TESTING,re.wpntTesting)
  }
  fwrite(re.WPNT.TESTING,paste0(results.dir,"/wpnt-testing-auc.txt"),row.names=F,col.names=T,quote=F,sep="\t")
}
# stop the clock
proc.time() - ptm # total process takes now < 1 seconds

}