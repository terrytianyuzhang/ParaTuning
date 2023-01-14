rm(list=ls()); gc()
options(stringsAsFactors = F)

### This is the master pipeline for the analysis of simulated data.
### The pipeline uses three sets of genotype files for two distinct ancestries:
###  1) training sets for calculating the GWAS, which after the first step will be referred to as TRN.ANC
###  2) tuning sets which are used for selection hyper parameters (TUNE.ANC)
###  3) testing sets which are used to calculated genetic scores and AUC given the selected hyoer parameters (TST.ANC)

# genetic scores will be calculated using thre transfer learning approach. It does this in threee steps
# 1) determine the different sets of betas using a TRaiNing lassomsum results and a TARGET GWAS and ld structure
# 2) find the set of betas that give the maximum AUC in a TUNEing set of the TARGET ancestry
# 3) apply this optimal set to a TeSTing dataset of the TARGET ancestry


### these are simulations that are currently available
# 600-609 h2 = 0.80 in both populations
# 700-709 h2=0.80 in CEU and 0.60 in YRI
# 800-809 H2 approx 0.8, betas the same in both populations


# load the general information
source("general-pipeline-information.R")

# load the pipeline functions
source("R-code/TL-pipeline-functions.R")

# some additional code  specific to transfer learning
source("R-code/TransferLearningFunctions/my_BasicFunctions.R")
source("R-code/TransferLearningFunctions/my_TL_PRS.R")

# create the sets of beta to be used

for(sim in c(601:609,701:709,800:809)){
  sim=paste0("Sim-",sim)
  TRN="CEU"
  TARGET="YRI"
  system.time(re<-tlBetas(TRN=TRN, TARGET=TARGET, WORK.DIR=paste0(work.dir,"/",sim),REF.DIR=ref.dir,ld.populations=ld.populations,plink=plink))
}

# find the best set of betas using the tuneing data from the TARGET ancestry
for(sim in c(600:609,700:709,800:809)){
  sim=paste0("Sim-",sim)
  results.dir=paste0(work.dir,"/",sim,"/Results")
  dir.create(results.dir,recursive = T,showWarnings = F)
  TRN="CEU"
  TARGET="YRI"
  
  system.time(re.TL.TUNING<-tlTuning(TRN, TARGET, paste0(work.dir,"/",sim), plink2))  
  
  fwrite(re.TL.TUNING,paste0(results.dir,"/tl-tuning-parameters.txt"),row.names=F,col.names=T,quote=F,sep="\t")
  
}

### perform the testing step for the transfer learning
for(sim in c(600:609,700:709,800:809)){
  sim=paste0("Sim-",sim)
  results.dir=paste0(work.dir,"/",sim,"/Results")
  dir.create(results.dir,recursive = T,showWarnings = F)
  
  tl.TuningParameters=fread(paste0(results.dir,"/tl-tuning-parameters.txt"),header=T,data.table=F)
  re.TL.TESTING=NULL
  TRN="CEU"
  TST="YRI"
  re.TL.TESTING<-tlTesting(TRN=TRN,TST=TST,combn=tl.TuningParameters[tl.TuningParameters$TARGET == TST,"combn"],
                              WORK.DIR=paste0(work.dir,"/",sim),plink2 = plink2)
  fwrite(re.TL.TESTING,paste0(results.dir,"/tl-testing-auc.txt"),row.names=F,col.names=T,quote=F,sep="\t")
}


