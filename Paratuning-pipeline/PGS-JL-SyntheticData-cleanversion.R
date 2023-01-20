rm(list=ls()); gc()
options(stringsAsFactors = F)

### This is the master pipeline for the analysis of simulated data.
### The pipeline uses for sets of genotype files for two distinct ancestries:
###  1) training sets for calculating the GWAS, which after the first step will be referred to as TRN.ANC
###  2) tuning sets which are used for selection hyper parameters (TUNE.ANC)
###  3) testing sets which are used to calculated genetic scores and AUC given the selected hyoer parameters (TST.ANC)
###  4) reference dataset of a set of random samples from a population (REF.ANC)

### these are simulations that are currently available
# 600-609 h2 = 0.80 in both populations
# 700-709 h2=0.80 in CEU and 0.60 in YRI
# 800-809 H2 approx 0.8, betas the same in both populations
if(!exists("i.sim")){
  i.sim <- 800
}

# load the general information
source("general_pipeline_parameters.R")
# load the pipeline functions
source("R-code-cleanversion/LS-pipeline-functions.R")
# Joint lassosum functions needed
source("R-code-cleanversion/jls-functions.R")  
Rcpp::sourceCpp(paste0("R-code-cleanversion/myfunctions.cpp"))

# create a directory to store the results
#results.dir=paste0(work.dir,"/",sim,"/Results")
#dir.create(results.dir,recursive = T,showWarnings = F)


#### perform joint lassomsum snp and beta selection
#TRN1="CEU"
#TRN2="YRI"
#REF1="CEU"
#REF2="YRI"
#sim="Sim-800"
#lambda=sort(0.021-exp(seq(log(0.001), log(0.02), length.out=20)))[4:8]
#gamma=seq(0.72,0.88,.04)
#shrink=.9
#chrs=1:22

load(paste0(main_simulation_pipeline_directory, "Work/Sim-",i.sim,"/simulation-params.RData"))
# main.dir <- params$run.info$main.dir #"/raid6/Tianyu/PRS/SimulationPipeline/"
work.dir <- params$run.info$work.dir #"/raid6/Tianyu/PRS/SimulationPipeline/Work/Sim-800/"
ParameterTuningDirectory <- paste0(work.dir, "/ParameterTuningData/")

# start the clock
ptm <- proc.time()
# for(sim in c(701:709,800:809)){
# sim=paste0("Sim-",i.sim)
# results.dir=paste0(work.dir,"/",sim,"/Results")
# dir.create(results.dir,recursive = T,showWarnings = F)

re<-jlsSelectSnp(TRN1="CEU", TRN2="YRI", 
                 REF1="CEU", REF2="YRI", 
                 # WORK.DIR = paste0(work.dir,"/",sim),
                 WORK.DIR = work.dir,
                 # REF.DIR = ref.dir,
                 REF.DIR = reference_genotype_directory,
                 ld.populations = ld.populations, 
                 # N = N,
                 N = sample_sizes,
                 # shrink = c(0.7,0.8,0.9),
                 shrink = c(0.9),
                 # lambda=sort(0.021-exp(seq(log(0.001), log(0.02), length.out=20)))[4:8], 
                 lambda = exp(seq(log(0.0025), log(0.025), length.out=10)),
                 # gamma=seq(0.76,0.86,.02),
                 gamma = c(0.2, 0.5, 0.8),
                 chrs=1:22)
# }

# stop the clock
proc.time() - ptm # total process takes now < 1 seconds


### perform the tuning on both datasets
# start the clock
ptm <- proc.time()
# for(sim in c(701:709,800:809)){
#   sim=paste0("Sim-",sim)
#   results.dir=paste0(work.dir,"/",sim,"/Results")
#   dir.create(results.dir,recursive = T,showWarnings = F)
#   
re.JLS.TUNING=NULL
for(TUNE in c("CEU","YRI")){
  re.jlsTuning<-jlsTuning(TRN1 = "CEU", TRN2 = "YRI", 
                          REF1= "CEU", REF2 = "YRI", 
                          TUNE = TUNE, 
                          WORK.DIR = paste0(work.dir,"/",sim), 
                          plink2 = plink2)
  re.JLS.TUNING=rbind(re.JLS.TUNING,re.jlsTuning)
}
# save the results from the tuning analysis
fwrite(re.JLS.TUNING,paste0(results.dir,"/jls-tuning-parameters.txt"),row.names=F,col.names=T,quote=F,sep="\t")
# }
# stop the clock
proc.time() - ptm # total process takes now < 1 seconds
# 
# 
# ### perform the testing step for the joint lassosum
# # start the clock
# ptm <- proc.time()
# for(sim in c(701:709,800:809)){
#   sim=paste0("Sim-",sim)
#   results.dir=paste0(work.dir,"/",sim,"/Results")
#   dir.create(results.dir,recursive = T,showWarnings = F)
#   
#   jls.TuningParameters=fread(paste0(results.dir,"/jls-tuning-parameters.txt"),header=T,data.table=F)
#   re.JLS.TESTING=NULL
#   for(TST in c("CEU","YRI")){
#     re.jlsTesting<-jlsTesting(TRN1="CEU",TRN2="YRI",REF1='CEU',REF2="YRI",TST,combn=jls.TuningParameters[jls.TuningParameters$TUNE == TST,"combn"],
#                               WORK.DIR=paste0(work.dir,"/",sim),plink2 = plink2)
#     re.JLS.TESTING=rbind.data.frame(re.JLS.TESTING,re.jlsTesting)
#   }
#   fwrite(re.JLS.TESTING,paste0(results.dir,"/jls-testing-auc.txt"),row.names=F,col.names=T,quote=F,sep="\t")
# }
# 
# # stop the clock
# proc.time() - ptm # total process takes now < 1 seconds
# 
# 
# 
# #### copy the reference data to the 
# # start the clock
# #ptm <- proc.time()
# #for(anc in c("CEU","YRI")){
# #  dir.create(paste0(ref.dir,"/",anc,".REF/CHR"),recursive = T,showWarnings = F)
# #  system.command=paste0("cp ",GNT.SETS[["REF"]][anc],"/CHR/*",
# #                        " ",paste0(ref.dir,"/",anc,".REF/CHR"))
# #  system(system.command)
# #}
# ## stop the clock
# #proc.time() - ptm # total process takes 34 seconds
# 
# # copy LD block information
# #ld.map=fread("/data3/Bert/PengLiu/SimulationCode-Feb2022/Data/chr1-22-qc-frq-ld.block.map",header=T,data.table=F)
# #ld.map=ld.map[,c("CHROM","ID","CEU.blk","YRI.blk")]
# ## write the ld.map by chromosome
# #dir.create(paste0(ref.dir,"/LD/CHR"),recursive = T,showWarnings = F)
# #for(chr in 1:22){
# #  fwrite(ld.map[ld.map$CHROM == paste0("chr",chr),],paste0(ref.dir,"/LD/CHR/ld.block-chr",chr,".map"),row.names=F,col.names=T,quote=F,sep="\t")
# #}
# #fwrite(ld.map,paste0(ref.dir,"/LD/ld.block",".map"),row.names=F,col.names=T,quote=F,sep="\t")
# 
# ##### start the analysis
# ######### START THE PIPELINE
# # select the simulation to use
# # select the betas and SNP to use for the lassosum
# # start the clock
# lambda=sort(0.021-exp(seq(log(0.001), log(0.02), length.out=20)))
# s=seq(0.8,1,0.05)
# 
# timestamp()
# ptm <- proc.time()
# for(sim in c(600:609,700:709,800:809)){
#   sim=paste("Sim",sim,sep="-")
#   
#   for(TRN in c("CEU","YRI")){
#     REF=TRN
#     re<-lsSelectSNP(TRN = TRN, REF = REF, WORK.DIR = paste0(work.dir,"/",sim), REF.DIR = ref.dir,ld.populations = ld.populations,
#                     N = N, s = s, lambda = lambda, chr = 1:22)
#   }
#   
# }
# # stop the clock
# proc.time() - ptm # total process takes 32 hours for all 4 * 30 combinations, ~ 30 minutes per combination
# timestamp()
# 
# # perform LassoSum Tuning
# timestamp()
# ptm <- proc.time()
# for(sim in c(600:609,700:709,800:809)){
#   sim=paste("Sim",sim,sep="-")
#   print(sim)
#   # create a directory to store the results
#   results.dir=paste0(work.dir,"/",sim,"/Results")
#   dir.create(results.dir,recursive = T,showWarnings = F)
#   
#   # store results
#   re.LS.TUNING=NULL
# 
#   for(TRN in c("CEU","YRI")){
#     REF=TRN
#     for(TUNE in c("CEU","YRI")){
#       re.lsTuning<-lsTuning(TRN = TRN, REF= REF, TUNE = TUNE, WORK.DIR = paste0(work.dir,"/",sim), plink2 = plink2)
#       re.LS.TUNING=rbind(re.LS.TUNING,re.lsTuning)
#     }
#   }
#   # save the results from the tuning analysis
#   fwrite(re.LS.TUNING,paste0(results.dir,"/ls-tuning-parameters.txt"),row.names=F,col.names=T,quote=F,sep="\t")
# }
# proc.time() - ptm # total process takes 5 minutes for 1 simulation. 1 minutes for each TUNE=CEU; 20 seconds for TUNE=YRI
# timestamp()
# 
# #### perform lassosum testing
# # start the clock
# ptm <- proc.time()
# # load the hyper parameters to use
# # create a directory to store the results
# for(sim in c(600:609,700:709,800:809)){
#   sim=paste0("Sim-",sim)
#   results.dir=paste0(work.dir,"/",sim,"/Results")
#   dir.create(results.dir,recursive = T,showWarnings = F)
#   ls.TuningParameters=fread(paste0(results.dir,"/ls-tuning-parameters.txt"),header=T,data.table=F)
#   re.LS.TESTING=NULL
#   for(TRN in c("CEU","YRI")){
#     REF=TRN
#     for(TST in c("CEU","YRI")){
#       re.lsTesting<-lsTesting(TRN,REF,TST,combn=ls.TuningParameters[ls.TuningParameters$TRN == TRN & ls.TuningParameters$REF == REF & ls.TuningParameters$TUNE == TST,"combn"],
#                                 WORK.DIR=paste0(work.dir,"/",sim),plink2 = plink2)
#       re.LS.TESTING=rbind.data.frame(re.LS.TESTING,re.lsTesting)
#     }
#   }
#   fwrite(re.LS.TESTING,paste0(results.dir,"/ls-testing-auc.txt"),row.names=F,col.names=T,quote=F,sep="\t")
# }
# # stop the clock
# proc.time() - ptm # total process takes 2 minutes
# 
# 
# #### this step will determine the optimal weights to combine the two LS PRS, I will refer to this as weighted LS
# #### I am using the selected cut-offs from step 3b
# # start the clock
# ptm <- proc.time()
# for(sim in c(700:709,800:809)){
#   sim=paste0("Sim-",sim)
#   results.dir=paste0(work.dir,"/",sim,"/Results")
#   dir.create(results.dir,recursive = T,showWarnings = F)
#   
#   ls.TuningParameters=fread(paste0(results.dir,"/ls-tuning-parameters.txt"),header=T,data.table=F)
#   re.WLS.TUNING=NULL
# #  for(REF in c("CEU","YRI")){
#     for(TUNE in c("CEU","YRI")){
#       re.wlsTuning=wlsTuning(TUNE=TUNE,WORK.DIR = paste0(work.dir,"/",sim), tuning.parameters = ls.TuningParameters)
#       re.WLS.TUNING=rbind.data.frame(re.WLS.TUNING,re.wlsTuning)
#     }
# #  }
#   fwrite(re.WLS.TUNING,paste0(results.dir,"/wls-tuning-parameters.txt"),row.names=F,col.names=T,quote=F,sep="\t")
# }
# # stop the clock
# proc.time() - ptm # total process takes 2 seconds
# 
# # apply the weights to the testing set
# ptm <- proc.time()
# for(sim in c(600:609,700:709,800:809)){
#   sim=paste0("Sim-",sim)
#   results.dir=paste0(work.dir,"/",sim,"/Results")
#   dir.create(results.dir,recursive = T,showWarnings = F)
#   
#   wls.TuningParameters=fread(paste0(results.dir,"/wls-tuning-parameters.txt"),header=T,data.table=F)
#   re.WLS.TESTING=NULL
# #  for(REF in c("CEU","YRI")){
#     for(TST in c("CEU","YRI")){
#       re.wlsTesting=wlsTesting(TST=TST, WORK.DIR = paste0(work.dir,"/",sim),ls.weights=wls.TuningParameters)
#       re.WLS.TESTING=rbind.data.frame(re.WLS.TESTING,re.wlsTesting)
#     }
# #  }
#   fwrite(re.WLS.TESTING,paste0(results.dir,"/wls-testing-auc.txt"),row.names=F,col.names=T,quote=F,sep="\t")
# }
# # stop the clock
# proc.time() - ptm # total process takes now < 1 seconds
# 
