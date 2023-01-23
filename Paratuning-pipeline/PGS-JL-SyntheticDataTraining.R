# rm(list=ls()); 
gc()
options(stringsAsFactors = F)

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


#### Functions to use
#Joint Lassosum Function Peng created
mylassosumFunction <- function(chr,gamma,lambda,shrink,
                               mem.limit,gwasANC,COR,
                               referenceFiles,LDblocks){
  re=mylassosum(cor1=COR[COR$CHR == chr,gwasANC[1]], 
                cor2=COR[COR$CHR == chr,gwasANC[2]],
                fileName1 = referenceFiles[[chr]][[gwasANC[1]]], 
                fileName2 = referenceFiles[[chr]][[gwasANC[2]]],  
                gamma = gamma, lambda = lambda, shrink=shrink,
                chunk=TRUE, mem.limit=mem.limit,
                trace=2,
                LDblocks1=as.data.frame(LDblocks[[gwasANC[1]]][LDblocks[[gwasANC[1]]]$chr == paste0("chr",chr),]), 
                LDblocks2=as.data.frame(LDblocks[[gwasANC[2]]][LDblocks[[gwasANC[2]]]$chr == paste0("chr",chr),]))
  return(re)
}

wrapperFunction <- function(i.combn, input.df, gwasANC, lambda, shrink, main.dir, work.dir, CHR=1:22, 
                            GWAS_file, output_directory,
                            mem.limit, mymc.cores = 8){
  
  # collect gamma
  gamma=input.df$gamma[i.combn]  
  
  # collectthe number of samples for the two populations
  gwasN=list()
  gwasN[["CEU"]]=input.df$N1[i.combn]
  gwasN[["YRI"]]=input.df$N2[i.combn]
  
  #### determine correlations
  ### Summary statistics ###
  # read a map for the first ancestry to set up a dataframe
  anc=gwasANC[1]
  map <- fread(paste0(work.dir, anc, '.TRN/', anc, '.TRN.pvar'), header=T, data.table=F)
  
  COR=data.frame(CHR=map$`#CHROM`,ID=map$ID)
  
  for(anc in gwasANC){
    # glm <- fread(paste0(work.dir, anc, '.TRN/Assoc/', anc, '.TRN.PHENO1.glm.logistic.hybrid'),
    #              header=T,data.table=F)
    #'/raid6/Ron/prs/data/bert_sample/CEU.TRN.PHENO1.glm.logistic.hybrid'
    glm <- fread(GWAS_file[[anc]], header = TRUE, data.table = FALSE)
    COR[,anc]=p2cor(p = glm$P, n = gwasN[[anc]], sign=log(glm$OR))
  }
  rownames(COR)=COR$ID
  
  # head(COR)
  # 
  # > head(COR)
  # CHR            ID          CEU          YRI
  # 1:1962845:T:C   1 1:1962845:T:C -0.008296110 0.0073444745
  # 1:1962899:A:C   1 1:1962899:A:C  0.013275142 0.0086428250
  # 1:1963406:G:A   1 1:1963406:G:A -0.006433539 0.0007753638
  # 1:1963538:T:C   1 1:1963538:T:C  0.004405373 0.0137580552
  # 1:1963738:C:T   1 1:1963738:C:T -0.006433539 0.0007753638
  # 1:1964101:A:G   1 1:1964101:A:G -0.006433539 0.0007753638
  
  # load the LD block information
  LDblocks=list()
  for(anc in gwasANC){
    if(anc == "CEU"){
      ld.anc="EUR"
    }else if(anc == "YRI"){
      ld.anc="AFR"
    }
    LDblocks[[anc]]=read.table2(system.file(paste0("data/Berisa.", 
                                                   paste0(ld.anc,".hg38"), ".bed"), 
                                            package="lassosum"), header=T)
    #  LDblocks[[anc]]$chr=gsub("chr","",LDblocks[[anc]]$chr)
  }
  
  # > head(LDblocks)
  # $CEU
  # chr     start      stop
  # 1     chr1   1961168   3666172
  # 2     chr1   3666172   4320751
  # 3     chr1   4320751   5853833
  # 4     chr1   5853833   7187275
  # 5     chr1   7187275   9305140
  # 6     chr1   9305140  10746927
  # 7     chr1  10746927  11717784
  # 8     chr1  11717784  12719464
  # 9     chr1  12719464  14565015
  
  
  # reference files by chromosome
  referenceFiles=list()
  for(chr in 1:22){
    referenceFiles[[chr]]=list()
    for(anc in gwasANC){
      # referenceFiles[[chr]][[anc]]=paste0("/data3/DownLoadedData/GWAS-Populations-SimulationInput/",anc,"_reference_LDblocks/CHR/",anc,"-chr",chr)
      referenceFiles[[chr]][[anc]] <- paste0(main.dir,"Data/Reference-LDblocks/", anc,"/CHR/", anc,"-chr",chr)
        # paste0(work.dir, "GWAS-Populations-SimulationInput/", anc, "_reference_LDblocks/CHR/", anc,"-chr",chr)
    }
  }
  
  # run lassosum
  print(CHR); flush.console()
  
  system.time(
    re.chr<-mclapply(CHR, mylassosumFunction,
                     gamma=gamma,lambda=lambda,
                     shrink=shrink,mem.limit = mem.limit,
                     gwasANC = gwasANC,COR=COR,
                     referenceFiles = referenceFiles,
                     LDblocks = LDblocks,
                     mc.cores = mymc.cores,
                     mc.preschedule = F)
  )
  print('loss in combined lassosum')
  print(re.chr$loss)
  print(re.chr$trainerror1)
  print(re.chr$trainerror2)
  
  # merge the results from the 22 chromosomes
  re.lasso = merge.mylassosum(re.chr)
  print('loss in combined lassosum, after merging')
  print(re.lasso$loss)
  print(re.lasso$trainerror1)
  print(re.lasso$trainerror2)
  # save ther results
  # dir.create(paste0(work.dir,"JointLassoSum/"),showWarnings = F,recursive = T)
  dir.create(output_directory, showWarnings = FALSE ,recursive = TRUE)
  # save(re.lasso,file=paste(work.dir,"JointLassoSum/JointLassosum-",sprintf("-gamma-%.2f",gamma),".Rdata",sep=""))
  save(re.lasso, file = paste0(output_directory, "JointLassosum-",
                               sprintf("-gamma-%.2f",gamma), ".Rdata"))
  return(i.combn)
}


#####

source("general_pipeline_parameters.R")
print('the directory of the main pipeline is')
print(main_simulation_pipeline_directory)

print('the directory of the parameter tuning pipeline is')
print(parameter_tuning_pipeline_directory)

####

# load the parameters for this simulation
load(paste0(main_simulation_pipeline_directory, "Work/Sim-",i.sim,"/simulation-params.RData"))
main.dir <- params$run.info$main.dir #"/raid6/Tianyu/PRS/SimulationPipeline/"
work.dir <- params$run.info$work.dir #"/raid6/Tianyu/PRS/SimulationPipeline/Work/Sim-800/"
ParameterTuningDirectory <- paste0(work.dir, "/ParameterTuningData/")
##

GWAS_file <- list(CEU = paste0(ParameterTuningDirectory, 
                               "Assoc/CEU-synthetic-train.PHENO1.glm.logistic.hybrid"),
                  YRI = paste0(ParameterTuningDirectory, 
                               "Assoc/YRI-synthetic-train.PHENO1.glm.logistic.hybrid"))
output_directory <- paste0(ParameterTuningDirectory, "JointLassosum/")

### ancestries
gwasANC <- c("CEU","YRI")
#### set the parameters

GAMMA <- c(0.2, 0.5, 0.8)
lambda <- exp(seq(log(0.0025), log(0.025), length.out=10))
shrink <- .9

##

chrs <- 1:22 #which chromosome did i use when training the model

#######

####sample sizes
N1 <- params$CEU.TRN$n.case + params$CEU.TRN$n.control
N2 <- params$YRI.TRN$n.case + params$YRI.TRN$n.control
input.df <- data.frame(gamma = GAMMA,
                       N1 = rep(N1, length(GAMMA)),
                       N2 = rep(N2, length(GAMMA)))
mem.limit=2*10e9

re.wrapper<-mclapply(1:nrow(input.df),wrapperFunction,
                     input.df=input.df,
                     gwasANC = gwasANC, 
                     lambda=lambda, 
                     shrink=shrink,
                     main.dir=main.dir,
                     work.dir=work.dir,
                     CHR = chrs, 
                     GWAS_file = GWAS_file,
                     output_directory = output_directory,
                     mem.limit=2e10,
                     mc.cores=3, mc.preschedule = F, mc.silent=F)



