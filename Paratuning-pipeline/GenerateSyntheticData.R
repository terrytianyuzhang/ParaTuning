###########generate boostrap samples and split them into training and validation######

# print('generate bootstrap data for cv')
# rm(list=ls()); 
gc()
options(stringsAsFactors = F)

if(!exists("plink") | !exists("plink2")){
  plink <- "/usr/local/bin/plink"
  plink2 <- "/usr/local/bin/plink2"
}

if(!exists("i.sim")){
  i.sim <- 800
}


# #### software needed
# plink <- "/usr/local/bin/plink"
# plink2 <- "/usr/local/bin/plink2"

# setwd("/home/tianyuz3/PRS/my_code")

library(data.table)
library(pryr) # check memory useage
library(lassosum) #transform p value to correlation
library(doParallel) # foreach
library(pROC) # for AUC 
library(R.utils)
library(snpStats)
library(wordspace)
library(foreach)
source("HyperparameterTuningFunctions.R")


# setting.title <- 'CEU0aYRI0a22Chr_lambda5'
# print(setting.title)
chrs <- 1:22 #which chromosome did i use when training the model

TrainTestNFold <- 5

gammaGenerateData <- 0.8
lambdaIndexGenerateData <- 5

# load the parameters for this simulation
load(paste0("/raid6/Tianyu/PRS/SimulationPipeline/Work/Sim-",i.sim,"/simulation-params.RData"))
main.dir <- params$run.info$main.dir #"/raid6/Tianyu/PRS/SimulationPipeline/"
work.dir <- params$run.info$work.dir #"/raid6/Tianyu/PRS/SimulationPipeline/Work/Sim-800/"

ParameterTuningDirectory <- paste0(work.dir, "/ParameterTuningData")
dir.create(ParameterTuningDirectory,
           showWarnings = F,recursive = T)
##########read in lasso results##########
# lasso.file <- paste0("/raid6/Tianyu/PRS/CombinedLassoSum/Tmp/GWAS-lasso-C20000-Y4000-gamma-0.50_", setting.title,".Rdata")
# lasso.file <- paste0("/raid6/Tianyu/PRS/CombinedLassoSum/Tmp/GWAS-lasso-C20000-Y4000-gamma-0.50_CEU1aYRI2a22Chr.Rdata")
TrainJLFile <- paste0(work.dir, 
                      'JointLassoSum/JointLassosum--gamma-', 
                      sprintf("%.2f",gammaGenerateData), 
                      '.Rdata')
TrainJLResult <- get(load(TrainJLFile))
AllBeta <- TrainJLResult$beta

# beta0.index <- 5 #this is for CEU1aYRI2a22Chr_lambda2
betaGenerateData <- AllBeta[, lambdaIndexGenerateData] #use a small lambda to generate bootstrap data
##########read in genotype and calculate score#############
# map <- fread('/raid6/Ron/prs/data/bert_sample/GWAS-Populations-SimulationInput/chr1-22-qc-frq-ld.block.map',header=T,data.table=F)[,c("CHROM","ID")]
map <- fread(paste0(main.dir,
                    'Data/chr1-22-qc-frq-ld.block.map'))
CHR <- gsub("chr","",map$CHROM)
####!!!!!!!!!!
# SNP <- map[CHR %in% 1:22, ]$ID
SNP <- map$ID
names(betaGenerateData) <- SNP

######SECTION 1: calculate PGS of the new individuals####
risk.score.list <- vector("list",2)
ancs <- c('CEU', 'YRI')
for(i.set in 1:2){
  anc <- ancs[i.set]
  re.pgss <- mclapply(chrs, PGS_bychr_bootstrap, anc = anc,
                      beta0 = betaGenerateData, mc.cores = 16, mc.preschedule = FALSE)
  ### sum the results
  pgs <- re.pgss[[1]] #this is the first chromosome
  for(i in 2:length(re.pgss)){
    pgs <- pgs+re.pgss[[i]]
  }

  risk.score.list[[i.set]] <- pgs

}
### save risk scores for each population

save(risk.score.list, file = paste0(ParameterTuningDirectory,
                                    'riskscore.Rdata'))
# ###pgs is the risk score for each subject in the reference panel
#

######SECTION 2: generate simulated Y####
risk.score.list <- get(load(paste0(ParameterTuningDirectory,
                                   'riskscore.Rdata')))
#####we need to read figure out what is the original data noise level
ancs <- c('CEU', 'YRI')
CEUSampleSize <- params$CEU.TRN$n.case + params$CEU.TRN$n.control
YRISampleSize <- params$YRI.TRN$n.case + params$YRI.TRN$n.control

# s.sizes <- c(20000, 4000)
s.sizes <- c(CEUSampleSize, YRISampleSize)
caseProportion <- params$CEU.TRN$n.case / CEUSampleSize ###the case proportion is the same for both populations


  
for(i.set in 1:2){
  TrainGWASFile <- paste0(work.dir, ancs[i.set], '.TRN/Assoc/',
                          ancs[i.set],'.TRN.PHENO1.glm.logistic.hybrid')
  SyntheticYFile <- paste0(ParameterTuningDirectory, '/', 
                           ancs[i.set],'-SyntheticY', '.RData')
  
  cv.generatey(TrainGWASFile = TrainGWASFile,
               SyntheticYFile = SyntheticYFile,
               anc = ancs[i.set], 
               s.size = s.sizes[i.set], 
               beta0 = betaGenerateData, 
               risk.score = risk.score.list[[i.set]], 
               case.prop = caseProportion)
}

print('generated booty')

######SECTION 3: split training and validation individuals####
###### 1/(nfold) left out for validation ####
##########split training and validation########
set.seed(2019)
# chrs <- 1:22
# ancs <- c('CEU', 'YRI')
# s.sizes <- c(20000, 4000)

for(i.set in 1:2){
  anc <- ancs[i.set]
  s.size <- s.sizes[i.set]
  
  # if(anc == 'CEU'){
  #   file.title <- 'CEU-20K'
  # }else{
  #   file.title <- 'YRI-4K'  
  # }
  # 
  # boost.data.folder <- paste0("/raid6/Tianyu/PRS/BootData/", file.title,"/CHR")
  # boost.data.train.index.R <- paste0(boost.data.folder, "/boost-train-index.RData")
  # boost.data.val.index.R <- paste0(boost.data.folder, "/boost-val-index.RData")
  TrainSampleIndexFile <- paste0(ParameterTuningDirectory, "/", anc, "-synthetic-train-index.txt")
  ValidationSampleIndexFile <- paste0(ParameterTuningDirectory, "/", anc, "-synthetic-validate-index.txt")
  TrainSampleIndexRFile <- paste0(ParameterTuningDirectory, "/", anc, "-synthetic-train-index.Rdata")
  ValidationSampleIndexRFile <- paste0(ParameterTuningDirectory, "/", anc, "-synthetic-validate-index.Rdata")
  
  # boost.data.train.index.R <- paste0(ParameterTuningDirectory, "/", anc, "-synthetic-train-index.RData")
  # boost.data.val.index.R <- paste0(ParameterTuningDirectory, "/", anc, "-synthetic-validate-index.RData")
  # 
  val.index <- sort(sample(1:s.size, floor(s.size/TrainTestNFold)))
  train.index <- (1:s.size)[-val.index] #this is in order
  
  # save(train.index, file = paste0("/raid6/Tianyu/PRS/BootData/",anc,".TUNE/",setting.title,"train_index.RData"))
  # save(val.index, file = paste0("/raid6/Tianyu/PRS/BootData/",anc,".TUNE/", setting.title,"val_index.RData"))
  save(val.index, file = ValidationSampleIndexRFile)
  save(train.index, file = TrainSampleIndexRFile)
  
  
  ####generate training and testing .fam files
  mclapply(chrs, splitTrainValidation, anc = anc,
          train.index = train.index, val.index = val.index,
          ParameterTuningDirectory = ParameterTuningDirectory,
          TrainSampleIndexFile = TrainSampleIndexFile,
          ValidationSampleIndexFile = ValidationSampleIndexFile,
          plink = plink, mc.cores = 22, mc.preschedule = FALSE)
}

print('finished sample splitting')

######SECTION 4: calculate pairwise gene/simulated Y association####
#########created simulated summary statistics########
######
# chrs <- 1:22
# ancs <- c('CEU', 'YRI')
# chrs <- 1:6
for(i.set in 1:2){
  anc <- ancs[i.set]
  
  if(anc == 'CEU'){
    file.title <- 'CEU-20K'
  }else{
    file.title <- 'YRI-4K'  
  }
  
  # boost.data.folder <- paste0("/raid6/Tianyu/PRS/BootData/", file.title,"/CHR")
  # boost.data.train.index.R <- paste0(boost.data.folder, "/boost-train-index.RData")
  TrainSampleIndexRFile <- paste0(ParameterTuningDirectory, "/", anc, "-synthetic-train-index.Rdata")
  train.index <- get(load(TrainSampleIndexRFile))
  
  SyntheticYFile <- paste0(ParameterTuningDirectory, '/', 
                           anc,'-SyntheticY', '.RData')
  SyntheticY <- get(load(SyntheticYFile))
  SyntheticY <- SyntheticY[train.index] #only keep the training samples' Y
  
  SyntheticY <- matrix(SyntheticY, ncol = 1)
  SyntheticY <- SyntheticY - rep(1, nrow(SyntheticY)) %*% t(colMeans(SyntheticY))
  SyntheticY <- normalize.cols(SyntheticY, method="euclidean",p=2)
  
  SyntheticGWAS <- mclapply(chrs, PairwiseCorrelationSyntheticData,
                            ParameterTuningDirectory = ParameterTuningDirectory,
                        anc = anc,
                        SyntheticY = SyntheticY,
                        mc.cores = 3, mc.preschedule = FALSE, mc.silent=F)
  # chr = 22
  # ParameterTuningDirectory = ParameterTuningDirectory
  # anc = anc
  # SyntheticY = SyntheticY
  # cl <- parallel::makeCluster(22, type = 'FORK')
  # 
  # system.time(
  #   SyntheticGWAS <- parallel::parLapplyLB(cl, chrs, PairwiseCorrelationSyntheticData,
  #                         ParameterTuningDirectory = ParameterTuningDirectory,
  #                         anc = anc,
  #                         SyntheticY = SyntheticY))
  # ##   user  system elapsed
  # ##  0.004   0.009   3.511
  # 
  # parallel::stopCluster(cl)


  # cl<-makeCluster(16, type = 'FORK')
  # registerDoParallel(cl)
  # 
  # SyntheticGWASForeach <- foreach(chr = chrs)  %dopar%  {
  #           PairwiseCorrelationSyntheticData(chr = chr,
  #           ParameterTuningDirectory = ParameterTuningDirectory,
  #           anc = anc,
  #           SyntheticY = SyntheticY)
  #           }
  # stopCluster(cl)
  
  for(ChromosomeIndex in 1:length(chrs)){
    if(ChromosomeIndex == 1){
      SyntheticGWASCombined <- SyntheticGWAS[[ChromosomeIndex]]
    }else{
      SyntheticGWASCombined <- rbind(SyntheticGWASCombined, SyntheticGWAS[[ChromosomeIndex]])
    }
  }
  # fwrite(one_fake_GWAS, file = paste0('/raid6/Tianyu/PRS/BootData/',anc,'_bootdata_repeat',i))
  write.table(SyntheticGWASCombined, 
         file = paste0(ParameterTuningDirectory,'/', anc,'-SyntheticGWAS'))
  
  # booty <- get(load(paste0("/raid6/Tianyu/PRS/BootData/",anc,"_bootY_",setting.title,".RData")))
  # train.index <- get(load(boost.data.train.index.R))
  # booty <- booty[train.index] #only use the training samples' Y
  
  # boot.Y <- matrix(booty, ncol = 1)
  # # center and normalize the boostrap Y
  # boot.Y <- boot.Y - rep(1, nrow(boot.Y)) %*% t(colMeans(boot.Y))
  # 
  # boot.Y <- normalize.cols(boot.Y, method="euclidean",p=2)
  
  # fake_GWAS <- mclapply(chrs, cor_bychr_bootstrap, anc = anc, boot.Y = boot.Y, 
  #                       B = B, mc.cores = 4, mc.preschedule = FALSE, mc.silent=F)
  # ##the length of fake_GWAS = number of chromosome
  # ##each component is a list, whose length = B, number of repeats
  # 
  # ###correctly merge the data and save it
  # for(i in 1:B){
  #   for(j in 1:length(chrs)){
  #     if(j == 1){
  #       one_fake_GWAS <- fake_GWAS[[j]][[i]]
  #     }else{
  #       one_fake_GWAS <- rbind(one_fake_GWAS, fake_GWAS[[j]][[i]])
  #     }
  #   }
  #   # fwrite(one_fake_GWAS, file = paste0('/raid6/Tianyu/PRS/BootData/',anc,'_bootdata_repeat',i))
  #   fwrite(one_fake_GWAS, file = paste0('/raid6/Tianyu/PRS/BootData/',anc,'_bootdata_',setting.title,'_repeat',i))
  #   
  # }
}  

print('pairwise association calculation')


