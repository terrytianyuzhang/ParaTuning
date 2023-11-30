###########generate boostrap samples and split them into training and validation######

gc()
options(stringsAsFactors = F)

if(!exists("i.sim")){
  i.sim <- 600
}

library(data.table)
library(pryr) # check memory useage
library(lassosum) #transform p value to correlation
library(doParallel) # foreach
library(pROC) # for AUC 
library(R.utils)
library(snpStats)
library(wordspace)
library(foreach)
#####
source("HyperparameterTuningFunctions.R")
#####
source("general_pipeline_parameters.R")
print('the directory of the main pipeline is')
print(main_simulation_pipeline_directory)

# print('the directory of the parameter tuning pipeline is')
# print(parameter_tuning_pipeline_directory)

# load the parameters for this simulation
# load(paste0(main_simulation_pipeline_directory, "Work/Sim-",i.sim,"/simulation-params.RData"))
# main.dir <- params$run.info$main.dir #"/raid6/Tianyu/PRS/SimulationPipeline/"
# work.dir <- params$run.info$work.dir #"/raid6/Tianyu/PRS/SimulationPipeline/Work/Sim-800/"
TimingResultDirectory <- '/raid6/Tianyu/PRS/parameter_tuning_from_GWAS'
timing_result_slow <- system.time({
ParameterTuningDirectory <- paste0(work.dir, "ParameterTuningData",
                                   "_gamma_", sprintf("%.2f",gammaGenerateData), 
                                   "_lambda_", sprintf("%.4f",lambda[lambdaIndexGenerateData]))
dir.create(ParameterTuningDirectory,
           showWarnings = F,recursive = T)
##########read in lasso results##########
TrainJLFile <- paste0(work.dir, 
                      'JointLassoSum/JointLassosum--gamma-', 
                      sprintf("%.2f",gammaGenerateData), 
                      '.Rdata')
TrainJLResult <- get(load(TrainJLFile))
AllBeta <- TrainJLResult$beta

betaGenerateData <- AllBeta[, lambdaIndexGenerateData] #use a small lambda to generate bootstrap data
##########read in genotype and calculate score#############
map <- fread(paste0(main.dir,
                    'Data/chr1-22-qc-frq-ld.block.map'))
CHR <- gsub("chr","",map$CHROM)
SNP <- map$ID
names(betaGenerateData) <- SNP


# ######SECTION 1: calculate PGS of the new individuals####
# risk.score.list <- vector("list",2)
####save the selected preliminary beta
beta_generate_data_file <- './temp_file_generate_data/beta_generate_data.txt'
A1 <- unlist(lapply(strsplit(SNP, ":"),`[[`,4))
effect_size_df <- data.frame(SNP = SNP, 
                             A1 = A1,
                             BETA = betaGenerateData)
write.table(effect_size_df, beta_generate_data_file, sep = "\t", quote = FALSE, row.names = FALSE)

ancs <- c('CEU', 'YRI')
for(i.set in 1:2){
  anc <- ancs[i.set]
  # re.pgss <- mclapply(chrs, PGS_by_chr, anc = anc,
  #                     beta0 = betaGenerateData, mc.cores = 8, mc.preschedule = FALSE)
  # 
  # 
  # ### sum the results
  # pgs <- re.pgss[[1]] #this is the first chromosome
  # for(i in 2:length(re.pgss)){
  #   pgs <- pgs+re.pgss[[i]]
  # }
  # 
  # risk.score.list[[i.set]] <- pgs
  
  ##calculate the PRS for each individual in the reference panel
  if(anc == 'CEU'){
    file.title <- 'CEU-20K'
  }else{
    file.title <- 'YRI-4K'  
  }
  
  reference_genotype_file <- paste0("/raid6/Tianyu/PRS/bert_sample/ReferencePopulation-Package/", file.title, "/", file.title)
  PRS_file <- "./temp_file_generate_data/PRS_generate_data"
  plink2.command = paste(plink2,"--nonfounders","--allow-no-sex","--threads",24,"--memory",25000,
                         "--bfile", reference_genotype_file ,
                         "--score", beta_generate_data_file,"header-read",1,2,
                         "--score-col-nums",3,
                         "--out",PRS_file,
                         sep=" ")
  
  system(plink2.command)
}
### save risk scores for each population

# save(risk.score.list, file = paste0(ParameterTuningDirectory,
                                    # '/riskscore.Rdata'))
})
timing_file_slow <- paste0(TimingResultDirectory,
                      '/timing_result_slow.rds')
saveRDS(timing_result_slow, file = timing_file_slow)
# ###pgs is the risk score for each subject in the reference panel
#
timing_result <- system.time({

######SECTION 2: generate simulated Y####
risk.score.list <- get(load(paste0(ParameterTuningDirectory,
                                   '/riskscore.Rdata')))
#####we need to read figure out what is the original data noise level
ancs <- c('CEU', 'YRI')
CEUSampleSize <- sample_sizes$CEU$n.case + sample_sizes$CEU$n.control
YRISampleSize <- sample_sizes$YRI$n.case + sample_sizes$YRI$n.control

s.sizes <- c(CEUSampleSize, YRISampleSize)
caseProportion <- sample_sizes$CEU$n.case / CEUSampleSize ###the case proportion is the same for both populations

for(i.set in 1:2){
  TrainGWASFile <- paste0(work.dir, 'TRN/',
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

print('generated synthetic outcome Y')


######SECTION 3: split training and validation individuals####
###### 1/(nfold) left out for validation ####
##########split training and validation########
set.seed(2019)

for(i.set in 1:2){
  anc <- ancs[i.set]
  s.size <- s.sizes[i.set]
  
  if(anc == 'CEU'){
    file.title <- 'CEU-20K'
  }else{
    file.title <- 'YRI-4K'
  }
  
  TrainSampleIndexFile <- paste0(ParameterTuningDirectory, "/", anc, "-synthetic-train-index.txt")
  ValidationSampleIndexFile <- paste0(ParameterTuningDirectory, "/", anc, "-synthetic-validate-index.txt")
  TrainSampleIndexRFile <- paste0(ParameterTuningDirectory, "/", anc, "-synthetic-train-index.Rdata")
  ValidationSampleIndexRFile <- paste0(ParameterTuningDirectory, "/", anc, "-synthetic-validate-index.Rdata")
  
  val.index <- sort(sample(1:s.size, floor(s.size/TrainTestNFold)))
  train.index <- (1:s.size)[-val.index] #this is in order
  
  save(val.index, file = ValidationSampleIndexRFile)
  save(train.index, file = TrainSampleIndexRFile)
  
  ####
  referece_panel_allSNP_name <- paste0("/raid6/Tianyu/PRS/bert_sample/ReferencePopulation-Package/", file.title,"/", file.title)
  all_sample_psam <- fread(paste0(referece_panel_allSNP_name, ".psam"))
  train_psam <- all_sample_psam[train.index, c(1,2)] #only keep family id and withtin family id
  fwrite(train_psam, TrainSampleIndexFile, col.names = F, sep = " ")
  val_psam <- all_sample_psam[val.index, c(1,2)]
  fwrite(val_psam, ValidationSampleIndexFile, col.names = F, sep = " ")
  
  ####generate training and testing .fam files
  mclapply(chrs, splitTrainValidation, anc = anc,
          train.index = train.index, val.index = val.index,
          ParameterTuningDirectory = ParameterTuningDirectory,
          TrainSampleIndexFile = TrainSampleIndexFile,
          ValidationSampleIndexFile = ValidationSampleIndexFile,
          plink = plink, mc.cores = 16, mc.preschedule = FALSE)
  print(paste(anc, 'is done splitting'))
}

print('finished sample splitting')

######SECTION 4: calculate GWAS on synthetic data######
chrs <- 1:22
ancs <- c('CEU', 'YRI')

####

dir.create(paste0(ParameterTuningDirectory, '/Assoc/'),
           showWarnings = F,recursive = T)

####

for(population_index in 1:2){
  ####
  ancestry <- ancs[population_index]
  
  if(ancestry == 'CEU'){
    file.title <- 'CEU-20K'
  }else{
    file.title <- 'YRI-4K'
  }
  ####
  
  TrainSampleIndexRFile <- paste0(ParameterTuningDirectory, "/", 
                                  ancestry, "-synthetic-train-index.Rdata")
  train.index <- get(load(TrainSampleIndexRFile))
  
  ####
  
  SyntheticYFile <- paste0(ParameterTuningDirectory, '/',
                           ancestry,'-SyntheticY', '.RData')
  SyntheticY <- get(load(SyntheticYFile))
  SyntheticY <- SyntheticY[train.index]
  
  ####
  
  for(chr in chrs){
    
    unlabeled_fam_file <- paste0(ParameterTuningDirectory,
                                 "/CHR/", ancestry, 
                                 "-chr", chr, "synthetic-train.fam")
    unlabeled_fam <- fread(unlabeled_fam_file)
    
    unlabeled_fam$V6 <- SyntheticY + 1
    fwrite(unlabeled_fam, file = unlabeled_fam_file, sep = ' ', col.names = FALSE)
    
    
    ####
    
    
    plink2.command=paste(plink2,"--nonfounders","--allow-no-sex",
                         "--bfile", paste0(ParameterTuningDirectory, "/CHR/", 
                                           ancestry,"-chr", chr,"synthetic-train"),
                         "--glm","allow-no-covars","omit-ref",
                         "--out", paste0(ParameterTuningDirectory, "/Assoc/", 
                                         ancestry,"-chr", chr,"synthetic-train"),
                         sep=" ")
    system(plink2.command)
    
  }
  ####
  
  all_GWAS <- data.table()
  for(chr in chrs){
    one_chr_GWAS <- fread(paste0(ParameterTuningDirectory, 
                               "/Assoc/", ancestry,"-chr", chr,
                               "synthetic-train.PHENO1.glm.logistic.hybrid"),
                        header = T)
    all_GWAS <- rbind(all_GWAS, one_chr_GWAS)
  }
  
  ####
  fwrite(all_GWAS, paste0(ParameterTuningDirectory, 
                          "/Assoc/", ancestry,"-",
                          "synthetic-train.PHENO1.glm.logistic.hybrid"),
         sep = ' ')
  
}
})
timing_file <- paste0(TimingResultDirectory,
                      '/timing_result.rds')
saveRDS(timing_result, file = timing_file)
