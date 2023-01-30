####load in the combined-lasso beta, and get testing auc
gc()
options(stringsAsFactors = F)

if(!exists("plink") | !exists("plink2")){
  plink <- "/usr/local/bin/plink"
  plink2 <- "/usr/local/bin/plink2"
}

if(!exists("i.sim")){
  i.sim <- 800
}

library(data.table)
library(pryr) # check memory useage
library(lassosum) #transform p value to correlation
library(doParallel) # foreach
library(pROC) # for AUC 
library(R.utils)
library(snpStats)
library(wordspace)

PGS_bychr_bootstrap <- function(chr, anc, beta, work.dir, genotype_prefix_list){
  print(paste0('claculating PGS of ',anc, ' using chr ', chr))
  ######next step is generating some risk score using reference genotype
  ######we can download these genotype from 1000 Genome Project
  ######for the purpose of thee paper I will just use the reference panel 
  ######from which the training samples are generated
  
  chr_loc <- as.numeric(gsub(':.*$','',rownames(beta)))
  snp <- rownames(beta)[chr_loc == chr]
  
  ####load genotype data without label
  genotype_prefix <- genotype_prefix_list[[chr]][[anc]]
  system.time(gnt<-read.plink(bed=paste0(genotype_prefix, ".bed"),
                              bim=paste0(genotype_prefix, ".bim"),
                              fam=paste0(genotype_prefix, ".fam"),
                              select.snps = snp)
  )
  #example : "/raid6/Tianyu/PRS/bert_sample/YRI.TUNE/CHR/YRI.TUNE-chr20.bed"
  
  # transform the genotypes to a matrix, and reverse the call count. This snpStats package counts the A2 allele, not the A1
  system.time(gnt<- 2 - as(gnt$genotypes,Class="numeric"))
  
  # center the columns
  system.time(gnt <-gnt - rep(1, nrow(gnt)) %*% t(colMeans(gnt)))
  # normalize the calls to the unit 1 norm
  system.time(gnt<-normalize.cols(gnt,method="euclidean",p=2))
  
  # apply the proper shrinkage
  # system.time(gnt<-gnt*sqrt(1-shrink))
  
  # calculate the pgs
  system.time(re.pgs<-gnt%*%beta[chr_loc == chr,])
  print(paste0('finished claculating PGS of ',anc, ' using chr ', chr))
  
  return(re.pgs)
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
ParameterTuningDirectory <- paste0(work.dir, "/ParameterTuningData",
                                   "_gamma_", sprintf("%.2f",gammaGenerateData), 
                                   "_lambda_", sprintf("%.4f",lambda[lambdaIndexGenerateData]), "/")
lassosum_directory <- paste0(ParameterTuningDirectory, "JointLassosum/")

####CREATE A LIST INDICATING GENOTYPE FILE NAME

validation_genotype_file_list <- list()
for(chr_index in 1:22){
  
  one_chr_genotype_file <- list()
  
  ####
  
  for(ancestry in c("CEU", "YRI")){
    one_chr_genotype_file[[ancestry]] <- paste0(ParameterTuningDirectory, "CHR/", ancestry,"-chr", chr_index,"synthetic-val")
  }
  
  ####
  
  validation_genotype_file_list[[chr_index]] <- one_chr_genotype_file
}

####FINISH:CREATE A LIST INDICATING GENOTYPE FILE NAME

GAMMA = c(0.2, 0.5, 0.8)
chrs <- 1:22


####
# map=fread("/data3/DownLoadedData/GWAS-Populations-SimulationInput/chr1-22-qc-frq-ld.block.map",header=T,data.table=F)[,c("CHROM","ID")]
map <- fread(paste0(main.dir, 'Data/chr1-22-qc-frq-ld.block.map'),header=T,data.table=F)[,c("CHROM","ID")]

for(gamma in GAMMA){
  lasso.file <- paste0(lassosum_directory, "JointLassosum-",
                       sprintf("-gamma-%.2f",gamma), ".Rdata")
  
  ####now we load the beta0
  load(lasso.file)
  beta <- re.lasso$beta
  shrink <- re.lasso$shrink 
  lambda <- re.lasso$lambda
  
  ####
  
  CHR=gsub("chr","",map$CHROM)
  SNP <- map[CHR %in% 1:22, ]$ID
  
  ####
  for(i.set in 1:2){
    
    ancestry <- if(i.set == 1){'CEU'} else 'YRI'
    
    ####read in the testing phenotype
    
    ValidationSampleIndexRFile <- paste0(ParameterTuningDirectory, 
                                    ancestry, "-synthetic-validate-index.Rdata")
    val.index <- get(load(ValidationSampleIndexRFile))
    
    ####
    
    SyntheticYFile <- paste0(ParameterTuningDirectory, '/',
                             ancestry,'-SyntheticY', '.RData')
    SyntheticY <- get(load(SyntheticYFile))
    SyntheticY <- SyntheticY[val.index]
    
    ####
    
    PGSnPHENO <- matrix(0, ncol = length(lambda) + 1, nrow = length(SyntheticY))
    PGSnPHENO[,length(lambda) + 1 ] <- SyntheticY
    
    ####
    rownames(beta) <- SNP ##should this be here?
    
    
    re.pgss <- mclapply(chrs, PGS_bychr_bootstrap, anc = ancestry, 
                        beta = beta, work.dir = work.dir, 
                        genotype_prefix_list = validation_genotype_file_list,
                        mc.cores = 4, mc.preschedule = F)
    
    ### sum the results
    pgs <- re.pgss[[1]] #this is the first chromosome
    for(i in 2:length(re.pgss)){
      pgs <- pgs+re.pgss[[i]]
    }
    PGSnPHENO[,1:length(lambda)] <- pgs
    
    PGS_file <- paste0(lassosum_directory, "JointLassosum-",sprintf("-gamma-%.2f",gamma),"-", ancestry, "-synthetic-validation-PGS.Rdata")
    save(PGSnPHENO, file = PGS_file)
  }
  
  
  
  #######load the PGS score and phenotype information, then calculate ROC
  for(i.set in 1:2){
    ancestry <- if(i.set == 1){'CEU'} else 'YRI'
    print(ancestry)
    ####
    
    PGS_file <- paste0(lassosum_directory, "JointLassosum-",sprintf("-gamma-%.2f",gamma),"-", ancestry, "-synthetic-validation-PGS.Rdata")
    PGSnPHENO <- get(load(PGS_file))
    
    ####
    
    nlambda <- NCOL(PGSnPHENO) - 1
    aucs <- rep(0, nlambda)
    ###AUC for each lambda
    for(i in 1:nlambda){
      aucs[i] <- auc(PGSnPHENO[, nlambda + 1], PGSnPHENO[,i])
    }
    save(aucs, file = paste0(lassosum_directory, "/JointLassosum-",sprintf("-gamma-%.2f",gamma),"-", ancestry, "-synthetic-validation-AUC.Rdata"))
  }
}





