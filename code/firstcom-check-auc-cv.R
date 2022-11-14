####load in the combined-lasso beta, and get testing auc
set.seed(2019)
rm(list=ls()); gc()
options(stringsAsFactors = F)

print('check AUC for CV data')
setting.title <- 'CEU0aYRI0a22Chr_lambda5'
setwd("/home/tianyuz3/PRS/my_code")

source("hyperpara-tuning-function.R")
library(data.table)
library(pryr) # check memory useage
library(lassosum) #transform p value to correlation
library(doParallel) # foreach
library(pROC) # for AUC 
library(R.utils)
library(snpStats)
library(wordspace)


main.dir <- "/raid6/Tianyu/PRS/"
work.dir <- "/raid6/Ron/prs/data/bert_sample/"

lasso.file <- paste0("/raid6/Tianyu/PRS/CombinedLassoSum/Tmp/GWAS-lasso-C20000-Y20000-gamma-0.50_boot_",setting.title, "_1.Rdata")

####now we load the beta0
load(lasso.file)
beta <- re.lasso$beta
shrink <- re.lasso$shrink 
lambda <- re.lasso$lambda

####
# map=fread("/data3/DownLoadedData/GWAS-Populations-SimulationInput/chr1-22-qc-frq-ld.block.map",header=T,data.table=F)[,c("CHROM","ID")]
map <- fread('/raid6/Ron/prs/data/bert_sample/GWAS-Populations-SimulationInput/chr1-22-qc-frq-ld.block.map',header=T,data.table=F)[,c("CHROM","ID")]

# > head(map)
# CHROM            ID
# 1  chr1 1:1962845:T:C
# 2  chr1 1:1962899:A:C
# 3  chr1 1:1963406:G:A
# 4  chr1 1:1963538:T:C
# 5  chr1 1:1963738:C:T
# 6  chr1 1:1964101:A:G

####!!!!!!
CHR=gsub("chr","",map$CHROM)
# > head(CHR)
# [1] "1" "1" "1" "1" "1" "1"
####!!!!!!!!!!
SNP <- map$ID
# SNP <- map[CHR %in% 1:22, ]$ID

chrs <- 1:22
ancs <- c('CEU', 'YRI')
rownames(beta) <- SNP

for(i.set in 1:2){

  anc <- ancs[i.set]
  
  booty <- get(load(file = paste0('/raid6/Tianyu/PRS/BootData/',anc,'_bootY_',setting.title,'.RData')))
  
  val.index <- get(load(paste0("/raid6/Tianyu/PRS/BootData/",anc,".TUNE/",setting.title, "val_index.RData")))
  
  pheno <- booty[val.index] #simulated true Y
  
  PGSnPHENO <- matrix(0, ncol = length(lambda) + 1, nrow = length(pheno))
  PGSnPHENO[,length(lambda) + 1 ] <- pheno
  
  re.pgss <- mclapply(chrs, PGS_bychr_bootstrap_allbeta, anc = anc, 
                      beta = beta, mc.cores = 8, mc.preschedule = F)
  
  ### sum the results
  pgs <- re.pgss[[1]] #this is the first chromosome
  for(i in 2:length(re.pgss)){
    pgs <- pgs+re.pgss[[i]]
  }
  PGSnPHENO[,1:length(lambda)] <- pgs
  ###pgs is the risk score for each subject in the reference panel
  save(PGSnPHENO, file = paste0('/raid6/Tianyu/PRS/AUC/',setting.title,'PGSnPHENO_',anc))
}



#######load the PGS score and phenotype information, then calculate ROC
library(data.table)
library(pROC)
setting.title <- 'CEU0aYRI0a22Chr_lambda5'
anc <- 'CEU'
PGSnPHENO <- get(load(paste0('/raid6/Tianyu/PRS/AUC/',setting.title,'PGSnPHENO_',anc)))

nlambda <- NCOL(PGSnPHENO) - 1
PGSnPHENO[PGSnPHENO[,nlambda + 1] >0, nlambda + 1] <- 0.5
PGSnPHENO[PGSnPHENO[,nlambda + 1] <0, nlambda + 1] <- -0.5
aucs <- rep(0, nlambda)
for(i in 1:nlambda){
  aucs[i] <- auc(PGSnPHENO[, nlambda + 1], PGSnPHENO[,i])
}
save(aucs, file = paste0('/raid6/Tianyu/PRS/AUC/',setting.title,'PGSnPHENO_auc_',anc, '.RData'))

anc <- 'YRI'
PGSnPHENO <- get(load(paste0('/raid6/Tianyu/PRS/AUC/',setting.title,'PGSnPHENO_',anc)))

nlambda <- NCOL(PGSnPHENO) - 1
PGSnPHENO[PGSnPHENO[,nlambda + 1] >0, nlambda + 1] <- 0.5
PGSnPHENO[PGSnPHENO[,nlambda + 1] <0, nlambda + 1] <- -0.5
aucs <- rep(0, nlambda)
for(i in 1:nlambda){
  aucs[i] <- auc(PGSnPHENO[, nlambda + 1], PGSnPHENO[,i])
}
save(aucs, file = paste0('/raid6/Tianyu/PRS/AUC/',setting.title,'PGSnPHENO_auc_',anc, '.RData'))




  
# }