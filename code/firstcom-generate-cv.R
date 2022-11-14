###########generate boostrap samples and split them into training and validation######
set.seed(2019)
print('generate bootstrap data for cv')
rm(list=ls()); gc()
options(stringsAsFactors = F)

#### software needed
plink <- "/usr/local/bin/plink"
plink2 <- "/usr/local/bin/plink2"

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

setting.title <- 'CEU0aYRI0a22Chr_lambda5'
print(setting.title)
chrs <- 1:22 #which chromosome did i use when training the model

nfold <- 5
##########read in lasso results##########
# lasso.file <- paste0("/raid6/Tianyu/PRS/CombinedLassoSum/Tmp/GWAS-lasso-C20000-Y4000-gamma-0.50_", setting.title,".Rdata")
lasso.file <- paste0("/raid6/Tianyu/PRS/CombinedLassoSum/Tmp/GWAS-lasso-C20000-Y4000-gamma-0.50_CEU1aYRI2a22Chr.Rdata")

re.lasso <- get(load(lasso.file))
beta <- re.lasso$beta

beta0.index <- 5 #this is for CEU1aYRI2a22Chr_lambda2
beta0 <- beta[, beta0.index] #use a small lambda to generate bootstrap data
##########read in genotype and calculate score#############
map <- fread('/raid6/Ron/prs/data/bert_sample/GWAS-Populations-SimulationInput/chr1-22-qc-frq-ld.block.map',header=T,data.table=F)[,c("CHROM","ID")]

CHR <- gsub("chr","",map$CHROM)
####!!!!!!!!!!
# SNP <- map[CHR %in% 1:22, ]$ID
SNP <- map$ID
names(beta0) <- SNP

######SECTION 1: calculate PGS of the new individuals####
risk.score.list <- vector("list",2)
ancs <- c('CEU', 'YRI')
for(i.set in 1:2){
  anc <- ancs[i.set]
  re.pgss <- mclapply(chrs, PGS_bychr_bootstrap, anc = anc,
                      beta0 = beta0, mc.cores = 8, mc.preschedule = FALSE)
  ### sum the results
  pgs <- re.pgss[[1]] #this is the first chromosome
  for(i in 2:length(re.pgss)){
    pgs <- pgs+re.pgss[[i]]
  }

  risk.score.list[[i.set]] <- pgs

}
### save risk scores for each population
save(risk.score.list, file = paste0('/raid6/Tianyu/PRS/RiskScore/', setting.title,'riskscore.Rdata'))
# ###pgs is the risk score for each subject in the reference panel
#

######SECTION 2: generate simulated Y####
risk.score.list <- get(load(paste0('/raid6/Tianyu/PRS/RiskScore/', setting.title,'riskscore.Rdata')))
#####we need to read figure out what is the original data noise level
ancs <- c('CEU', 'YRI')
s.sizes <- c(20000, 4000)
for(i.set in 1:2){
  cv.generatey(anc = ancs[i.set], 
               s.size = s.sizes[i.set], 
               beta0 = beta0, 
               risk.score = risk.score.list[[i.set]], 
               case.prop = 0.5)
}

print('generated booty')

######SECTION 3: split training and validation individuals####
###### 1/(nfold) left out for validation ####
##########split training and validation########
set.seed(2019)
chrs <- 1:22
ancs <- c('CEU', 'YRI')
s.sizes <- c(20000, 4000)

for(i.set in 1:2){
  anc <- ancs[i.set]
  s.size <- s.sizes[i.set]

  val.index <- sort(sample(1:s.size, floor(s.size/nfold)))
  save(val.index, file = paste0("/raid6/Tianyu/PRS/BootData/",anc,".TUNE/", setting.title,"val_index.RData"))

  train.index <- (1:s.size)[-val.index] #this is in order
  save(train.index, file = paste0("/raid6/Tianyu/PRS/BootData/",anc,".TUNE/",setting.title,"train_index.RData"))

  ####generate training and testing .fam files
  mclapply(chrs, split.train.val, anc = anc,
          train.index = train.index, val.index = val.index,
          plink = plink, mc.cores = 22, mc.preschedule = FALSE)
}

print('finished sample splitting')

######SECTION 4: calculate pairwise gene/simulated Y association####
#########created simulated summary statistics########
######
B <- 1
chrs <- 1:22
ancs <- c('CEU', 'YRI')

for(i.set in 1:2){
  anc <- ancs[i.set]
  
  booty <- get(load(paste0("/raid6/Tianyu/PRS/BootData/",anc,"_bootY_",setting.title,".RData")))
  train.index <- get(load(paste0("/raid6/Tianyu/PRS/BootData/",anc,".TUNE/",setting.title, "train_index.RData")))
  booty <- booty[train.index] #only use the training samples' Y
  
  boot.Y <- matrix(booty, ncol = 1)
  # center and normalize the boostrap Y
  boot.Y <- boot.Y - rep(1, nrow(boot.Y)) %*% t(colMeans(boot.Y))
  
  boot.Y <- normalize.cols(boot.Y, method="euclidean",p=2)
  
  fake_GWAS <- mclapply(chrs, cor_bychr_bootstrap, anc = anc, boot.Y = boot.Y, 
                        B = B, mc.cores = 4, mc.preschedule = FALSE, mc.silent=F)
  ##the length of fake_GWAS = number of chromosome
  ##each component is a list, whose length = B, number of repeats
  
  ###correctly merge the data and save it
  for(i in 1:B){
    for(j in 1:length(chrs)){
      if(j == 1){
        one_fake_GWAS <- fake_GWAS[[j]][[i]]
      }else{
        one_fake_GWAS <- rbind(one_fake_GWAS, fake_GWAS[[j]][[i]])
      }
    }
    # fwrite(one_fake_GWAS, file = paste0('/raid6/Tianyu/PRS/BootData/',anc,'_bootdata_repeat',i))
    fwrite(one_fake_GWAS, file = paste0('/raid6/Tianyu/PRS/BootData/',anc,'_bootdata_',setting.title,'_repeat',i))
    
  }
}  

print('pairwise association calculation')


