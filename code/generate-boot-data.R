#####there are 129,897 SNPs in chr20
#####dim(numeric_ref_genotypes)
#####[1]   5000 129897
set.seed(2019)
print('generate bootstrap data')
rm(list=ls()); gc()
options(stringsAsFactors = F)


#### software needed
plink <- "/usr/local/bin/plink"
plink2 <- "/usr/local/bin/plink2"

#### load the functions that are needed
source("R-code/simulation-functions.R")

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
library(snpStats)
library(wordspace)

#### load the parameter information
# i.sim="NC8000"
# load(paste0("Work/Sim-C0.80-Y0.60-rep",i.sim,"/simulation-params.RData"))
# load("/raid6/Ron/prs/data/bert_sample/simulation-params.RData")

# set directories for this simulation
# main.dir=params$run.info$main.dir
# work.dir=params$run.info$work.dir
main.dir <- "/raid6/Tianyu/PRS/"
work.dir <- "/raid6/Ron/prs/data/bert_sample/"

# lasso.files=list.files(path=paste0(main.dir,"CombinedLassoSum/Tmp/"),pattern = ".Rdata",full.names = T)
# work.df=data.frame(trn.set=rep(c("CEU.TRN","YRI.TRN"),length(lasso.files)),trn.n=20000,lasso.file=rep(lasso.files,each=2))

lasso.file <- "/raid6/Tianyu/PRS/CombinedLassoSum/Tmp/GWAS-lasso-C20000-Y20000-gamma-0.50.Rdata"

####now we load the beta0
load(lasso.file)
beta=re.lasso$beta
shrink <- re.lasso$shrink 

beta0 <- beta[,ceiling(NCOL(beta)/2)] #use the second smallest lambda to generate bootstrap data

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
# SNP=map$ID
SNP <- map[CHR %in% 20:21, ]$ID
names(beta0) <- SNP

B <- 2 ###B is the bootstrap repeats
chrs <- 20:21 #which chromosome did i use when training the model

PGS_bychr_bootstrap <- function(chr, anc, beta0, shrink){
  
  ######next step is generating some risk score using reference genotype
  ######we can download these genotype from 1000 Genome Project
  ######for the purpose of thee paper I will just use the reference panel 
  ######from which the training samples are generated
  
  chr_loc <- as.numeric(gsub(':.*$','',names(beta0)))
  snp <- names(beta0)[chr_loc == chr]
  
  ####load genotype data without label
  
  system.time(gnt<-read.plink(bed=paste0("/raid6/Tianyu/PRS/bert_sample/",anc,".TUNE/CHR/",anc,".TUNE-chr",chr,".bed"),
                              bim=paste0("/raid6/Tianyu/PRS/bert_sample/",anc,".TUNE/CHR/",anc,".TUNE-chr",chr,".bim"),
                              fam=paste0("/raid6/Tianyu/PRS/bert_sample/",anc,".TUNE/CHR/",anc,".TUNE-chr",chr,".fam"),
                              select.snps = snp)
  )
  
  # system.time(gnt<-read.plink(bed=paste0("/raid6/Ron/prs/data/bert_sample/GWAS-Populations-SimulationInput/",anc,"_reference_LDblocks/CHR/",anc,"-chr",chr, ".bed"),
  #                             bim=paste0("/raid6/Ron/prs/data/bert_sample/GWAS-Populations-SimulationInput/",anc,"_reference_LDblocks/CHR/",anc,"-chr",chr, ".bim"),
  #                             fam=paste0("/raid6/Ron/prs/data/bert_sample/GWAS-Populations-SimulationInput/",anc,"_reference_LDblocks/CHR/",anc,"-chr",chr, ".fam"),
  #                             select.snps = snp)
  # )
  gnt_map <- gnt$map
  # transform the genotypes to a matrix, and reverse the call count. This snpStats package counts the A2 allele, not the A1
  system.time(gnt<- 2 - as(gnt$genotypes,Class="numeric"))
  
  # center the columns
  system.time(gnt <-gnt - rep(1, nrow(gnt)) %*% t(colMeans(gnt)))
  # normalize the calls to the unit 1 norm
  system.time(gnt<-normalize.cols(gnt,method="euclidean",p=2))
  
  # apply the proper shrinkage
  # system.time(gnt<-gnt*sqrt(1-shrink))
  
  # calculate the pgs
  system.time(re.pgs<-gnt%*%beta0[chr_loc == chr])
  return(re.pgs)
}

cor_bychr_bootstrap <- function(chr, anc, bootY, beta0, shrink, B){
  
  ##calculate pairwise correlation
  
  chr_loc <- as.numeric(gsub(':.*$','',names(beta0)))
  snp <- names(beta0)[chr_loc == chr]
  
  ####load genotype data without label
  system.time(gnt<-read.plink(bed=paste0("/raid6/Tianyu/PRS/bert_sample/",anc,".TUNE/CHR/",anc,".TUNE-chr",chr,".bed"),
                              bim=paste0("/raid6/Tianyu/PRS/bert_sample/",anc,".TUNE/CHR/",anc,".TUNE-chr",chr,".bim"),
                              fam=paste0("/raid6/Tianyu/PRS/bert_sample/",anc,".TUNE/CHR/",anc,".TUNE-chr",chr,".fam"),
                              select.snps = snp)
  )
  
  # system.time(gnt<-read.plink(bed=paste0("/raid6/Ron/prs/data/bert_sample/GWAS-Populations-SimulationInput/",anc,"_reference_LDblocks/CHR/",anc,"-chr",chr, ".bed"),
  #                             bim=paste0("/raid6/Ron/prs/data/bert_sample/GWAS-Populations-SimulationInput/",anc,"_reference_LDblocks/CHR/",anc,"-chr",chr, ".bim"),
  #                             fam=paste0("/raid6/Ron/prs/data/bert_sample/GWAS-Populations-SimulationInput/",anc,"_reference_LDblocks/CHR/",anc,"-chr",chr, ".fam"),
  #                             select.snps = snp)
  # )
  
  gnt_map <- gnt$map
  # transform the genotypes to a matrix, and reverse the call count. This snpStats package counts the A2 allele, not the A1
  system.time(gnt<-2-as(gnt$genotypes,Class="numeric"))
  
  # center the columns
  system.time(gnt <-gnt - rep(1, nrow(gnt)) %*% t(colMeans(gnt)))
  # normalize the calls to the unit 1 norm
  system.time(gnt<-normalize.cols(gnt,method="euclidean",p=2))
  
  ###calculate pairwise correlation
  print('calculating pairwise correaltion')
  boot.cor <- matrix(0, nrow = NCOL(gnt), ncol = B)
  for(i in 1:B){
    print(paste0('the No.',i,'boostrap data'))
    for(j in 1:NCOL(gnt)){
      ###after the proper normalization, correlation is the same as inner product
      boot.cor[j,i] <- crossprod(boot.Y[,i], gnt[,j] )
    }
  }
  
  #####make some "fake" GWAS results
  print('translate correlation into p-values')
  fake_GWAS_list <- list()
  sample_size <- NROW(gnt)
  for(i in 1:B){
    t_stat <- boot.cor[,i] * sqrt(sample_size - 2)/sqrt(1 - (boot.cor[,i])^2)
    fake_GWAS <- data.frame(P = 2*pt(q = abs(t_stat), df = sample_size - 2, lower.tail = F), #two-sided p-value
                            OR = 2*as.numeric(boot.cor[,i]>0) + 0.5*as.numeric(boot.cor[,i]<0))
    fake_GWAS$n <- sample_size
    fake_GWAS_list[[i]] <- fake_GWAS
    # fwrite(fake_GWAS, file = paste0('/raid6/Tianyu/PRS/BootData/',anc,'_bootdata_repeat',i,'_chr',chr))
  }
  
  return(fake_GWAS_list)
}

set.seed(2019)

for(i.set in 1:2){
    
  # lasso.file=work.df$lasso.file[i.set]
  # trn.set=work.df$trn.set[i.set]
  # trn.n=work.df$trn.n[i.set]
  
  anc <- if(i.set == 1) 'CEU' else 'YRI'
  
  re.pgss <- mclapply(chrs, PGS_bychr_bootstrap, anc = anc, 
           beta0 = beta0, shrink = shrink)
  
  ### sum the results
  pgs <- re.pgss[[1]] #this is the first chromosome
  for(i in 2:length(re.pgss)){
    pgs <- pgs+re.pgss[[i]]
  }
  
  ###pgs is the risk score for each subject in the reference panel
  
  #####we need to read figure out what is the original data noise level
  
  pheno <- fread(file = paste0("/raid6/Ron/prs/data/bert_sample/",anc,".TRN/",anc,".TRN.psam"))
  summary(as.factor(pheno$PHENO1))
  
  
  ##the summary stat of phenotype deteremine how we calibrate the score into probability
  
  # > summary(as.factor(pheno$PHENO1))
  # 1    2 
  # 2000 2000 
  #####calibrate the PRS to a probability, this is 1:1 case control setting
  re.pgs.prob <- (pgs - min(pgs))/ (max(pgs) - min(pgs))
  #!!!!!!!!!
  re.pgs.prob <- 0.5 + (re.pgs.prob - 0.5)*0
  #### SD of the noise variable
  # sd <- sqrt(abs(re.pgs))
  
  ### generate bootstrap Y
  
  boot.Y <- matrix(0, nrow = length(re.pgs.prob), ncol = B)
  for(i in 1:length(re.pgs.prob)){
    boot.Y[i,] <- rbinom(B, size = 1, prob = re.pgs.prob[i])
  }
  
  # center and normalize the boostrap Y
  boot.Y <- boot.Y - rep(1, nrow(boot.Y)) %*% t(colMeans(boot.Y))
  
  boot.Y <- normalize.cols(boot.Y, method="euclidean",p=2)
  
  # save(boot.Y, 
  #      file = paste0('/raid6/Tianyu/PRS/BootData/',anc,'_bootY'))
  
  save(boot.Y, 
       file = paste0('/raid6/Tianyu/PRS/BootData/',anc,'_bootY_lowsignal'))
  
  # cor_bychr_bootstrap(chr = chrs[1], anc = anc, bootY = bootY, 
  #                     beta0 = beta0, shrink = shrink, B = B)
  fake_GWAS <- mclapply(chrs, cor_bychr_bootstrap, anc = anc, bootY = bootY, 
                      beta0 = beta0, shrink = shrink, B = B)
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
    
    fwrite(one_fake_GWAS, file = paste0('/raid6/Tianyu/PRS/BootData/',anc,'_bootdata_lowsignal_repeat',i))
  }
  # '/raid6/Tianyu/PRS/BootData/CEU_bootdata_repeat1'
  
  # 
  # ######next step is generating some risk score using reference genotype
  # ######we can download these genotype from 1000 Genome Project
  # ######for the purpose of thee paper I will just use the reference panel 
  # ######from which the training samples are generated
  # chr_loc <- as.numeric(gsub(':.*$','',names(beta0)))
  # snp <- names(beta0)[chr_loc == chr]
  # 
  # ####load genotype data without label
  # system.time(gnt<-read.plink(bed=paste0("/raid6/Ron/prs/data/bert_sample/GWAS-Populations-SimulationInput/",anc,"_reference_LDblocks/CHR/",anc,"-chr",chr, ".bed"),
  #                             bim=paste0("/raid6/Ron/prs/data/bert_sample/GWAS-Populations-SimulationInput/",anc,"_reference_LDblocks/CHR/",anc,"-chr",chr, ".bim"),
  #                             fam=paste0("/raid6/Ron/prs/data/bert_sample/GWAS-Populations-SimulationInput/",anc,"_reference_LDblocks/CHR/",anc,"-chr",chr, ".fam"),
  #                             select.snps = snp)
  # )
  # gnt_map <- gnt$map
  # # transform the genotypes to a matrix, and reverse the call count. This snpStats package counts the A2 allele, not the A1
  # system.time(gnt<-2-as(gnt$genotypes,Class="numeric"))
  # 
  # # center the columns
  # system.time(gnt <-gnt - rep(1, nrow(gnt)) %*% t(colMeans(gnt)))
  # # normalize the calls to the unit 1 norm
  # system.time(gnt<-normalize.cols(gnt,method="euclidean",p=2))
  # # apply the proper shrinkage
  # system.time(gnt<-gnt*sqrt(1-shrink))
  # # calculate the pgs
  # system.time(re.pgs<-gnt%*%beta0[chr_loc == chr])
  # 
  # #####now re.pgs is our preliminary estimate of the true risk score
  # #####add some noise to construct Y
  # 
  
  
  
  # ### calculate pairwise p value of correlation coefficient
  # boot.cor <- matrix(0, nrow = NCOL(gnt), ncol = B)
  # for(i in 1:B){
  #   print(paste0("generating the No.",i ," bootstrap data"))
  #   for(j in 1:NCOL(gnt)){
  #     boot.cor[j,i] <- cor(boot.Y[,i], gnt[,j])
  #   }
  # }
  # 
  # fwrite(boot.cor, file = paste0('/raid6/Tianyu/PRS/BootData/',anc,'_TRN_',anc,'_reference_chr', chr,'_bootcor' ))
  # 
  # #####make some "fake" GWAS results
  # sample_size <- NROW(gnt)
  # for(i in 1:B){
  #   t_stat <- boot.cor[,i] * sqrt(sample_size - 2)/sqrt(1 - (boot.cor[,i])^2)
  #   fake_GWAS <- data.frame(P = pt(q = abs(t_stat), df = sample_size - 2, lower.tail = F),
  #                           OR = 2*as.numeric(boot.cor[,i]>0) + 0.5*as.numeric(boot.cor[,i]<0))
  #   fake_GWAS$n <- sample_size
  #   fwrite(fake_GWAS, file = paste0('/raid6/Tianyu/PRS/BootData/',anc,'_bootdata_repeat',i,'_chr',chr))
  # }
}

# 
# org_GWAS_results <- fread(file = "YRI.TRN.PHENO1.glm.logistic.hybrid", header=T, data.table=F)
# 
# reference_id <- fread(paste0("/raid6/Ron/prs/data/bert_sample/GWAS-Populations-SimulationInput/YRI_reference_LDblocks/CHR/YRI-chr",chr, ".fam"))
# reference_id[,6] <- boot.Y[,1] + 1 ## '1' = control, '2' = case for fam files
# fwrite(reference_id, file = paste0('/raid6/Tianyu/PRS/BootData/YRI_TRN_YRI_reference_chr', chr,'.fam' ))


# #####create a fake \hat \beta_0
# set.seed(1)
# chr20snpN <- 129897
# beta0 <- rep(0, chr20snpN)
# non_zero_snp <- sample(1:chr20snpN, 100)
# beta0[non_zero_snp] <- 1
# 
# ######use a loop to calculate tau and C\beta #####
# # install.packages('RSpectra')
# library(RSpectra)
# library(lassosum)
# library(data.table)
# library(snpStats)
# setwd('/home/tianyuz3/PRS/my_code/')
# 
# chr_index <- 'chr20'
# ld.anc <- 'AFR' ###ancestor
# sample_size <- 5000
# if(ld.anc == 'AFR'){
#   anc_initial <- 'YRI'
# }
# 
# SVD_big_list <- get(load(paste0("/raid6/Tianyu/PRS/SVDdata/", anc_initial, "reference_LDblocks_", chr_index, "_SVD.RData")))
# 
# #nothing more than doing a matrix, vector multiplication
# pointer_begin <- pointer_end <- 1
# Cbeta <- NULL
# for(i in 1:length(SVD_big_list)){
#   tempC <- SVD_big_list[[i]]$C
#   
#   pointer_end <- pointer_begin + dim(tempC)[1] - 1
#   Cbeta <- c(Cbeta, tempC %*% beta0[pointer_begin:pointer_end])
#   
#   pointer_begin <- pointer_end + 1
# }
# 
# ###these quantities will be used during boostrap
# nCbeta <- sample_size * Cbeta
# tau <- t(beta0) %*% as.matrix(Cbeta, ncol = 1)
# 
# ####save files
# pre_cal <- list(tau = tau, nCbeta = nCbeta)
# save(pre_cal, file = paste0("/raid6/Tianyu/PRS/SVDdata/", anc_initial, "pre_cal_for_bootstrap.RData"))
