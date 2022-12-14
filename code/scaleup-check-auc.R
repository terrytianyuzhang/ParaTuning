####load in the combined-lasso beta, and get testing auc
set.seed(2019)
rm(list=ls()); gc()
options(stringsAsFactors = F)

# #### software needed
# plink <- "/usr/local/bin/plink"
# plink2 <- "/usr/local/bin/plink2"
# 
# #### load the functions that are needed
# source("R-code/simulation-functions.R")
# 
# ### Load functions ###
# source(paste0("R-code/LassosumFunctions/parseselect.R"))
# source(paste0("R-code/LassosumFunctions/parseblocks.R"))
# source(paste0("R-code/LassosumFunctions/ncol.bfile.R"))
# source(paste0("R-code/LassosumFunctions/nrow.bfile.R"))
# source(paste0("R-code/LassosumFunctions/read.table2.R"))
# source(paste0("R-code/LassosumFunctions/selectregion.R"))
# source(paste0("R-code/LassosumFunctions/parse.pheno.covar.R"))
# source(paste0("R-code/LassosumFunctions/myelnet.R"))
# source(paste0("R-code/LassosumFunctions/mylassosum.R"))
# source(paste0("R-code/LassosumFunctions/splitgenome.R"))
# source(paste0("R-code/LassosumFunctions/validation.R"))
# source(paste0("R-code/LassosumFunctions/merge.mylassosum.R"))
# Rcpp::sourceCpp(paste0("R-code/LassosumFunctions/myfunctions.cpp"))

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
setting.title <- 'CEU1aYRI2a22Chr'
main.dir <- "/raid6/Tianyu/PRS/"
work.dir <- "/raid6/Ron/prs/data/bert_sample/"

# lasso.files=list.files(path=paste0(main.dir,"CombinedLassoSum/Tmp/"),pattern = ".Rdata",full.names = T)
# work.df=data.frame(trn.set=rep(c("CEU.TRN","YRI.TRN"),length(lasso.files)),trn.n=20000,lasso.file=rep(lasso.files,each=2))

lasso.file <- paste0("/raid6/Tianyu/PRS/CombinedLassoSum/Tmp/GWAS-lasso-C20000-Y4000-gamma-0.50_", setting.title,".Rdata")

####now we load the beta0
load(lasso.file)
beta=re.lasso$beta
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
# SNP=map$ID
SNP <- map[CHR %in% 1:22, ]$ID


# for(i.set in 1:2){

# lasso.file=work.df$lasso.file[i.set]
# trn.set=work.df$trn.set[i.set]
# trn.n=work.df$trn.n[i.set]
PGS_bychr_bootstrap <- function(chr, anc, beta){
  print(paste0('claculating PGS of ',anc, 'using chr ', chr))
  ######next step is generating some risk score using reference genotype
  ######we can download these genotype from 1000 Genome Project
  ######for the purpose of thee paper I will just use the reference panel 
  ######from which the training samples are generated
  
  chr_loc <- as.numeric(gsub(':.*$','',rownames(beta)))
  snp <- rownames(beta)[chr_loc == chr]
  
  ####load genotype data without label
  
  system.time(gnt<-read.plink(bed=paste0("/raid6/Tianyu/PRS/bert_sample/",anc,".TUNE/CHR/",anc,".TUNE-chr",chr,".bed"),
                              bim=paste0("/raid6/Tianyu/PRS/bert_sample/",anc,".TUNE/CHR/",anc,".TUNE-chr",chr,".bim"),
                              fam=paste0("/raid6/Tianyu/PRS/bert_sample/",anc,".TUNE/CHR/",anc,".TUNE-chr",chr,".fam"),
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
  print(paste0('finished claculating PGS of ',anc, 'using chr ', chr))
  
  return(re.pgs)
}

for(i.set in 1:2){
  
  anc <- if(i.set == 1){'CEU'} else 'YRI'
  
  pheno <- fread(file = paste0("/raid6/Tianyu/PRS/bert_sample/",anc,".TUNE/CHR/",anc,".TUNE-chr20.fam"))
  
  PGSnPHENO <- matrix(0, ncol = length(lambda) + 1, nrow = length(pheno$V6))
  PGSnPHENO[,length(lambda) + 1 ] <- pheno$V6
  
  rownames(beta) <- SNP
  chrs <- 1:22 
  
  re.pgss <- mclapply(chrs, PGS_bychr_bootstrap, anc = anc, 
                      beta = beta, mc.cores = 4)
  
  ### sum the results
  pgs <- re.pgss[[1]] #this is the first chromosome
  for(i in 2:length(re.pgss)){
    pgs <- pgs+re.pgss[[i]]
  }
  PGSnPHENO[,1:length(lambda)] <- pgs
  
  save(PGSnPHENO, file = paste0('/raid6/Tianyu/PRS/trash/', setting.title,'PGSnPHENO_org_',anc, '.RData'))
}
#######load the PGS score and phenotype information, then calculate ROC
library(data.table)
library(pROC)
setting.title <- 'CEU1aYRI2a22Chr'
anc <- 'CEU'
PGSnPHENO <- get(load(paste0('/raid6/Tianyu/PRS/trash/', setting.title,'PGSnPHENO_org_',anc, '.RData')))
print(anc)
nlambda <- NCOL(PGSnPHENO) - 1
aucs <- rep(0, nlambda)
for(i in 1:nlambda){
  aucs[i] <- auc(PGSnPHENO[, nlambda + 1], PGSnPHENO[,i])
}
save(aucs, file = paste0('/raid6/Tianyu/PRS/trash/', setting.title,'PGSnPHENO_org_auc_',anc, '.RData'))

anc <- 'YRI'
PGSnPHENO <- get(load(paste0('/raid6/Tianyu/PRS/trash/', setting.title,'PGSnPHENO_org_',anc, '.RData')))
print(anc)
nlambda <- NCOL(PGSnPHENO) - 1
aucs <- rep(0, nlambda)
for(i in 1:nlambda){
  aucs[i] <- auc(PGSnPHENO[, nlambda + 1], PGSnPHENO[,i])
}
save(aucs, file = paste0('/raid6/Tianyu/PRS/trash/', setting.title,'PGSnPHENO_org_auc_',anc, '.RData'))




  

  
# }