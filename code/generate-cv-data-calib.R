###########generate boostrap samples and split them into training and validation######
set.seed(2019)
print('generate bootstrap data for cv')
rm(list=ls()); gc()
options(stringsAsFactors = F)

#### software needed
plink <- "/usr/local/bin/plink"
plink2 <- "/usr/local/bin/plink2"

setwd("/home/tianyuz3/PRS/my_code")

library(data.table)
library(pryr) # check memory useage
library(lassosum) #transform p value to correlation
library(doParallel) # foreach
library(pROC) # for AUC 
library(R.utils)
library(snpStats)
library(wordspace)

setting.title <- '1in2signal'
chrs <- 20:21 #which chromosome did i use when training the model

nfold <- 5
##########read in lasso results##########
lasso.file <- "/raid6/Tianyu/PRS/CombinedLassoSum/Tmp/GWAS-lasso-C20000-Y4000-gamma-0.50.Rdata"

re.lasso <- get(load(lasso.file))
beta <- re.lasso$beta

beta0 <- beta[,ceiling(NCOL(beta)/2)] #use a small lambda to generate bootstrap data
##########read in genotype and calculate score#############
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
CHR <- gsub("chr","",map$CHROM)
# > head(CHR)
# [1] "1" "1" "1" "1" "1" "1"
####!!!!!!!!!!
# SNP=map$ID
SNP <- map[CHR %in% 20:21, ]$ID
names(beta0) <- SNP

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
  #example : "/raid6/Tianyu/PRS/bert_sample/YRI.TUNE/CHR/YRI.TUNE-chr20.bed"
  
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

risk.score.list <- vector("list",2)

####uncomment from here
for(i.set in 1:2){
  
  if(i.set == 1){
    anc <- 'CEU'
  }else{
    anc <- 'YRI'
  }
  
  re.pgss <- mclapply(chrs, PGS_bychr_bootstrap, anc = anc,
                      beta0 = beta0, shrink = shrink, mc.cores = 4)
  
  ### sum the results
  pgs <- re.pgss[[1]] #this is the first chromosome
  for(i in 2:length(re.pgss)){
    pgs <- pgs+re.pgss[[i]]
  }
  
  risk.score.list[[i.set]] <- pgs
  
}
### save risk scores for each population
save(risk.score.list, file = '/raid6/Tianyu/PRS/CombinedLassoSum/Tmp/GWAS-lasso-C20000-Y4000-gamma-0.50-riskscore.Rdata')
###pgs is the risk score for each subject in the reference panel
######uncomment


###########generate Y##############
risk.score.list <- get(load('/raid6/Tianyu/PRS/CombinedLassoSum/Tmp/GWAS-lasso-C20000-Y4000-gamma-0.50-riskscore.Rdata'))
#####we need to read figure out what is the original data noise level
anc <- 'CEU'
pheno <- fread(file = paste0("/raid6/Ron/prs/data/bert_sample/",anc,".TRN/",anc,".TRN.psam"))
#example: "/raid6/Ron/prs/data/bert_sample/CEU.TRN/CEU.TRN.psam"
summary(as.factor(pheno$PHENO1))

s.size <- length(pheno$PHENO1)
s2 <- sum(risk.score.list[[1]]^2)
#this is the calibration slope a
ceuorg <- fread('/raid6/Ron/prs/data/bert_sample/CEU.TRN.PHENO1.glm.logistic.hybrid')

names(ceuorg)[1] <- 'chr'
ceuorg <- ceuorg[chr %in% c(20,21),]
summary(ceucor.org <- p2cor(p = ceuorg$P, n = 20000,sign = log(ceuorg$OR)))

a <- ceucor.org %*% re.lasso$beta[,5] * sqrt(s.size)/(2 * s2)
a <- as.numeric(a)
a <- a/2
mu <-  a * risk.score.list[[1]] + 0.5

print(paste0('number of beyond 0/1 range', sum(mu > 1 | mu<0)))
mu[mu > 1] <- 1
mu[mu < 0] <- 0

epsilon <- rbinom(s.size,1,mu)
epsilon <- as.numeric(epsilon == 1)*(1-mu) + as.numeric(epsilon == 0)*(-mu)
booty <- a*risk.score.list[[1]] + epsilon

##this value should be close to the original r*beta
save(booty, file = paste0('/raid6/Tianyu/PRS/BootData/',anc,'_bootY_',setting.title,'.RData'))
##"/raid6/Tianyu/PRS/BootData/CEU_bootY_calib.RData"

#####we need to read figure out what is the original data noise level
anc <- 'YRI'
pheno <- fread(file = paste0("/raid6/Ron/prs/data/bert_sample/",anc,".TRN/",anc,".TRN.psam"))
#example: "/raid6/Ron/prs/data/bert_sample/CEU.TRN/CEU.TRN.psam"
summary(as.factor(pheno$PHENO1))

s.size <- length(pheno$PHENO1)
s2 <- sum(risk.score.list[[2]]^2)

yriorg <- fread('/raid6/Ron/prs/data/bert_sample/YRI.TRN.PHENO1.glm.logistic.hybrid')

names(yriorg)[1] <- 'chr'
yriorg <- yriorg[chr %in% c(20,21),]
summary(yricor.org <- p2cor(p = yriorg$P, n = 4000,sign = log(yriorg$OR)))

#this is the calibration slope a
a <- yricor.org %*% re.lasso$beta[,5] * sqrt(s.size)/(2 * s2)
a <- as.numeric(a)
a <- a/2
mu <-  a * risk.score.list[[2]] + 0.5

print(paste0('number of beyond 0/1 range', sum(mu > 1 | mu<0)))
mu[mu > 1] <- 1
mu[mu < 0] <- 0

epsilon <- rbinom(s.size,1,mu)
epsilon <- as.numeric(epsilon == 1)*(1-mu) + as.numeric(epsilon == 0)*(-mu)
booty <- a*risk.score.list[[2]] + epsilon

##this value should be close to the original r*beta
save(booty, file = paste0('/raid6/Tianyu/PRS/BootData/',anc,'_bootY_',setting.title,'.RData'))



##########split training and testing########
set.seed(2019)
for(anc in c("CEU", "YRI")){
  if(anc == "CEU"){
  s.size <-  20000
  }else{
  s.size <- 4000
  }
  
  val.index <- sort(sample(1:s.size, floor(s.size/nfold)))
  save(val.index, file = paste0("/raid6/Tianyu/PRS/BootData/",anc,".TUNE/val_index.RData"))
  
  train.index <- (1:s.size)[-val.index] #this is in order
  save(train.index, file = paste0("/raid6/Tianyu/PRS/BootData/",anc,".TUNE/train_index.RData"))

  ####generate training and testing .fam files

  for(chr in c(20,21)){
    full.fam <- fread(paste0("/raid6/Tianyu/PRS/bert_sample/",anc,".TUNE/CHR/",anc,".TUNE-chr",chr,".fam"))
    train.fam <- full.fam[train.index, c(1,2)] #only keep family id and withtin family id
    fwrite(train.fam, paste0("/raid6/Tianyu/PRS/BootData/",anc,".TUNE/CHR/",anc,".TUNE-boost-train-index.txt"),
           col.names = F, sep = " ")
    
    plink.command <- paste(plink, "--bfile", paste0("/raid6/Tianyu/PRS/bert_sample/",anc,".TUNE/CHR/",anc,".TUNE-chr",chr),
                           "--allow-no-sex",
                           "--keep", paste0("/raid6/Tianyu/PRS/BootData/",anc,".TUNE/CHR/",anc,".TUNE-boost-train-index.txt"),
                           "--make-bed", "--out", paste0("/raid6/Tianyu/PRS/BootData/", anc,".TUNE/CHR/",anc,".TUNE-chr",chr, "boost-train"),
                           "--noweb",
                           sep = " ")
    system(plink.command)
    
    val.fam <- full.fam[val.index, c(1,2)]
    fwrite(val.fam, paste0("/raid6/Tianyu/PRS/BootData/",anc,".TUNE/CHR/",anc,".TUNE-boost-val-index.txt"),
           col.names = F, sep = " ")
    
    plink.command <- paste(plink, "--bfile", paste0("/raid6/Tianyu/PRS/bert_sample/",anc,".TUNE/CHR/",anc,".TUNE-chr",chr),
                           "--allow-no-sex",
                          "--keep", paste0("/raid6/Tianyu/PRS/BootData/",anc,".TUNE/CHR/",anc,".TUNE-boost-val-index.txt"),
                          "--make-bed", "--out", paste0("/raid6/Tianyu/PRS/BootData/", anc,".TUNE/CHR/",anc,".TUNE-chr",chr, "boost-val"),
                          "--noweb",
                          sep = " ")
    system(plink.command)
  }
}

#########save training association########
cor_bychr_bootstrap <- function(chr, anc, boot.Y, B){
  
  ##calculate pairwise correlation
  
  # chr_loc <- as.numeric(gsub(':.*$','',names(beta0)))
  # snp <- names(beta0)[chr_loc == chr]
  
  ####load genotype data without label
  system.time(gnt<-read.plink(bed=paste0("/raid6/Tianyu/PRS/BootData/",anc,".TUNE/CHR/",anc,".TUNE-chr",chr,"boost-train.bed"),
                              bim=paste0("/raid6/Tianyu/PRS/BootData/",anc,".TUNE/CHR/",anc,".TUNE-chr",chr,"boost-train.bim"),
                              fam=paste0("/raid6/Tianyu/PRS/BootData/",anc,".TUNE/CHR/",anc,".TUNE-chr",chr,"boost-train.fam"))
  )
  
  gnt_map <- gnt$map
  # transform the genotypes to a matrix, and reverse the call count. This snpStats package counts the A2 allele, not the A1
  system.time(gnt<-2-as(gnt$genotypes,Class="numeric"))
  
  # center the columns
  system.time(gnt <-gnt - rep(1, nrow(gnt)) %*% t(colMeans(gnt)))
  # normalize the calls to the unit 1 norm
  system.time(gnt<-normalize.cols(gnt,method="euclidean",p=2))
  
  ###calculate pairwise correlation
  print('calculating pairwise correaltion')
  print(paste0('B is ', B))
  boot.cor <- matrix(0, nrow = NCOL(gnt), ncol = B)
  for(i in 1:B){
    print(paste0('the No.',i,'boostrap data of ', anc))
    for(j in 1:NCOL(gnt)){
      ###after the proper normalization, correlation is the same as inner product
      boot.cor[j,i] <- crossprod(boot.Y[,i], gnt[,j] )
    }
  }
  # save(boot.cor, file = '/raid6/Tianyu/PRS/trash/beta0is5-YRI-TRN-chr20-cor.RData')
  
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

######
B <- 1
chrs <- 20:21
for(i.set in 1:2){
  if(i.set == 1){
    anc <- 'CEU'
  }else{
    anc <- 'YRI'
  }
  
  booty <- get(load(paste0("/raid6/Tianyu/PRS/BootData/",anc,"_bootY_",setting.title,".RData")))
  train.index <- get(load(paste0("/raid6/Tianyu/PRS/BootData/",anc,".TUNE/train_index.RData")))
  booty <- booty[train.index]
  # ceubooty <- get(load("/raid6/Tianyu/PRS/BootData/CEU_bootY_calib.RData"))
  # yribooty <- get(load("/raid6/Tianyu/PRS/BootData/YRI_bootY_calib.RData"))
  
  boot.Y <- matrix(booty, ncol = 1)
  # center and normalize the boostrap Y
  boot.Y <- boot.Y - rep(1, nrow(boot.Y)) %*% t(colMeans(boot.Y))
  
  boot.Y <- normalize.cols(boot.Y, method="euclidean",p=2)
  
  fake_GWAS <- mclapply(chrs, cor_bychr_bootstrap, anc = anc, boot.Y = boot.Y, 
                        B = B)
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



