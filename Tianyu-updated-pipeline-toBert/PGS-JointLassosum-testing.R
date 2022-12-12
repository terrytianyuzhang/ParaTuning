####load in the combined-lasso beta, and get testing auc
gc()
options(stringsAsFactors = F)

library(data.table)
library(pryr) # check memory useage
library(lassosum) #transform p value to correlation
library(doParallel) # foreach
library(pROC) # for AUC 
library(R.utils)
library(snpStats)
library(wordspace)

PGS_bychr_bootstrap <- function(chr, anc, beta, work.dir){
  print(paste0('claculating PGS of ',anc, ' using chr ', chr))
  ######next step is generating some risk score using reference genotype
  ######we can download these genotype from 1000 Genome Project
  ######for the purpose of thee paper I will just use the reference panel 
  ######from which the training samples are generated
  
  chr_loc <- as.numeric(gsub(':.*$','',rownames(beta)))
  snp <- rownames(beta)[chr_loc == chr]
  
  ####load genotype data without label
  
  system.time(gnt<-read.plink(bed=paste0(work.dir, anc,".TST/CHR/",anc,".TST-chr",chr,".bed"),
                              bim=paste0(work.dir, anc,".TST/CHR/",anc,".TST-chr",chr,".bim"),
                              fam=paste0(work.dir, anc,".TST/CHR/",anc,".TST-chr",chr,".fam"),
                              select.snps = snp)
  )
  print(paste0('loaded plink file for ',anc, ' chr ', chr))
  
  #example : "/raid6/Tianyu/PRS/bert_sample/YRI.TUNE/CHR/YRI.TUNE-chr20.bed"
  
  # transform the genotypes to a matrix, and reverse the call count. This snpStats package counts the A2 allele, not the A1
  system.time(gnt <- 2 - as(gnt$genotypes,Class="numeric"))
  print(paste0('1', ' chr ', chr))
  
  # center the columns
  system.time(gnt <-gnt - rep(1, nrow(gnt)) %*% t(colMeans(gnt)))
  # colmean.gnt <- t(colMeans(gnt))
  # for(rowi in 1:NROW(gnt)){
  #   gnt[rowi,] <- gnt[rowi,] - colmean.gnt
  # }
  
  print(paste0('2', ' chr ', chr))
  # normalize the calls to the unit 1 norm
  system.time(gnt <- normalize.cols(gnt,method="euclidean",p=2))
  print(paste0('3', ' chr ', chr))
  # apply the proper shrinkage
  # system.time(gnt<-gnt*sqrt(1-shrink))
  
  # calculate the pgs
  system.time(re.pgs<-gnt %*% beta[chr_loc == chr,])
  print(paste0('finished claculating PGS of ',anc, ' using chr ', chr))
  # rm(gnt)
  return(re.pgs)
}

#### load the parameter information
load(paste0("Work/Sim-",i.sim,"/simulation-params.RData"))
main.dir=params$run.info$main.dir #"/raid6/Tianyu/PRS/SimulationPipeline/"
work.dir=params$run.info$work.dir #"/raid6/Tianyu/PRS/SimulationPipeline/Work/Sim-800/"

GAMMA = c(0.2, 0.5, 0.8)

####
# map=fread("/data3/DownLoadedData/GWAS-Populations-SimulationInput/chr1-22-qc-frq-ld.block.map",header=T,data.table=F)[,c("CHROM","ID")]
map <- fread(paste0(main.dir, 'Data/chr1-22-qc-frq-ld.block.map'),header=T,data.table=F)[,c("CHROM","ID")]

for(gamma in GAMMA){
  lasso.file <- paste(work.dir,"JointLassoSum/JointLassosum-",sprintf("-gamma-%.2f",gamma),".Rdata",sep="")
  
  ####now we load the beta0
  load(lasso.file)
  beta <- re.lasso$beta
  shrink <- re.lasso$shrink 
  lambda <- re.lasso$lambda
  
  
  # > head(map)
  # CHROM            ID
  # 1  chr1 1:1962845:T:C
  # 2  chr1 1:1962899:A:C
  # 3  chr1 1:1963406:G:A
  # 4  chr1 1:1963538:T:C
  # 5  chr1 1:1963738:C:T
  # 6  chr1 1:1964101:A:G
  
  
  CHR=gsub("chr","",map$CHROM)
  # > head(CHR)
  # [1] "1" "1" "1" "1" "1" "1"
  
  # SNP=map$ID
  SNP <- map[CHR %in% 1:22, ]$ID
  chrs <- 1:22 
  
  for(i.set in 1:2){
    
    anc <- if(i.set == 1){'CEU'} else 'YRI'
    ####read in the testing phenotype
    pheno <- fread(file = paste0(work.dir,anc,".TST/CHR/",anc,".TST-chr22.fam"))
    
    PGSnPHENO <- matrix(0, ncol = length(lambda) + 1, nrow = length(pheno$V6))
    PGSnPHENO[,length(lambda) + 1 ] <- pheno$V6
    
    rownames(beta) <- SNP
    
    
    re.pgss <- mclapply(chrs, PGS_bychr_bootstrap, anc = anc, 
                        beta = beta, work.dir = work.dir, mc.cores = 2,
                        mc.preschedule = F, mc.silent=F)
    
    ### sum the results
    pgs <- re.pgss[[1]] #this is the first chromosome
    for(i in 2:length(re.pgss)){
      pgs <- pgs+re.pgss[[i]]
    }
    PGSnPHENO[,1:length(lambda)] <- pgs
    
    save(PGSnPHENO, file = paste0(work.dir, "JointLassoSum/JointLassosum-",sprintf("-gamma-%.2f",gamma),"-", anc, ".TST-PGS.Rdata"))
  }
  
  
  #######load the PGS score and phenotype information, then calculate ROC
  # library(data.table)
  # library(pROC)
  for(i.set in 1:2){
    anc <- if(i.set == 1){'CEU'} else 'YRI'
    
    PGSnPHENO <- get(load(paste0(work.dir, "JointLassoSum/JointLassosum-",sprintf("-gamma-%.2f",gamma),"-", anc, ".TST-PGS.Rdata")))
    print(anc)
    nlambda <- NCOL(PGSnPHENO) - 1
    aucs <- rep(0, nlambda)
    ###AUC for each lambda
    for(i in 1:nlambda){
      aucs[i] <- auc(PGSnPHENO[, nlambda + 1], PGSnPHENO[,i])
    }
    save(aucs, file = paste0(work.dir, "JointLassoSum/JointLassosum-",sprintf("-gamma-%.2f",gamma),"-", anc, ".TST-AUC.Rdata"))
  }
}

rm.list=ls()
rm.list=rm.list[!(rm.list %in% c("i.sim","sims", "plink", "plink2"))]
rm(list = rm.list); flush.console()
