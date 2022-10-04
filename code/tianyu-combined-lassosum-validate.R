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
i.sim="NC8000"
# load(paste0("Work/Sim-C0.80-Y0.60-rep",i.sim,"/simulation-params.RData"))
load("/raid6/Ron/prs/data/bert_sample/simulation-params.RData")

# set directories for this simulation
# main.dir=params$run.info$main.dir
# work.dir=params$run.info$work.dir
main.dir <- "/raid6/Tianyu/PRS/"
work.dir <- "/raid6/Ron/prs/data/bert_sample/"



### load one of the lassosum runs
#load("Work/Sim-C0.80-Y0.60-rep1/CombinedLassoSum/Tmp/GWAS-lasso-C12000-Y12000-gamma-0.20.Rdata")

## next step is validation this needs to be done using the validate command, due to the transformation performed on the genotypes within lasso
## this transformation involves setting the mean to zero and the makingeach genotype vector have unit length (2-norm, sum of squares is 1)
## the approach will be to do this by chromosome, this means that each of the training/validation genotypes needs to be split by chromosome
## ughhh.



#system.time(re.validate <- validate(re.lasso, "Work/Sim-C0.80-Y0.60-rep1/CEU.TRN/CEU.TRN-4000", extract=NULL, keep=NULL, distance=NULL, mem.limit=mem.limit))

#save(re.validate,file="Tmp/re.validate-CEU-4000.RData")




pgsCalculationByCHR <- function(chr,work.dir,trn.set,trn.n,shrink,CHR,beta){
  # find the beta of interest
  snp=rownames(beta)[CHR == chr & rowSums(beta != 0) > 0]
  
  # /raid6/Ron/prs/data/bert_sample/GWAS-Populations-SimulationInput/CEU_reference_LDblocks/CHR
  #####!!!!!
  #I don't have these by chromsome data
  # system.time(gnt<-read.plink(bed=paste0(work.dir,trn.set,"/",trn.set,"-",trn.n,"-chr",chr,".bed"),
  #                             bim=paste0(work.dir,trn.set,"/",trn.set,"-",trn.n,"-chr",chr,".bim"),
  #                             fam=paste0(work.dir,trn.set,"/",trn.set,"-",trn.n,"-chr",chr,".fam"),select.snps = snp))
  system.time(gnt<-read.plink(bed=paste0("/raid6/Ron/prs/data/bert_sample/GWAS-Populations-SimulationInput/CEU_reference_LDblocks/CHR/CEU-chr",chr, ".bed"),
                              bim=paste0("/raid6/Ron/prs/data/bert_sample/GWAS-Populations-SimulationInput/CEU_reference_LDblocks/CHR/CEU-chr",chr, ".bim"),
                              fam=paste0("/raid6/Ron/prs/data/bert_sample/GWAS-Populations-SimulationInput/CEU_reference_LDblocks/CHR/CEU-chr",chr, ".fam"),
                              select.snps = snp))
  
  # > gnt$
  #   gnt$genotypes  gnt$fam        gnt$map
  
  # beforeflip = as(gnt$genotypes,Class="numeric")
  # > beforeflip[1:10,1:5]
  # 20:81154:T:G 20:82590:T:G 20:82603:A:C 20:83158:C:T 20:85729:G:A
  # CEU.1.REF             0            2            2            2            0
  # CEU.2.REF             0            2            2            2            0
  # CEU.3.REF             1            2            1            1            1
  # CEU.4.REF             0            2            2            2            0
  # CEU.5.REF             0            2            2            2            0
  # CEU.6.REF             0            2            2            2            0
  # CEU.7.REF             1            2            1            1            1
  # CEU.8.REF             1            2            1            1            1
  # CEU.9.REF             1            2            2            1            0
  # CEU.10.REF            2            2            2            0            0
  # 
  # transform the genotypes to a matrix, and reverse the call count. This snpStats package counts the A2 allele, not the A1
  system.time(gnt<-2-as(gnt$genotypes,Class="numeric"))
  
  # > gnt[1:10,1:5]
  # 20:81154:T:G 20:82590:T:G 20:82603:A:C 20:83158:C:T 20:85729:G:A
  # CEU.1.REF             2            0            0            0            2
  # CEU.2.REF             2            0            0            0            2
  # CEU.3.REF             1            0            1            1            1
  # CEU.4.REF             2            0            0            0            2
  # CEU.5.REF             2            0            0            0            2
  # CEU.6.REF             2            0            0            0            2
  # CEU.7.REF             1            0            1            1            1
  # CEU.8.REF             1            0            1            1            1
  # CEU.9.REF             1            0            0            1            2
  # CEU.10.REF            0            0            0            2            2
  # 
  # center the columns
  system.time(gnt <-gnt - rep(1, nrow(gnt)) %*% t(colMeans(gnt)))
  # normalize the calls to the unit 1 norm
  system.time(gnt<-normalize.cols(gnt,method="euclidean",p=2))
  # apply the proper shrinkage
  system.time(gnt<-gnt*sqrt(1-shrink))
  # calculate the pgs
  system.time(re.pgs<-gnt%*%beta[snp,])
  
  # > head(re.pgs)
  # lambda:0.0010 lambda:0.0014 lambda:0.0020 lambda:0.0029 lambda:0.0042
  # CEU.1.REF   0.001482600   0.001569183   0.001859870   0.001761504   0.001437728
  # CEU.2.REF  -0.003057378  -0.002774722  -0.002689785  -0.002222400  -0.001795640
  # CEU.3.REF   0.022345202   0.020316570   0.017967079   0.015036058   0.011865565
  # CEU.4.REF  -0.006479571  -0.006203878  -0.006129419  -0.005414084  -0.004386692
  # CEU.5.REF   0.007059047   0.006331442   0.005877239   0.005245275   0.004738714
  # CEU.6.REF  -0.006357727  -0.006249774  -0.006322660  -0.006122497  -0.005636031
  # lambda:0.0060 lambda:0.0085 lambda:0.0122 lambda:0.0175 lambda:0.0250
  # CEU.1.REF  0.0008671290  1.380483e-04 -0.0002583220 -2.549219e-04 -1.264908e-07
  # CEU.2.REF -0.0009455641 -6.904135e-05  0.0005535979  9.338254e-05 -1.264908e-07
  # CEU.3.REF  0.0087989511  5.199179e-03  0.0017507500  7.293846e-05 -1.264908e-07
  # CEU.4.REF -0.0030225186 -1.792922e-03 -0.0008500103 -1.800310e-04  4.363618e-06
  # CEU.5.REF  0.0037210673  2.610239e-03  0.0013154248  1.944960e-05 -1.264908e-07
  # CEU.6.REF -0.0046372074 -3.045054e-03 -0.0003247607  1.240125e-04  4.363618e-06
  
  # return the results
  return(re.pgs)
}




### some type of function  i.set,work.df,work.dir,SNP,CHR
pgsCalculation <- function(i.set,work.df=work.df,work.dir=work.dir,SNP=SNP,CHR=CHR){
  # set parameters
  lasso.file=work.df$lasso.file[i.set]
  trn.set=work.df$trn.set[i.set]
  trn.n=work.df$trn.n[i.set]
  
  # > lasso.file
  # [1] "/raid6/Tianyu/PRS/CombinedLassoSum/Tmp//GWAS-lasso-C20000-Y20000-gamma-0.50.Rdata"
  # > trn.set
  # [1] "CEU.TRN"
  # > trn.n
  # [1] 20000

  # load the lasso file and collect important information
  load(lasso.file)
  beta=re.lasso$beta
  rownames(beta)=SNP
  colnames(beta)=sprintf("lambda:%.4f",re.lasso$lambda)
  shrink=re.lasso$shrink
  rm(re.lasso)

  # run the validation
  system.time(re<-mclapply(1:22,pgsCalculationByCHR,work.dir=work.dir,trn.set=trn.set,trn.n=trn.n,CHR=CHR,shrink=shrink,beta=beta,mc.cores=10,mc.preschedule = F))
  
  
  # > head(re)
  # lambda:0.0010 lambda:0.0014 lambda:0.0020 lambda:0.0029 lambda:0.0042
  # CEU.1.REF   0.001482600   0.001569183   0.001859870   0.001761504   0.001437728
  # CEU.2.REF  -0.003057378  -0.002774722  -0.002689785  -0.002222400  -0.001795640
  # CEU.3.REF   0.022345202   0.020316570   0.017967079   0.015036058   0.011865565
  # CEU.4.REF  -0.006479571  -0.006203878  -0.006129419  -0.005414084  -0.004386692
  # CEU.5.REF   0.007059047   0.006331442   0.005877239   0.005245275   0.004738714
  # CEU.6.REF  -0.006357727  -0.006249774  -0.006322660  -0.006122497  -0.005636031
  # lambda:0.0060 lambda:0.0085 lambda:0.0122 lambda:0.0175 lambda:0.0250
  # CEU.1.REF  0.0008671290  1.380483e-04 -0.0002583220 -2.549219e-04 -1.264908e-07
  # CEU.2.REF -0.0009455641 -6.904135e-05  0.0005535979  9.338254e-05 -1.264908e-07
  # CEU.3.REF  0.0087989511  5.199179e-03  0.0017507500  7.293846e-05 -1.264908e-07
  # CEU.4.REF -0.0030225186 -1.792922e-03 -0.0008500103 -1.800310e-04  4.363618e-06
  # CEU.5.REF  0.0037210673  2.610239e-03  0.0013154248  1.944960e-05 -1.264908e-07
  # CEU.6.REF -0.0046372074 -3.045054e-03 -0.0003247607  1.240125e-04  4.363618e-06
  # 
  
  ### sum the results
  pgs=re[[1]] #this is the first chromosome
  for(i in 2:length(re)){
    pgs=pgs+re[[i]]
  }

  # write the results
  gwasset=unlist(strsplit(lasso.file,"/"))
  gwasset=gsub(".Rdata","",gsub("GWAS","",gwasset[length(gwasset)]))
  # > gwasset
  # [1] "-lasso-C20000-Y20000-gamma-0.50"
  fwrite(as.data.frame(pgs),paste0(work.dir,"CombinedLassoSum/PGS-",trn.set,"-",trn.n,gwasset,".pgs"),row.names=T,col.names=T,quote=F,sep="\t")

  return(i.set)
}





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
CHR=gsub("chr","",map[CHR == 20, ]$CHROM)
# > head(CHR)
# [1] "1" "1" "1" "1" "1" "1"
####!!!!!!!!!!
SNP=map$ID
SNP <- map[CHR == 20, ]$ID
rm(map)


# lasso.files=list.files(path=paste0(work.dir,"CombinedLassoSum/Tmp/"),pattern = ".Rdata",full.names = T)
lasso.files=list.files(path=paste0(main.dir,"CombinedLassoSum/Tmp/"),pattern = ".Rdata",full.names = T)


work.df=data.frame(trn.set=rep(c("CEU.TRN","YRI.TRN"),length(lasso.files)),trn.n=20000,lasso.file=rep(lasso.files,each=2))

# > work.df
# trn.set trn.n
# 1 CEU.TRN 20000
# 2 YRI.TRN 20000
# lasso.file
# 1 /raid6/Tianyu/PRS/CombinedLassoSum/Tmp//GWAS-lasso-C20000-Y20000-gamma-0.50.Rdata
# 2 /raid6/Tianyu/PRS/CombinedLassoSum/Tmp//GWAS-lasso-C20000-Y20000-gamma-0.50.Rdata

system.time(re.out<-mclapply(1:nrow(work.df),pgsCalculation,work.df=work.df,work.dir=work.dir,SNP=SNP,CHR=CHR,mc.cores=2,mc.preschedule = F))

