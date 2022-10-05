rm(list=ls()); 
gc()
options(stringsAsFactors = F)
#####I am trying to fit combined lassosum algorithm

#########
### The strategy to determine the combined lassosum is to use the architecture of the computer to determine the solutions
### It works from the original Peng code without any modifications to increase the speed.
### In this setting we do a grid search on 5 values of gamma and 5 for lambda
### The code accepts multiple lambdas in pone run, but swill only use a single gamma
### For each complete run we need to run the program, 5 times.
### By calculating the lassosum by chromososome there is an additional level of parallelization that can be achieved.
### After some experimentation I came to the conclusion that the quickest approach is to use 12 processors for the 22 chromosomes.
### For this I used a memory limit per processor of 20Gb
### This allows the smaller chromosomes to run on the same processor and be done at the same times the smaller ones.
### A single run can then be completed in 20 minutes using 1+12 processors at a memory cost of ~260Gb. On our 2 Tb computer with 192 cores 
### we can run 7 gammas simultaneously. I opted to use 5, to allow me some memory to do additional work.

#### software needed
# plink="/data3/Software/Plink/plink"
# plink2="/data3/Software/Plink2/plink2"
plink <- "/usr/local/bin/plink"
plink2 <- "/usr/local/bin/plink2"

setwd("/home/tianyuz3/PRS/my_code")
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





#### Functions to use
mylassosumFunction <- function(chr,gamma,lambda,shrink,
                               mem.limit,gwasANC,COR,
                               referenceFiles,LDblocks){
  re=mylassosum(cor1=COR[COR$CHR == chr,gwasANC[1]], 
                cor2=COR[COR$CHR == chr,gwasANC[2]],
                fileName1 = referenceFiles[[chr]][[gwasANC[1]]], 
                fileName2 = referenceFiles[[chr]][[gwasANC[2]]],  
                gamma = gamma, lambda = lambda, shrink=shrink,
                chunk=TRUE, mem.limit=mem.limit,
                trace=2,
                LDblocks1=as.data.frame(LDblocks[[gwasANC[1]]][LDblocks[[gwasANC[1]]]$chr == paste0("chr",chr),]), 
                LDblocks2=as.data.frame(LDblocks[[gwasANC[2]]][LDblocks[[gwasANC[2]]]$chr == paste0("chr",chr),]))
  return(re)
}

##### general setup
# i.sim="NC8000"
# # load(paste0("Work/Sim-C0.80-Y0.60-rep",i.sim,"/simulation-params.RData"))
# load("/raid6/Ron/prs/data/bert_sample/simulation-params.RData")

# set directories for this simulation
# main.dir=params$run.info$main.dir
# work.dir=params$run.info$work.dir

main.dir <- "/raid6/Tianyu/PRS/"
work.dir <- "/raid6/Ron/prs/data/bert_sample/"

### ancestries
gwasANC=c("CEU","YRI")
#### set the parameters
GAMMA=0.5
#!!!!!! need to change this back
# lambda=exp(seq(log(0.001), log(0.025), length.out=10))
lambda=exp(seq(log(0.001), log(0.025), length.out=2))
shrink=.9

N=seq(20000,20000,4000)
input.df=data.frame(N1=rep(N,each=length(N)),N2=rep(N,length(N)))
input.df=data.frame(gamma=rep(GAMMA,each=nrow(input.df)),N1=rep(input.df$N1,length(GAMMA)),N2=rep(input.df$N2,length(GAMMA)))

# memory limit
mem.limit=2*10e9

wrapperFunction <- function(i.combn, input.df, gwasANC, lambda, shrink, main.dir, work.dir, CHR=1:22, mem.limit, mc.cores = 12){

  # collect gamma
  gamma=input.df$gamma[i.combn]  

  # collectthe number of samples for the two populations
  gwasN=list()
  gwasN[["CEU"]]=input.df$N1[i.combn]
  gwasN[["YRI"]]=input.df$N2[i.combn]
  
  # > gwasN
  # $CEU
  # [1] 20000
  # 
  # $YRI
  # [1] 20000


  #### determine correlations
  ### Summary statistics ###
  # read a map for the first ancestry to set up a dataframe
  anc=gwasANC[1]
  map <- fread(paste0(work.dir, anc, '.TRN/', anc, '.TRN.pvar'), header=T, data.table=F)
  # map=fread(paste0(work.dir,anc,".GWAS/",anc,".GWAS-",gwasN[[anc]],".pvar"),header=T,data.table=F)
  # > head(map)
  # #CHROM     POS            ID REF ALT
  # 1      1 1962845 1:1962845:T:C   T   C
  # 2      1 1962899 1:1962899:A:C   A   C
  # 3      1 1963406 1:1963406:G:A   G   A
  # 4      1 1963538 1:1963538:T:C   T   C
  # 5      1 1963738 1:1963738:C:T   C   T
  # 6      1 1964101 1:1964101:A:G   A   G
  
  
  COR=data.frame(CHR=map$`#CHROM`,ID=map$ID)
  
  #####!!!!!
  COR <- COR[COR$CHR == 20,]
  # ##> head(COR)
  # CHR            ID
  # 1   1 1:1962845:T:C
  # 2   1 1:1962899:A:C
  # 3   1 1:1963406:G:A
  # 4   1 1:1963538:T:C
  # 5   1 1:1963738:C:T
  # 6   1 1:1964101:A:G
  
  for(anc in gwasANC){
    # in bert's computer
    # glm=fread(paste0(work.dir,anc,".GWAS/Assoc/",anc,".GWAS-",gwasN[[anc]],".PHENO1.glm.logistic.hybrid"),header=T,data.table=F)
    # when training the original data
    # glm <- fread(paste0(work.dir, anc, '.TRN.PHENO1.glm.logistic.hybrid'), header=T, data.table=F)
    # COR[,anc]=p2cor(p = glm$P, n = gwasN[[anc]], sign=log(glm$OR))

    # training bootstrap data
    glm <- fread(paste0("/raid6/Tianyu/PRS/BootData/", anc,"_bootdata_repeat1_chr20"), header=T, data.table=F)
    ###### !!!! n needs to be able to change
    COR[,anc]=p2cor(p = glm$P, n = 5000, sign=log(glm$OR))
    
  }
  rownames(COR)=COR$ID
  
  # > head(COR)
  # CHR           ID         CEU         YRI
  # 5350016  20 20:80457:C:T  0.01006106 -0.01018116
  # 5350017  20 20:81154:T:G -0.01517010  0.01512512
  # 5350018  20 20:82590:T:G  0.01196787  0.01373924
  # 5350019  20 20:82603:A:C -0.01182462 -0.01355234
  # 5350020  20 20:83158:C:T  0.01517010 -0.01229402
  # 5350021  20 20:85259:G:A -0.01006106 -0.01048339
  
  # load the LD block information
  LDblocks=list()
  for(anc in gwasANC){
    if(anc == "CEU"){
      ld.anc="EUR"
    }else if(anc == "YRI"){
      ld.anc="AFR"
    }
    LDblocks[[anc]]=read.table2(system.file(paste0("data/Berisa.", 
                                                   paste0(ld.anc,".hg38"), ".bed"), 
                                            package="lassosum"), header=T)
  #  LDblocks[[anc]]$chr=gsub("chr","",LDblocks[[anc]]$chr)
  }
  
  # > head(LDblocks)
  # $CEU
  # chr     start      stop
  # 1     chr1   1961168   3666172
  # 2     chr1   3666172   4320751
  # 3     chr1   4320751   5853833
  # 4     chr1   5853833   7187275
  # 5     chr1   7187275   9305140
  # 6     chr1   9305140  10746927
  # 7     chr1  10746927  11717784
  # 8     chr1  11717784  12719464
  # 9     chr1  12719464  14565015
  
  
  # reference files by chromosome
  referenceFiles=list()
  for(chr in 1:22){
    referenceFiles[[chr]]=list()
    for(anc in gwasANC){
      # referenceFiles[[chr]][[anc]]=paste0("/data3/DownLoadedData/GWAS-Populations-SimulationInput/",anc,"_reference_LDblocks/CHR/",anc,"-chr",chr)
      referenceFiles[[chr]][[anc]] <- paste0(work.dir, "GWAS-Populations-SimulationInput/", anc, "_reference_LDblocks/CHR/", anc,"-chr",chr)
    }
  }
  
  # > referenceFiles
  # [[1]]
  # [[1]]$CEU
  # [1] "/raid6/Ron/prs/data/bert_sample/GWAS-Populations-SimulationInput/CEU_reference_LDblocks/CHR/CEU-chr1"
  # 
  # [[1]]$YRI
  # [1] "/raid6/Ron/prs/data/bert_sample/GWAS-Populations-SimulationInput/YRI_reference_LDblocks/CHR/YRI-chr1"
  # 
  # 
  # [[2]]
  # [[2]]$CEU
  # [1] "/raid6/Ron/prs/data/bert_sample/GWAS-Populations-SimulationInput/CEU_reference_LDblocks/CHR/CEU-chr2"
  # 
  # run lassosum
  print(CHR); flush.console()

  system.time(re.chr<-mclapply(CHR,mylassosumFunction,gamma=gamma,lambda=lambda,shrink=shrink,mem.limit = mem.limit,gwasANC = gwasANC,COR=COR,referenceFiles = referenceFiles,LDblocks = LDblocks,
                   mc.cores = mc.cores,mc.preschedule = F))
  print('loss in combined lassosum')
  print(re.chr[[1]]$loss)
  print(re.chr[[1]]$trainerror1)
  print(re.chr[[1]]$trainerror2)
  
  save(re.chr, file = '/raid6/Tianyu/PRS/trash/re_chr.RData')
  # merge the results from the 22 chromosomes
  re.lasso = merge.mylassosum(re.chr)
  print('loss in combined lassosum, after merging')
  print(re.lasso$loss)
  print(re.lasso$trainerror1)
  print(re.lasso$trainerror2)

  # save ther results
  # dir.create(paste0(main.dir,"CombinedLassoSum/Tmp/"),showWarnings = F,recursive = T)
  # when training the original data
  # save(re.lasso,file=paste(main.dir,"CombinedLassoSum/Tmp/GWAS-lasso-C",gwasN$CEU,"-Y",gwasN$YRI,sprintf("-gamma-%.2f",gamma),".Rdata",sep=""))
  # training the bootstrap data
  save(re.lasso,file=paste(main.dir,"CombinedLassoSum/Tmp/GWAS-lasso-C",gwasN$CEU,"-Y",gwasN$YRI,sprintf("-gamma-%.2f",gamma),"_boot1.Rdata",sep=""))
  return(i.combn)
}
# end of the wrapper function

system.time(re.wrapper<-mclapply(1:nrow(input.df),wrapperFunction,input.df=input.df,
                                 gwasANC = gwasANC, lambda=lambda, shrink=shrink,
                                 main.dir=main.dir,work.dir=work.dir,CHR= 20,mem.limit=2e10,
                                 mc.cores=16,mc.preschedule = F, mc.silent=F))
#####print the above variables
# > input.df
# gamma    N1    N2
# 1   0.5 20000 20000
# > gwasANC
# [1] "CEU" "YRI"
# > lambda
# [1] 0.001000000 0.001429969 0.002044812 0.002924018 0.004181255 0.005979066
# [7] 0.008549880 0.012226064 0.017482895 0.025000000




