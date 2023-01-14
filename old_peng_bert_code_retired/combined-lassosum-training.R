rm(list=ls()); 
gc()
options(stringsAsFactors = F)


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
plink="/data3/Software/Plink/plink"
plink2="/data3/Software/Plink2/plink2"

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
i.sim="NC8000"
load(paste0("Work/Sim-C0.80-Y0.60-rep",i.sim,"/simulation-params.RData"))

# set directories for this simulation
main.dir=params$run.info$main.dir
work.dir=params$run.info$work.dir

### ancestries
gwasANC=c("CEU","YRI")
#### set the parameters
GAMMA=0.5
lambda=exp(seq(log(0.001), log(0.025), length.out=10))
shrink=.9

N=seq(20000,20000,4000)
input.df=data.frame(N1=rep(N,each=length(N)),N2=rep(N,length(N)))
input.df=data.frame(gamma=rep(GAMMA,each=nrow(input.df)),N1=rep(input.df$N1,length(GAMMA)),N2=rep(input.df$N2,length(GAMMA)))

# memory limit
mem.limit=2*10e9

wrapperFunction <- function(i.combn, input.df, gwasANC, lambda, shrink, main.dir, work.dir, CHR=1:22, mem.limit){

  # collect gamma
  gamma=input.df$gamma[i.combn]  

  # collectthe number of samples for the two populations
  gwasN=list()
  gwasN[["CEU"]]=input.df$N1[i.combn]
  gwasN[["YRI"]]=input.df$N2[i.combn]


  #### determine correlations
  ### Summary statistics ###
  # read a map for the first ancestry to set up a dataframe
  anc=gwasANC[1]
  map=fread(paste0(work.dir,anc,".GWAS/",anc,".GWAS-",gwasN[[anc]],".pvar"),header=T,data.table=F)
  COR=data.frame(CHR=map$`#CHROM`,ID=map$ID)
  for(anc in gwasANC){
    glm=fread(paste0(work.dir,anc,".GWAS/Assoc/",anc,".GWAS-",gwasN[[anc]],".PHENO1.glm.logistic.hybrid"),header=T,data.table=F)
    COR[,anc]=p2cor(p = glm$P, n = gwasN[[anc]], sign=log(glm$OR))
  }
  rownames(COR)=COR$ID

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
  
  # reference files by chromosome
  referenceFiles=list()
  for(chr in 1:22){
    referenceFiles[[chr]]=list()
    for(anc in gwasANC){
      referenceFiles[[chr]][[anc]]=paste0("/data3/DownLoadedData/GWAS-Populations-SimulationInput/",anc,"_reference_LDblocks/CHR/",anc,"-chr",chr)
    }
  }

  # run lassosum
  print(CHR); flush.console()

  system.time(re.chr<-mclapply(CHR,mylassosumFunction,gamma=gamma,lambda=lambda,shrink=shrink,mem.limit = mem.limit,gwasANC = gwasANC,COR=COR,referenceFiles = referenceFiles,LDblocks = LDblocks,
                   mc.cores = 12,mc.preschedule = F))
  # merge the results from the 22 chromosomes
  re.lasso=merge.mylassosum(re.chr)

  # save ther results
  dir.create(paste0(work.dir,"CombinedLassoSum/Tmp"),showWarnings = F,recursive = T)
  save(re.lasso,file=paste(work.dir,"CombinedLassoSum/Tmp/GWAS-lasso-C",gwasN$CEU,"-Y",gwasN$YRI,sprintf("-gamma-%.2f",gamma),".Rdata",sep=""))
  return(i.combn)
}
# end of the wrapper function

system.time(re.wrapper<-mclapply(1:nrow(input.df),wrapperFunction,input.df=input.df,gwasANC = gwasANC, lambda=lambda, shrink=shrink,main.dir=main.dir,work.dir=work.dir,CHR=1:22,mem.limit=2e10,
                                 mc.cores=5,mc.preschedule = F, mc.silent=F))
     
########### stop her for now


