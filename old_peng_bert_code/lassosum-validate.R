rm(list=ls()); gc()
options(stringsAsFactors = F)


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

#### load the parameter information
mem.limit=5*20e10
i.sim=1
load("Work/Sim-C0.80-Y0.60-rep1/simulation-params.RData")

# set directories for this simulation
main.dir=params$run.info$main.dir
work.dir=params$run.info$work.dir


### load one of the lassosum runs
load("Work/Sim-C0.80-Y0.60-rep1/CombinedLassoSum/Tmp/GWAS-lasso-C12000-Y12000-gamma-0.20.Rdata")

## next step is validation
system.time(re.validate <- validate(re.lasso, "Work/Sim-C0.80-Y0.60-rep1/CEU.TRN/CEU.TRN-4000", extract=NULL, keep=NULL, distance=NULL, mem.limit=mem.limit))

