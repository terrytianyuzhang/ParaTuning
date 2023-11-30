gc()
options(stringsAsFactors = F)

if(!exists("i.sim")){
  i.sim <- 609
}
#### load the functions that are needed
source("simulation-functions.R")
source("general_pipeline_parameters.R")
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

source('fit_JLS_function.R')

# setting.title <- 'CEU1aYRI2a22Chr'
##### general setup
# i.sim="NC8000"
# load(paste0("Work/Sim-C0.80-Y0.60-rep",i.sim,"/simulation-params.RData"))
# load("/raid6/Ron/prs/data/bert_sample/simulation-params.RData")

# set directories for this simulation
# main.dir=params$run.info$main.dir
# work.dir=params$run.info$work.dir

# load the parameters for this simulation
# load(paste0("Work/Sim-",i.sim,"/simulation-params.RData"))
# main.dir=params$run.info$main.dir #"/raid6/Tianyu/PRS/SimulationPipeline/"
# work.dir=params$run.info$work.dir #"/raid6/Tianyu/PRS/SimulationPipeline/Work/Sim-800/"

### ancestries
gwasANC=c("CEU","YRI")
# #### set the parameters
# GAMMA = c(0.2, 0.5, 0.8)
# ###!!!!!
# lambda=exp(seq(log(0.0025), log(0.025), length.out=10))
# lambda[1] <- 0.00025
# lambda=exp(seq(log(0.007), log(0.05), length.out=10))  ##try larger lambda

# lambda=exp(seq(log(0.001), log(0.025), length.out=10))
# shrink=.9

####sample sizes
# N1 <- params$CEU.TRN$n.case + params$CEU.TRN$n.control
# N2 <- params$YRI.TRN$n.case + params$YRI.TRN$n.control
N1 <- sample_sizes$CEU$n.case + sample_sizes$CEU$n.control
N2 <- sample_sizes$YRI$n.case + sample_sizes$YRI$n.control
input.df <- data.frame(gamma = GAMMA,
                       N1 = rep(N1, length(GAMMA)),
                       N2 = rep(N2, length(GAMMA)))
###!!!!this is not correct
# N=seq(20000,20000,4000)
# input.df=data.frame(N1=rep(N,each=length(N)),N2=rep(N,length(N)))
# input.df=data.frame(gamma=rep(GAMMA,each=nrow(input.df)),N1=rep(input.df$N1,length(GAMMA)),N2=rep(input.df$N2,length(GAMMA)))
# 
# ###this is the correct sample size
# input.df <- data.frame(gamma = GAMMA, N1 = 20000, N2 = 4000)
# memory limit
mem.limit <- 2*10e9

# end of the wrapper function

system.time(re.wrapper<-mclapply(1:nrow(input.df),wrapperFunction,input.df=input.df,
                                 gwasANC = gwasANC, lambda=lambda, shrink=shrink,
                                 main.dir=main.dir,work.dir=work.dir,CHR= 1:22, mem.limit=2e10,
                                 mc.cores=3, mc.preschedule = F, mc.silent=F))

rm.list=ls()
rm.list=rm.list[!(rm.list %in% c("i.sim","sims", "plink", "plink2"))]
rm(list = rm.list); flush.console()


