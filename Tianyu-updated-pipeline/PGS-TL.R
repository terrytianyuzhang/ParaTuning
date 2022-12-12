gc()
options(stringsAsFactors = F)

library(data.table)
library(lassosum)
library(TLPRS)

if(!exists("plink") | !exists("plink2")){
  plink <- "/usr/local/bin/plink"
  plink2 <- "/usr/local/bin/plink2"
}

if(!exists("i.sim")){
  i.sim <- 800
}

# load the parameters for this simulation
load(paste0("Work/Sim-",i.sim,"/simulation-params.RData"))
main.dir=params$run.info$main.dir #"/raid6/Tianyu/PRS/SimulationPipeline/"
work.dir=params$run.info$work.dir #"/raid6/Tianyu/PRS/SimulationPipeline/Work/Sim-800/"

#########prepare the files in the required format######
####target population validation set outcome
val_outcome <- fread(paste0(work.dir, '/YRI.TUNE/YRI.TUNE.fam'), data.table = F,
                     header = F)
val_outcome <- val_outcome[,c(1,2,5,6)] ##FID,IID,SEX and Y
colnames(val_outcome) <- c('FID' , 'IID', 'sex', 'Y')
fwrite(val_outcome, paste0(work.dir, '/YRI.TUNE/YRI-TUNE-fam-forTL.txt'))

####target population summary stat
target.gwas <- fread(paste0(work.dir, "/YRI.TRN/Assoc/YRI.TRN.PHENO1.glm.logistic.hybrid"), data.table = F, header = T)
target.gwas.print <- data.frame(SNP = target.gwas$ID,
                                A1 = target.gwas$A1,
                                beta = log(target.gwas$OR),
                                N = target.gwas$OBS_CT,
                                p = target.gwas$P)
fwrite(target.gwas.print, paste0(work.dir, '/YRI.TRN/YRI-TRN-GWAS-forTL.txt'))

#####
READ AUC.weighted
lassosum.auc <- get(load(paste0(work.dir,"/YRI.TST/Lassosum/YRI.TST/CEU.TRN-YRI.TST-lassosum.AUC.RData")))
FIND THE BEST INDEX
best.auc <- which.max(lassosum.auc$auc)
READ lassosum
full.df.beta <- paste0(test.dir,"Lassosum/",train.set,"-lassosum-",test.set,".RData")
best.beta <- full.df.beta$beta$`0.9`[,best.auc]
READ MAP
COMBINE SNP ID AND beta(
SAVE




#########Step 1. Run TL-PRS using example data#########
# setwd("/raid6/Tianyu/PRS/LearnTLPRS/ExampleData_TL-PRS")  ###setup your work path.
# ped_file="ped_validation.txt" ##need to create this?
ped_file <- paste0(work.dir, '/YRI.TUNE/YRI-TUNE-fam-forTL.txt')
Y_name <- "Y" ##what is Y name
Ytype <- "B" ##binary outcome

# train_file="ExampleData_1000G_African_train"
train_file <- paste0(work.dir, '/YRI.TRN/YRI.TRN')
# validate_file="ExampleData_1000G_African_validate"
validate_file <- paste0(work.dir, '/YRI.TUNE/YRI.TUNE')
sum_stats_file="Beta_lassosum_usingEuropeanSummaryStatistics.txt"
# target_sumstats_file="African_SummaryStatistics.txt"
target_sumstats_file <- paste0(work.dir, '/YRI.TRN/YRI-TRN-GWAS-forTL.txt')

LDblocks <- "AFR.hg19"
outfile="Output_TLPRS"

system.time(out.beta <- TL_PRS(ped_file,Covar_name,Y_name, 
                               Ytype= Ytype, 
                               train_file = train_file,
                               validate_file = validate_file,
                               sum_stats_file = sum_stats_file,
                               target_sumstats_file = target_sumstats_file,
                               LDblocks = LDblocks,
                               outfile = outfile,
                               cluster=NULL,
                               plink = plink))
summary(out.beta)
write.table(out.beta, file='out.beta.txt',quote=F,row.names=F,col.names=F)

#########Step 2. Compare TL-PRS with its baseline method using example testing data##########
test_file="ExampleData_1000G_African_test"
ped_test_file="ped_test.txt"

Correlation_betweenPRSandY <- function( test_file, stats_file, ped_test_file,Y_name){
  Y_test=read.table(ped_test_file,header=T)
  cmd=paste0("/usr/local/bin/plink --bfile ",test_file,"  --score ",stats_file, " sum", " --out ",stats_file,".test.PRS")
  system(cmd)
  temp=read.table(paste0(stats_file,".test.PRS.profile"),header=T)
  merged=merge(temp,Y_test,by=c("FID","IID"))
  return(cor(merged[,Y_name],merged$SCORESUM))
}

###Correlation between lassosum PRS and Y in the testing data
Correlation_betweenPRSandY( test_file,sum_stats_file, ped_test_file, Y_name)  

###Correlation between TL-PRS and Y in the testing data
Correlation_betweenPRSandY( test_file, paste0(outfile, "_best.beta.txt"), ped_test_file,Y_name)

rm.list=ls()
rm.list=rm.list[!(rm.list %in% c("i.sim","sims", "plink", "plink2"))]
rm(list = rm.list); flush.console()