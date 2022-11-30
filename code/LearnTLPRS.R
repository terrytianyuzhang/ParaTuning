##plink-1.9 needs to be pre-installed.
# library(devtools)
# install_github("ZhangchenZhao/TLPRS")
#########Step 1. Run TL-PRS using example data#########
library(data.table)
library(lassosum)
library(TLPRS)
setwd("/raid6/Tianyu/PRS/LearnTLPRS/ExampleData_TL-PRS")  ###setup your work path.
ped_file="ped_validation.txt";
Covar_name="sex"
Y_name="Y"
Ytype="C"
train_file="ExampleData_1000G_African_train"
validate_file="ExampleData_1000G_African_validate"
sum_stats_file="Beta_lassosum_usingEuropeanSummaryStatistics.txt"
target_sumstats_file="African_SummaryStatistics.txt"
LDblocks="AFR.hg19"
outfile="Output_TLPRS"
system.time({out.beta=TL_PRS(ped_file,Covar_name,Y_name, Ytype="C",train_file,validate_file,sum_stats_file,target_sumstats_file, LDblocks,outfile,cluster=NULL)})
summary(out.beta)
write.table(out.beta, file=,quote=F,row.names=F,col.names=F)

#########Step 2. Compare TL-PRS with its baseline method using example testing data##########
test_file="ExampleData_1000G_African_test"
ped_test_file="ped_test.txt"

Correlation_betweenPRSandY <- function( test_file, stats_file, ped_test_file,Y_name){
  Y_test=read.table(ped_test_file,header=T)
  cmd=paste0("plink-1.9 --bfile ",test_file,"  --score ",stats_file, " sum", " --out ",stats_file,".test.PRS")
  system(cmd)
  temp=read.table(paste0(stats_file,".test.PRS.profile"),header=T)
  merged=merge(temp,Y_test,by=c("FID","IID"))
  return(cor(merged[,Y_name],merged$SCORESUM))
}

###Correlation between lassosum PRS and Y in the testing data
Correlation_betweenPRSandY( test_file,sum_stats_file, ped_test_file, Y_name)  

###Correlation between TL-PRS and Y in the testing data
Correlation_betweenPRSandY( test_file, paste0(outfile, "_best.beta.txt"), ped_test_file,Y_name)
