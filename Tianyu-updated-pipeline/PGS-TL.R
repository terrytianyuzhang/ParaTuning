gc()
# library(devtools)
# install_github("terrytianyuzhang/TLPRS")

options(stringsAsFactors = F)

library(data.table)
library(lassosum)
library(TLPRS)
library(parallel)
library(pROC)

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

###folder that contains all the data formulated in the required format
TL.data.dir <- paste0(work.dir, 'TLPGS-data/')
dir.create(TL.data.dir,recursive = T,showWarnings = F)

if(1){
#########prepare the files in the required format######
####target population validation set outcome
val_outcome <- fread(paste0(work.dir, '/YRI.TUNE/YRI.TUNE.fam'), data.table = F,
                     header = F)
val_outcome <- val_outcome[,c(1,2,5,6)] ##FID,IID,SEX and Y
colnames(val_outcome) <- c('FID' , 'IID', 'sex', 'Y')
val_outcome$Y <- val_outcome$Y - 1 ##TL requires the outcome to be between 0 and 1
fwrite(val_outcome, paste0(TL.data.dir, '/YRI-TUNE-fam-forTL.txt'), sep = ' ')

####target population test set outcome
test_outcome <- fread(paste0(work.dir, '/YRI.TST/YRI.TST.fam'), data.table = F,
                      header = F)
test_outcome <- test_outcome[,c(1,2,5,6)] ##FID,IID,SEX and Y
colnames(test_outcome) <- c('FID' , 'IID', 'sex', 'Y')
test_outcome$Y <- test_outcome$Y - 1 ##TL requires the outcome to be between 0 and 1
fwrite(test_outcome, paste0(TL.data.dir, '/YRI-TST-fam-forTL.txt'), sep = ' ')

####target population summary stat
target.gwas <- fread(paste0(work.dir, "/YRI.TRN/Assoc/YRI.TRN.PHENO1.glm.logistic.hybrid"), data.table = F, header = T)
target.gwas.print <- data.frame(SNP = target.gwas$ID,
                                A1 = target.gwas$A1,
                                beta = log(target.gwas$OR),
                                N = target.gwas$OBS_CT,
                                p = target.gwas$P)
fwrite(target.gwas.print, paste0(TL.data.dir, '/YRI-TRN-GWAS-forTL.txt'))


#####
#READ AUC.weighted
lassosum.auc <- get(load(paste0(work.dir,"/YRI.TST/Lassosum/YRI.TST/CEU.TRN-YRI.TST-lassosum.AUC.RData")))

#FIND THE BEST INDEX
best.auc <- which.max(lassosum.auc$auc)
#READ lassosum
full.df.beta <- get(load(paste0(work.dir,"/YRI.TST/Lassosum/CEU.TRN-lassosum-YRI.TST.RData")))
best.beta <- full.df.beta$beta$`0.9`[,best.auc]
# head(full.df.beta$sumstats)
# chr     pos A1 A2           cor order
# 1   1 1962845  C  T  0.0008701268     1
# 2   1 1962899  C  A -0.0040755079     2
# 3   1 1963406  A  G  0.0021235751     3
# 4   1 1963538  C  T -0.0002843339     4
# 5   1 1963738  T  C  0.0021235751     5
# 6   1 1964101  G  A  0.0021235751     6

#READ MAP
map <- fread(paste0(work.dir, 'YRI.TRN/YRI.TRN.pvar'), header=T, data.table=T)
names(map)[3] <- 'SNP'
# head(map)
# #CHROM     POS            ID REF ALT
# 1      1 1962845 1:1962845:T:C   T   C
# 2      1 1962899 1:1962899:A:C   A   C
# 3      1 1963406 1:1963406:G:A   G   A
# 4      1 1963538 1:1963538:T:C   T   C
# 5      1 1963738 1:1963738:C:T   C   T
# 6      1 1964101 1:1964101:A:G   A   G

####best.beta length 4734473, but map has 5630745. Joint lassosum's beta has 5630745.
#COMBINE SNP ID AND beta
sumstat.beta <- cbind(full.df.beta$sumstats, matrix(best.beta, ncol = 1))
names(sumstat.beta)[7] <- 'Beta'
sumstat.beta <- data.table(sumstat.beta)
sumstat.beta[,chr := as.numeric(chr)]

sumstat.beta.id <- merge(sumstat.beta,
                         map,
                         by.x = c('chr','pos'),
                         by.y = c('#CHROM','POS'))

# > sumstat.beta.id
# chr      pos A1 A2           cor order        Beta              ID REF
# 1:   1  1962845  C  T  0.0008701268     1 0.000000000   1:1962845:T:C   T
# 2:   1  1962899  C  A -0.0040755079     2 0.000000000   1:1962899:A:C   A
# 3:   1  1963406  A  G  0.0021235751     3 0.000000000   1:1963406:G:A   G
# 4:   1  1963538  C  T -0.0002843339     4 0.000000000   1:1963538:T:C   T
# 5:   1  1963738  T  C  0.0021235751     5 0.000000000   1:1963738:C:T   C
# ---
#   4734469:  22 50797551  A  G  0.0147442549 74344 0.001745122 22:50797551:G:A   G
# 4734470:  22 50798021  G  A -0.0028225845 74345 0.000000000 22:50798021:A:G   A
# 4734471:  22 50798635  C  T  0.0054757089 74346 0.000000000 22:50798635:T:C   T
# 4734472:  22 50801260  C  T  0.0006475620 74347 0.000000000 22:50801260:T:C   T
# 4734473:  22 50802392  T  C  0.0021817060 74348 0.000000000 22:50802392:C:T   C

sumstat.beta.id.simple <- sumstat.beta.id[, .(SNP, A1, Beta)]
fwrite(sumstat.beta.id.simple,
       paste0(TL.data.dir, '/CEU-TRN-Lassosum-forTL.txt'), sep = ' ')

# fwrite(sumstat.beta.id.simple,
#        paste0(work.dir, '/CEU.TRN/CEU-TRN-Lassosum-forTL.txt'))

# target.gwas.print <- fread(paste0(work.dir, '/YRI.TRN/YRI-TRN-GWAS-forTL.txt'),
#                            header = T)
# two.txt <- merge(sumstat.beta.id.simple, target.gwas.print,
#                  by = 'SNP')
# > identical(two.txt$A1.x, two.txt$A1.y)
# [1] TRUE
}

#########Step 1. Run TL-PRS using example data#########
# setwd("/raid6/Tianyu/PRS/LearnTLPRS/ExampleData_TL-PRS")  ###setup your work path.
# ped_file="ped_validation.txt" ##need to create this?
ped_file <- paste0(TL.data.dir, 'YRI-TUNE-fam-forTL.txt')
Y_name <- "Y" ##what is Y name
Ytype <- "B" 

# train_file="ExampleData_1000G_African_train"
train_file <- paste0(work.dir, 'YRI.TRN/YRI.TRN')
# validate_file="ExampleData_1000G_African_validate"
validate_file <- paste0(work.dir, 'YRI.TUNE/YRI.TUNE')
# sum_stats_file="Beta_lassosum_usingEuropeanSummaryStatistics.txt"
sum_stats_file <- paste0(TL.data.dir, 'CEU-TRN-Lassosum-forTL.txt')
# target_sumstats_file="African_SummaryStatistics.txt"
target_sumstats_file <- paste0(TL.data.dir, 'YRI-TRN-GWAS-forTL.txt')

LDblocks <- "AFR.hg19" ##it used to be "AFR.hg19"
# outfile="Output_TLPRS"
outfile <- paste0(TL.data.dir, 'TL')

cl <- makeCluster(24, type="FORK")

system.time(out.beta <- TL_PRS(ped_file = ped_file,
                               Covar_name = NULL,
                               Y_name = Y_name, 
                               Ytype= Ytype, 
                               train_file = train_file,
                               test_file = validate_file,
                               sum_stats_file = sum_stats_file,
                               target_sumstats_file = target_sumstats_file,
                               LDblocks = LDblocks,
                               outfile = outfile,
                               cluster= cl,
                               plink = plink))
stopCluster(cl) 

AUC.from.beta.candidates <- function(beta.can.file,
                                     test_file, stats_file, ped_test_file,Y_name){
  beta.all <- fread(beta.can.file) ##read in all the beta
  
  ##normalize and make sure the signs are correct
  beta.all[,11:ncol(beta.all)] <- beta.all[,11:ncol(beta.all)] / beta.all$sd
  not.same.A1 <- which(beta.all$V5 != beta.all$A1)
  beta.all[not.same.A1,11:ncol(beta.all)] <- -beta.all[not.same.A1,11:ncol(beta.all)]
  names(beta.all)[1:2] <- c("SNP","A1")
  
  Y_test=read.table(ped_test_file,header=T, sep = ' ')
  
  beta.all <- as.data.frame(beta.all)
  aucs <- rep(0, length(12:ncol(beta.all)))
  k <- 1
  for(i in 12:ncol(beta.all)){
    temp.beta <- beta.all[, c(1,2,i)]
    fwrite(temp.beta,
           file = stats_file, sep = ' ',
           col.names = FALSE)
    cmd <- paste0("/usr/local/bin/plink --bfile ",test_file,
                  "  --score ",stats_file, " sum",
                  " --out ",stats_file,".test.PRS",
                  " --allow-no-sex")
    system(cmd)
    print("0")
    temp=read.table(paste0(stats_file,".test.PRS.profile"),header=T)
    print("2")
    merged=merge(temp,Y_test,by=c("FID","IID"))
    print("3")
    aucs[k] <- auc(merged[,Y_name],merged$SCORESUM)
    k <- k+1
    print(aucs)
  }
  
  return(aucs)
  # > aucs
  # [1] 0.7124512 0.7124523 0.7124569 0.7124571 0.7124580 0.7124589 0.7124601
  # [8] 0.7124647 0.7124666 0.7124677 0.7124701 0.7124725 0.7124764 0.7124796
  # [15] 0.7124797 0.7124677 0.7124936 0.7125233 0.7125501 0.7125733 0.7125966
  # [22] 0.7126212 0.7126488 0.7126644 0.7126870 0.7127111 0.7127335 0.7127590
  # [29] 0.7127785 0.7128108 0.7126876 0.7129319 0.7131436 0.7133729 0.7136186
  # [36] 0.7138386 0.7140674 0.7142802 0.7145073 0.7147295 0.7149423 0.7151689
  # [43] 0.7153975 0.7155980 0.7158003 0.7147240 0.7168759 0.7187756 0.7206672
  # [50] 0.7223919 0.7239199 0.7254027 0.7266540 0.7278350 0.7288401 0.7297561
  # [57] 0.7305475 0.7312015 0.7318841 0.7323994
}
AUC.best.beta <- function(beta.can.file,
                          test_file, stats_file, ped_test_file,Y_name,
                          best.beta.index){
  beta.all <- fread(beta.can.file) ##read in all the beta
  
  ##normalize and make sure the signs are correct
  beta.all[,11:ncol(beta.all)] <- beta.all[,11:ncol(beta.all)] / beta.all$sd
  not.same.A1 <- which(beta.all$V5 != beta.all$A1)
  beta.all[not.same.A1,11:ncol(beta.all)] <- -beta.all[not.same.A1,11:ncol(beta.all)]
  names(beta.all)[1:2] <- c("SNP","A1")
  
  Y_test=read.table(ped_test_file,header=T, sep = ' ')
  
  beta.all <- as.data.frame(beta.all)
  best.auc <- -1
  
  k <- 1
  for(i in 12:ncol(beta.all)){
    if(k == best.beta.index){
      temp.beta <- beta.all[, c(1,2,i)]
      fwrite(temp.beta,
             file = stats_file, sep = ' ',
             col.names = FALSE)
      cmd <- paste0("/usr/local/bin/plink --bfile ",test_file,
                    "  --score ",stats_file, " sum",
                    " --out ",stats_file,".test.PRS",
                    " --allow-no-sex")
      system(cmd)

      temp=read.table(paste0(stats_file,".test.PRS.profile"),header=T)
      
      merged=merge(temp,Y_test,by=c("FID","IID"))
      
      best.auc <- auc(merged[,Y_name],merged$SCORESUM)
      print(best.auc)
    }##if best beta
    k <- k+1
  }
  
  return(best.auc)
}


test_file <- paste0(work.dir, 'YRI.TUNE/YRI.TUNE')
ped_test_file <- paste0(TL.data.dir, 'YRI-TUNE-fam-forTL.txt')
Y_name <- "Y"
stats_file <- paste0(TL.data.dir,'TL_temp.beta.txt')
beta.can.file <- paste0(TL.data.dir,"TL_beta.candidates.fresh.txt")

tune.auc <- AUC.from.beta.candidates(beta.can.file,
                                     test_file, stats_file, ped_test_file, Y_name)
save(tune.auc, file = paste0(TL.data.dir,'TL_tune-auc.RData'))

# tune.auc <- get(load(paste0(main.dir,'firstTL_tune-auc.RData')))

test_file <- paste0(work.dir, 'YRI.TST/YRI.TST')
ped_test_file <- paste0(TL.data.dir, 'YRI-TST-fam-forTL.txt')
Y_name <- "Y"
stats_file <- paste0(TL.data.dir,'TL_temp.beta.txt')
beta.can.file <- paste0(TL.data.dir,"TL_beta.candidates.fresh.txt")

test.auc <- AUC.best.beta(beta.can.file,
                          test_file, stats_file, ped_test_file, Y_name,
                          which.max(tune.auc))
save(test.auc, file = paste0(TL.data.dir,'TL_test-auc.RData'))
# > test.aucs
# [1] 0.7124512 0.7124523 0.7124569 0.7124571 0.7124580 0.7124589 0.7124601
# [8] 0.7124647 0.7124666 0.7124677 0.7124701 0.7124725 0.7124764 0.7124796
# [15] 0.7124797 0.7147240 0.7168759 0.7187756 0.7206672 0.7223919 0.7239199
# [22] 0.7254027 0.7266540 0.7278350 0.7288401 0.7297561 0.7305475 0.7312015
# [29] 0.7318841 0.7323994 0.7289631 0.7335943 0.7304489 0.7235622 0.7153189
# [36] 0.7070783 0.6993625 0.6923570 0.6860880 0.6805120 0.6755146 0.6710677
# [43] 0.6670146 0.6633822 0.6601632 0.6785349 0.6471623 0.6324021 0.6242296
# [50] 0.6185533 0.6144268 0.6112639 0.6086980 0.6064754 0.6046272 0.6030144
# [57] 0.6015476 0.6002817 0.5991710 0.5981460


# 
# 
# summary(out.beta)
# save(out.beta, file = 'TL.first.out.RData')
# write.table(out.beta, file='out.beta.txt',quote=F,row.names=F,col.names=F)
# 
# 
# 
# all.beta <- fread(file = '/raid6/Tianyu/PRS/SimulationPipeline/firstTL_beta.candidates.txt')
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #########Step 2. Compare TL-PRS with its baseline method using example testing data##########
# 
# 
# # test_file="ExampleData_1000G_African_test"
# test_file <- paste0(work.dir, 'YRI.TST/YRI.TST')
# # ped_test_file="ped_test.txt"
# ped_test_file <- paste0(work.dir, '/YRI.TST/YRI-TST-fam-forTL.txt')
# 
# # Y_name <- "Y"
# Y_name <- "Y"
# 
# 
# outfile <- paste0(main.dir, 'firstTL')
# 
# 
# Correlation_betweenPRSandY <- function( test_file, stats_file, ped_test_file,Y_name){
#   Y_test=read.table(ped_test_file,header=T, sep = ',')
#   cmd=paste0("/usr/local/bin/plink --bfile ",test_file,"  --score ",stats_file, " sum", " --out ",stats_file,".test.PRS",
#              " --allow-no-sex")
#   system(cmd)
#   temp=read.table(paste0(stats_file,".test.PRS.profile"),header=T)
#   merged=merge(temp,Y_test,by=c("FID","IID"))
#   return(auc(merged[,Y_name],merged$SCORESUM))
#   # return(cor(merged[,Y_name],merged$SCORESUM))
# }
# 
# ###Correlation between lassosum PRS and Y in the testing data
# Correlation_betweenPRSandY( test_file,
#                             stats_file = sum_stats_file, ped_test_file, Y_name)  
# ###0.71
# ###Correlation between TL-PRS and Y in the testing data
# Correlation_betweenPRSandY( test_file, 
#                             stats_file = paste0(outfile, "_best.beta.txt"), ped_test_file,Y_name)
# ###0.6189
# 
# ####check the performance of this beta of validation set
# Correlation_betweenPRSandY( validate_file, 
#                             stats_file = paste0(outfile, "_best.beta.txt"), ped_file,Y_name)
# ###0.6122
# 
# ####check if the best PGS is working
# ##read in PGS
# TL.best.val.PGS <- fread(paste0(main.dir, '/firstTL_best.PRS.txt'))
# ##read in outcome
# TL.val.outcome <- fread(paste0(work.dir, '/YRI.TUNE/YRI-TUNE-fam-forTL.txt'))
# ##calculate AUC
# auc(TL.val.outcome$Y,TL.best.val.PGS$PRS.NULL) ##0.604
# auc(TL.val.outcome$Y,TL.best.val.PGS$PRS.TL) ##0.6122


rm.list=ls()
rm.list=rm.list[!(rm.list %in% c("i.sim","sims", "plink", "plink2"))]
rm(list = rm.list); flush.console()