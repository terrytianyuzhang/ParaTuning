#####there are 129,897 SNPs in chr20
#####dim(numeric_ref_genotypes)
#####[1]   5000 129897

#####create a fake \hat \beta_0
set.seed(1)
chr20snpN <- 129897
beta0 <- rep(0, chr20snpN)
non_zero_snp <- sample(1:chr20snpN, 100)
beta0[non_zero_snp] <- 1

######use a loop to calculate tau and C\beta #####
# install.packages('RSpectra')
library(RSpectra)
library(lassosum)
library(data.table)
library(snpStats)
setwd('/home/tianyuz3/PRS/my_code/')

chr_index <- 'chr20'
ld.anc <- 'AFR' ###ancestor
sample_size <- 5000
if(ld.anc == 'AFR'){
  anc_initial <- 'YRI'
}

SVD_big_list <- get(load(paste0("/raid6/Tianyu/PRS/SVDdata/", anc_initial, "reference_LDblocks_", chr_index, "_SVD.RData")))

#nothing more than doing a matrix, vector multiplication
pointer_begin <- pointer_end <- 1
Cbeta <- NULL
for(i in 1:length(SVD_big_list)){
  tempC <- SVD_big_list[[i]]$C
  
  pointer_end <- pointer_begin + dim(tempC)[1] - 1
  Cbeta <- c(Cbeta, tempC %*% beta0[pointer_begin:pointer_end])
  
  pointer_begin <- pointer_end + 1
}

###these quantities will be used during boostrap
nCbeta <- sample_size * Cbeta
tau <- t(beta0) %*% as.matrix(Cbeta, ncol = 1)

####save files
pre_cal <- list(tau = tau, nCbeta = nCbeta)
save(pre_cal, file = paste0("/raid6/Tianyu/PRS/SVDdata/", anc_initial, "pre_cal_for_bootstrap.RData"))
