#####load correlation vectors
library(data.table)
library(lassosum)

setting.title <- '1in2signal'
 ##suppose this is the population we are interested in

anc <- 'YRI'
anc <- 'CEU'
#####load estimated beta
re.lasso <- get(load("/raid6/Tianyu/PRS/CombinedLassoSum/Tmp/GWAS-lasso-C20000-Y4000-gamma-0.50.Rdata"))
trainerror1 <- re.lasso$trainerror1
trainerror2 <- re.lasso$trainerror2

#####load bootstrap Y
booty <- get(load(file = paste0('/raid6/Tianyu/PRS/BootData/',anc,'_bootY_',setting.title,'.RData')))

#####calculate bootstrap score
boot.lasso <- get(load(paste0("/raid6/Tianyu/PRS/CombinedLassoSum/Tmp/GWAS-lasso-C20000-Y20000-gamma-0.50_boot_",setting.title,"_",b,".Rdata")))
boot.beta <- boot.lasso$beta

####check the norms of betas
sqrt(colSums(re.lasso$beta^2))
sqrt(colSums(boot.beta^2))
sqrt(colSums(re.lasso$beta^2))/sqrt(colSums(boot.beta^2)) * 2*cov 
sqrt(colSums(re.lasso$beta^2))/sqrt(colSums(boot.beta^2)) * 2*cov  + trainerror1
#####loop through all chromosomes
#####load the genotype of one chromosome

#####calculate the predicted risk score 



B <- 1
cov <- rep(0, length(re.lasso$lambda))
cov <- rep(0, 1)
for(b in 1:B){
  boot.lasso <- get(load(paste0("/raid6/Tianyu/PRS/CombinedLassoSum/Tmp/GWAS-lasso-C20000-Y20000-gamma-0.50_boot_",setting.title,"_",b,".Rdata")))
  boot.beta <- boot.lasso$beta
  
  ##read correlation
  boot.summary <- fread(paste0("/raid6/Tianyu/PRS/BootData/",anc,"_bootdata_",setting.title,"_repeat", b),
                        header = T, data.table =F)
  
  boot.cor <- p2cor(p = boot.summary$P, n = boot.summary$n[1],sign = log(boot.summary$OR))

  ###estimated covariance
  cov.temp <- t(crossprod(boot.beta, boot.cor))
  print(cov.temp)
  # cov.temp <- sapply(boot.beta, cor, y = boot.cor)
  cov <- cov + cov.temp
}

cov <- cov/B

trainerror1 + 2*cov ##selection criteria

######BELOW CAN BE USED TO LOOK AT THE BOOT P-VALUES###
lassoorg <- get(load("/raid6/Tianyu/PRS/CombinedLassoSum/Tmp/GWAS-lasso-C20000-Y20000-gamma-0.50.Rdata"))
lassoboot <- get(load("/raid6/Tianyu/PRS/CombinedLassoSum/Tmp/GWAS-lasso-C20000-Y20000-gamma-0.50_boot_lowsignal_1.Rdata"))
#ceu
ceubooty <- get(load('/raid6/Tianyu/PRS/BootData/CEU_bootY_repeat20000'))
ceuboot <- fread('/raid6/Tianyu/PRS/BootData/CEU_bootdata_repeat1')
ceuorg <- fread('/raid6/Ron/prs/data/bert_sample/CEU.TRN.PHENO1.glm.logistic.hybrid')

names(ceuorg)[1] <- 'chr'
ceuorg <- ceuorg[chr %in% c(20,21),]
summary(ceucor.org <- p2cor(p = ceuorg$P, n = 20000,sign = log(ceuorg$OR)))
summary(ceucor.boot <- p2cor(p = ceuboot$P, n = 20000,sign = log(ceuboot$OR)))
crossprod(lassoboot$beta[,1], ceucor.boot)
crossprod(lassoorg$beta[,1], ceucor.org)
#yri
yribooty <- get(load('/raid6/Tianyu/PRS/BootData/YRI_bootY_repeat4000'))

yriboot <- fread('/raid6/Tianyu/PRS/BootData/YRI_bootdata_nosignal_repeat1')
yriorg <- fread('/raid6/Ron/prs/data/bert_sample/YRI.TRN.PHENO1.glm.logistic.hybrid')

names(yriorg)[1] <- 'chr'
yriorg <- yriorg[chr %in% c(20,21),]
summary(yricor.org <- p2cor(p = yriorg$P, n = 4000,sign = log(yriorg$OR)))
summary(yricor.boot <- p2cor(p = yriboot$P, n = 4000,sign = log(yriboot$OR)))
# 


##########
