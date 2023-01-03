######Tianyu is the original creater of this file
######tianyuz3@andrew.cmu.edu####

####calculate the risk score of some simulated individuals
####the genotype information is simulated
####the regression coefficient vector beta0 is calcuated from the original fit of combined lasso
####it should be one corresponding to a smallish lambda

PGS_bychr_bootstrap <- function(chr, anc, beta0){
  print(paste0("calculate PGS with chr", chr))
  
  chr_loc <- as.numeric(gsub(':.*$','',names(beta0)))
  sub.beta0 <- beta0[chr_loc == chr & beta0 != 0] #subset relevant coordinates of beta
  snp <- names(sub.beta0)
  
  print(head(snp))
  print(head(sub.beta0))
  ####load genotype data without label
  if(anc == 'CEU'){
    file.title <- 'CEU-20K'
  }else{
    file.title <- 'YRI-4K'  
  }
  
  ###this is the reference panel
  system.time(gnt<-read.plink(bed=paste0("/raid6/Tianyu/PRS/bert_sample/ReferencePopulation-Package/", file.title,"/CHR/",file.title,"-chr", chr,".bed"),
                              bim=paste0("/raid6/Tianyu/PRS/bert_sample/ReferencePopulation-Package/", file.title,"/CHR/",file.title,"-chr", chr,".bim"),
                              fam=paste0("/raid6/Tianyu/PRS/bert_sample/ReferencePopulation-Package/", file.title,"/CHR/",file.title,"-chr", chr,".fam"),
                              select.snps = snp)
  )
  # system.time(gnt<-read.plink(bed=paste0("/raid6/Tianyu/PRS/bert_sample/",anc,".TUNE/CHR/",anc,".TUNE-chr",chr,".bed"),
  #                             bim=paste0("/raid6/Tianyu/PRS/bert_sample/",anc,".TUNE/CHR/",anc,".TUNE-chr",chr,".bim"),
  #                             fam=paste0("/raid6/Tianyu/PRS/bert_sample/",anc,".TUNE/CHR/",anc,".TUNE-chr",chr,".fam"),
  #                             select.snps = snp)
  # )
  #example : "/raid6/Tianyu/PRS/bert_sample/YRI.TUNE/CHR/YRI.TUNE-chr20.bed"
  
  # transform the genotypes to a matrix, and reverse the call count. This snpStats package counts the A2 allele, not the A1
  system.time(gnt<- 2 - as(gnt$genotypes,Class="numeric"))
  
  print(gnt[1:5,1:5])
  
  # center the columns
  system.time(gnt <-gnt - rep(1, nrow(gnt)) %*% t(colMeans(gnt)))
  # normalize the calls to the unit 1 norm
  system.time(gnt<-normalize.cols(gnt,method="euclidean",p=2))
  
  # calculate the pgs
  system.time(re.pgs<-gnt %*% sub.beta0)
  print(paste0("finish calculating PGS with chr", chr))
  
  return(re.pgs)
}

PGS_bychr_bootstrap_allbeta <- function(chr, anc, beta){
  print(paste0('claculating PGS of ',anc, 'using chr ', chr))
  
  ######next step is generating some risk score using reference genotype
  ######we can download these genotype from 1000 Genome Project
  ######for the purpose of thee paper I will just use the reference panel 
  ######from which the training samples are generated
  
  chr_loc <- as.numeric(gsub(':.*$','',rownames(beta)))
  snp <- rownames(beta)[chr_loc == chr]
  
  ####load genotype data without label
  
  if(anc == 'CEU'){
    file.title <- 'CEU-20K'
  }else{
    file.title <- 'YRI-4K'  
  }
  
  boost.data.folder <- paste0("/raid6/Tianyu/PRS/BootData/", file.title,"/CHR")
  boost.data.val.data <- paste0(boost.data.folder, "/", file.title,"-chr", chr, "boost-val")
  
  ###this is the reference panel
  system.time(gnt<-read.plink(bed=paste0(boost.data.val.data,".bed"),
                              bim=paste0(boost.data.val.data,".bim"),
                              fam=paste0(boost.data.val.data,".fam"),
                              select.snps = snp)
  )
  #example : "/raid6/Tianyu/PRS/bert_sample/YRI.TUNE/CHR/YRI.TUNE-chr20.bed"
  
  # transform the genotypes to a matrix, and reverse the call count. This snpStats package counts the A2 allele, not the A1
  system.time(gnt<- 2 - as(gnt$genotypes,Class="numeric"))
  
  # center the columns
  system.time(gnt <-gnt - rep(1, nrow(gnt)) %*% t(colMeans(gnt)))
  # normalize the calls to the unit 1 norm
  system.time(gnt<-normalize.cols(gnt,method="euclidean",p=2))
  
  # apply the proper shrinkage
  # system.time(gnt<-gnt*sqrt(1-shrink))
  
  # calculate the pgs
  system.time(re.pgs<-gnt%*%beta[chr_loc == chr,])
  print(paste0('finished claculating PGS of ',anc, 'using chr ', chr))
  
  return(re.pgs)
}

####generate boostrap Y
####mean value of Y should be similar to the case proportion
####Y is assumed to be generated from the model
####Y = binomial(mu), mu = a*riskscore + case proportion
####need to estimate the calibration slope a (at least the scale should be correct).

cv.generatey <- function(anc, s.size, beta0, risk.score, case.prop = 0.5, extra.scaling = 1){
  ###case.prop > 0.5 means there is more case than control 
  
  s2 <- sum(risk.score^2) #this is the (rescaled) variance of the risk scores
  
  ##load the original GWAS results
  p.org <- fread(paste0('/raid6/Ron/prs/data/bert_sample/', anc,'.TRN.PHENO1.glm.logistic.hybrid'))
  
  # names(p.org)[1] <- 'chr'
  # p.org <- p.org[chr %in% 1:22,]
  cor.org <- p2cor(p = p.org$P, n = s.size, sign = log(p.org$OR)) #original correlation vector
  
  design.factor <- sqrt(case.prop * (1 - case.prop)) 
  a <- cor.org %*% beta0 * sqrt(s.size) * design.factor / s2 ###model-based formula, I should inlcude the derivation in the paper
  a <- as.numeric(a) * extra.scaling ###sometimes may need a extra.scaling != 1 to generate better cv data sets
  
  mu <-  a * risk.score + case.prop
  
  #make sure mu is between 0 and 1 since it is a probability
  print(paste0('number of beyond 0/1 range', sum(mu > 1 | mu<0)))
  mu[mu > 1] <- 1
  mu[mu < 0] <- 0
  
  # epsilon <- rbinom(s.size,1,mu)
  # epsilon <- as.numeric(epsilon == 1)*(1-mu) + as.numeric(epsilon == 0)*(-mu)
  # booty <- mu + epsilon
  ###
  booty <- rbinom(s.size, 1, mu)
  
  save(booty, file = paste0('/raid6/Tianyu/PRS/BootData/',anc,'_bootY_',setting.title,'.RData'))
  ##"/raid6/Tianyu/PRS/BootData/CEU_bootY_calib.RData"
  
}

split.train.val <- function(chr, anc, train.index, val.index, plink){
  
  if(anc == 'CEU'){
    file.title <- 'CEU-20K'
  }else{
    file.title <- 'YRI-4K'  
  }
  
  ##this is the reference panel genotype data befpre splitting
  referece.panel.name <- paste0("/raid6/Tianyu/PRS/bert_sample/ReferencePopulation-Package/", file.title,"/CHR/",file.title,"-chr", chr)
  boost.data.folder <- paste0("/raid6/Tianyu/PRS/BootData/", file.title,"/CHR")
  boost.data.train.index <- paste0(boost.data.folder, "/boost-train-index.txt")
  boost.data.val.index <- paste0(boost.data.folder, "/boost-val-index.txt")
  boost.data.train.data <- paste0(boost.data.folder, "/", file.title,"-chr", chr, "boost-train")
  boost.data.val.data <- paste0(boost.data.folder, "/", file.title,"-chr", chr, "boost-val")
  
  
  #complete panel information
  # full.fam <- fread(paste0("/raid6/Tianyu/PRS/bert_sample/",anc,".TUNE/CHR/",anc,".TUNE-chr",chr,".fam"))
  full.fam <- fread(paste0(referece.panel.name, ".fam"))
  
  train.fam <- full.fam[train.index, c(1,2)] #only keep family id and withtin family id
  #write down which are in the training set
  # fwrite(train.fam, paste0("/raid6/Tianyu/PRS/BootData/",anc,".TUNE/CHR/",anc,".TUNE-boost-train-index.txt"),
  #        col.names = F, sep = " ")
  fwrite(train.fam, boost.data.train.index,
         col.names = F, sep = " ")
  
  #split the reference panel into training and validation sets
  plink.command <- paste(plink, "--bfile", referece.panel.name,
                         "--allow-no-sex",
                         "--keep", boost.data.train.index,
                         "--make-bed", "--out", boost.data.train.data,
                         "--noweb", "--keep-allele-order",
                         sep = " ")
  system(plink.command)
  
  val.fam <- full.fam[val.index, c(1,2)]
  #write down which are in the validation set
  fwrite(val.fam, boost.data.val.index,
         col.names = F, sep = " ")
  
  plink.command <- paste(plink, "--bfile", referece.panel.name,
                         "--allow-no-sex",
                         "--keep", boost.data.val.index,
                         "--make-bed", "--out", boost.data.val.data,
                         "--noweb", "--keep-allele-order",
                         sep = " ")
  system(plink.command)
  
  print(paste0('finished sample splitting for chr ', chr))
}

cor_bychr_bootstrap <- function(chr, anc, boot.Y, B){
  print(paste0('chr is', chr))
  ##calculate pairwise correlation
  
  # chr_loc <- as.numeric(gsub(':.*$','',names(beta0)))
  # snp <- names(beta0)[chr_loc == chr]
  
  if(anc == 'CEU'){
    file.title <- 'CEU-20K'
  }else{
    file.title <- 'YRI-4K'  
  }
  
  boost.data.folder <- paste0("/raid6/Tianyu/PRS/BootData/", file.title,"/CHR")
  boost.data.train.data <- paste0(boost.data.folder, "/", file.title,"-chr", chr, "boost-train")

  ####load training genotype data 
  system.time(gnt<-read.plink(bed=paste0(boost.data.train.data,".bed"),
                              bim=paste0(boost.data.train.data,".bim"),
                              fam=paste0(boost.data.train.data,".fam"))
  )
  
  # transform the genotypes to a matrix, and reverse the call count. This snpStats package counts the A2 allele, not the A1
  system.time(gnt<-2-as(gnt$genotypes,Class="numeric"))
  
  # center the columns
  system.time(gnt <-gnt - rep(1, nrow(gnt)) %*% t(colMeans(gnt)))
  # normalize the calls to the unit 1 norm
  system.time(gnt<-normalize.cols(gnt,method="euclidean",p=2))
  
  ###calculate pairwise correlation
  print('calculating pairwise correlation')
  print(paste0('B is ', B))
  boot.cor <- matrix(0, nrow = NCOL(gnt), ncol = B)
  for(i in 1:B){
    print(paste0('the No.',i,'boostrap data of ', anc))
    for(j in 1:NCOL(gnt)){
      ###after the proper normalization, correlation is the same as inner product
      boot.cor[j,i] <- crossprod(boot.Y[,i], gnt[,j] )
    }
  }
  # save(boot.cor, file = '/raid6/Tianyu/PRS/trash/beta0is5-YRI-TRN-chr20-cor.RData')
  
  #####make some "fake" GWAS results
  print('translate correlation into p-values')
  fake_GWAS_list <- list()
  sample_size <- NROW(gnt)
  for(i in 1:B){
    t_stat <- boot.cor[,i] * sqrt(sample_size - 2)/sqrt(1 - (boot.cor[,i])^2)
    fake_GWAS <- data.frame(P = 2*pt(q = abs(t_stat), df = sample_size - 2, lower.tail = F), #two-sided p-value
                            OR = 2*as.numeric(boot.cor[,i]>0) + 0.5*as.numeric(boot.cor[,i]<0))
    fake_GWAS$n <- sample_size
    fake_GWAS_list[[i]] <- fake_GWAS
    # fwrite(fake_GWAS, file = paste0('/raid6/Tianyu/PRS/BootData/',anc,'_bootdata_repeat',i,'_chr',chr))
  }
  
  return(fake_GWAS_list)
}
