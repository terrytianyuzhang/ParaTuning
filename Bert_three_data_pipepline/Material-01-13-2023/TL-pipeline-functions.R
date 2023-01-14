#### this contains a list of functions that are used for the Transfer Learning part of the pipeline

### step 1 - Determining the lasso sum betas from the TRaiNing data and the LD pattern and GWAS from the GWAS data
tlBetas <- function(TRN, TARGET, WORK.DIR,REF.DIR,ld.populations,plink){
  # this function performs the calculation of the the betas used for the transfer learning
  # 
  
  # TRN - TRaiNing dataset that provides the lassosum betas using the optimal combination of the hyper parameters
  # TARGET - TARGET population for which PRS are to be calculated.
  # WORK.DIR - directory to which results are written WORK.DIR=paste0(work.dir,"/",sim)
  # REF.DIR - directory with the reference population genotypes
  # ld.populations - are used by lassosum to find ld blocks
  # plink - location of plink 1.9 used for calculating the correlations

  # create a directory to store the results
  TL.DIR=paste0(WORK.DIR,"/",TARGET,".TUNE/TL-BETAS")
  dir.create(paste0(TL.DIR,"/WORK"),recursive = T,showWarnings = F)
  
  ### TRAINING DATA
  # Genotypes to use for calculating LD in the TARGET population
  # train_file="ExampleData_1000G_African_train"
  target.ref.file <- paste0(REF.DIR,"/",TARGET,".REF/",TARGET,".REF")
  
  
  # read the GWAS to obtain information on SNP
  map <- fread(paste0(WORK.DIR,"/",TRN,".TRN/GWAS/",TRN,".PHENO1.glm.logistic.hybrid"), header=T, data.table=F)[,c("ID","A1")]
  colnames(map)=c("SNP","A1")
  rownames(map)=map$SNP
  
  # read the lassosum betas for the optimal 
  # read in the lassosum betas
  ls.betas=fread(paste0(WORK.DIR,"/",TRN,".TRN/LS-SELECT-SNP/",TRN,".REF/",TRN,"-selected-snp.wghts"),header=T,data.table=F)
  # find the betas for the optimal set of parameters
  ls.results=fread(paste0(WORK.DIR,"/Results/ls-tuning-parameters.txt"),header=T,data.table=F)
  combn=unlist(ls.results[ls.results$TRN == TRN & ls.results$REF == TRN & ls.results$TUNE == TARGET,"combn"])
  # merge with the map
  map$Beta=0
  map[ls.betas$ID,"Beta"]=ls.betas[,combn]
  
  # write the information to a file
  trn.lasso.file=paste0(TL.DIR,"/",TRN,".TRN-LS.txt")
  fwrite(map,trn.lasso.file, row.names=F,col.names=T,quote=F,sep="\t")
  
  rm(ls.results,map); gc()
  
  # target population gwas results
  target.gwas <- fread(paste0(WORK.DIR,"/",TARGET,".TRN/GWAS/",TARGET,".PHENO1.glm.logistic.hybrid"), data.table = F, header = T)
  # keep required fields
  target.gwas <- data.frame(SNP = target.gwas$ID,
                                  A1 = target.gwas$A1,
                                  beta = round(log(target.gwas$OR),7),
                                  N = target.gwas$OBS_CT,
                                  p = target.gwas$P)
  # save the results
  target.gwas.file=paste0(TL.DIR,"/",TARGET,"-GWAS.txt")
  fwrite(target.gwas,target.gwas.file, row.names=F,col.names=T,quote=F,sep="\t")
  
  # run the transfer learnging beta step
 
    
  cl <- makeCluster(96, type="FORK")
  
  system.time(out.beta <- TL_BETAS(train_file = target.ref.file,
                                 sum_stats_file = trn.lasso.file,
                                 target_sumstats_file = target.gwas.file,
                                 LDblocks = ld.populations[TARGET],
                                 outfile = paste0(TL.DIR,"/WORK/TL"),
                                 cluster= cl,
                                 plink = plink))
  stopCluster(cl)
  
  return(sim)
}


tlTuning <- function(TRN, TARGET, WORK.DIR, plink2){
  # this function takes the beta sets from the tlBetas function and applies them to a tuning dataset of the TARGET ancestry. It then 
  # selects which of the betas provides the best PRS based on AUC.
  # TRN - TRaiNing ancestry
  # TARGET - TARGET ancestry
  # WORK.DIR - location of the working directory to which results are written # WORK.DIR=paste0(work.dir,"/",sim)
  # plink2 - location of the plink 2.0 program
  
  # create the directory to which the results will be written
  TL.TUNING.DIR=paste0(WORK.DIR,"/",TARGET,".TUNE/TL-TUNING")
  dir.create(paste0(TL.TUNING.DIR),recursive = T,showWarnings = F)
  
  
  # read the betas
  beta.all <- fread(paste0(WORK.DIR,"/",TARGET,".TUNE/TL-BETAS/WORK/TL_beta.candidates.fresh.txt"), header=T,data.table=F) ##read in all the beta sets
  
  ##normalize and make sure the signs are correct
  i.beta=grep("betatemp",colnames(beta.all))
  beta.all[,i.beta] <- beta.all[,i.beta] / beta.all$sd
  not.same.A1 <- which(beta.all$V5 != beta.all$A1)
  beta.all[not.same.A1,i.beta] <- -beta.all[not.same.A1,i.beta]
  colnames(beta.all)[c(1:2,i.beta)] <- c("SNP","A1",paste0("beta",i.beta-i.beta[1]))
  # write this file to be used for the PRS calculations
  fwrite(beta.all[,c(1:2,i.beta)],paste0(WORK.DIR,"/",TARGET,".TUNE/TL-BETAS/TL-beta-sets.txt"),row.names=F,col.names=T,quote=F,sep="\t")
  rm(beta.all); gc()
  
  # calculate the PRS for the samples in the TUNEING set of TARGET ancestry
  tmp=fread(paste0(WORK.DIR,"/",TARGET,".TUNE/TL-BETAS/TL-beta-sets.txt"),header=T,data.table=F,nrows = 1)
  plink2.command = paste(plink2,"--nonfounders","--allow-no-sex","--threads",8,"--memory",25000,
                         "--bfile",paste0(WORK.DIR,"/",TARGET,".TUNE/GNT/",TARGET,".TUNE"),
                         "--score",paste0(WORK.DIR,"/",TARGET,".TUNE/TL-BETAS/TL-beta-sets.txt"),"header-read",1,2,
                         "--score-col-nums",paste0("3","-",ncol(tmp)),
                         "--out",paste0(TL.TUNING.DIR,"/",TRN,"-",TARGET,"-TUNING"),
                         sep=" ")
  system(plink2.command)
  
  ### read the scores
  sscore=fread(paste0(TL.TUNING.DIR,"/",TRN,"-",TARGET,"-TUNING",".sscore"),header=T,data.table=F)
  
  # determine the AUC for each of the scores
  auc.tuning=NULL
  scores=colnames(sscore)[grep("AVG",colnames(sscore))]
  for(score in scores){
    auc.tuning=rbind.data.frame(auc.tuning,data.frame(score=score,auc=auc(sscore$PHENO1-1~sscore[,score])[1]))
  }
  # save these results
  fwrite(auc.tuning,paste0(TL.TUNING.DIR,"/",TRN,".TRN-",TARGET,".TUNE-TUNING.txt"),row.names=F,col.names=T,quote=F,sep="\t")
  
  # determine the maximum AUC
  i.max=which.max(auc.tuning$auc)
  
  return(data.frame(TRN,TARGET,auc=auc.tuning$auc[i.max],i=as.numeric(gsub("beta","",unlist(strsplit(auc.tuning$score[i.max],"_"))[1])),
                    combn=gsub("_AVG","",auc.tuning$score[i.max])))
} 


tlTesting <- function(TRN, TST, combn, WORK.DIR, plink2){
  # this function takes the beta sets from the tlBetas function and applies them to a tuning dataset of the TARGET ancestry. It then 
  # selects which of the betas provides the best PRS based on AUC.
  # TRN - TRaiNing ancestry
  # TST - TeSTing dataset 
  # combn - best combination of beta values based on tuning
  # WORK.DIR - location of the working directory to which results are written # WORK.DIR=paste0(work.dir,"/",sim)
  # plink2 - location of the plink 2.0 program
  
  # create the directory to which the results will be written
  TL.TESTING.DIR=paste0(WORK.DIR,"/",TST,".TST/TL-TESTING")
  dir.create(paste0(TL.TESTING.DIR),recursive = T,showWarnings = F)
  

  # calculate the PRS for the samples in the TESTING 
  tmp=fread(paste0(WORK.DIR,"/",TST,".TUNE/TL-BETAS/TL-beta-sets.txt"),header=T,data.table=F,nrows = 1)
  # find the location of the combination
  i.comb=grep(combn,colnames(tmp))
  plink2.command = paste(plink2,"--nonfounders","--allow-no-sex","--threads",8,"--memory",25000,
                         "--bfile",paste0(WORK.DIR,"/",TST,".TST/GNT/",TST,".TST"),
                         "--score",paste0(WORK.DIR,"/",TST,".TUNE/TL-BETAS/TL-beta-sets.txt"),"header-read",1,2,i.comb,
                         "--out",paste0(TL.TESTING.DIR,"/",TRN,"-",TST,"-TESTING"),
                         sep=" ")
  system(plink2.command)
  
  ### read the scores
  sscore=fread(paste0(TL.TESTING.DIR,"/",TRN,"-",TST,"-TESTING.sscore"),header=T,data.table=F)
  
  # determine the AUC for each of the scores
  auc.tuning=NULL
  scores=colnames(sscore)[grep("AVG",colnames(sscore))]
  for(score in scores){
    auc.tuning=rbind.data.frame(auc.tuning,data.frame(score=score,auc=auc(sscore$PHENO1-1~sscore[,score])[1]))
  }
  # save these results
  fwrite(auc.tuning,paste0(TL.TESTING.DIR,"/",TRN,".TRN-",TST,".TST-TESTING.txt"),row.names=F,col.names=T,quote=F,sep="\t")
  
  # determine the maximum AUC
  i.max=which.max(auc.tuning$auc)
  
  return(data.frame(TRN,TST,auc=auc.tuning$auc[i.max],i=as.numeric(gsub("beta","",unlist(strsplit(auc.tuning$score[i.max],"_"))[1])),
                    combn=gsub("_AVG","",auc.tuning$score[i.max])))
} 
 