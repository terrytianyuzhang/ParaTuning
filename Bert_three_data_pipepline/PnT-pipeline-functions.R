#### this contains a list of functions that are used in the PNT part of the  pipeline

##### Used in step 4b. Here the weights from step 4a are applied to our TeSTing set
wpntTesting <- function(TST,WORK.DIR,pnt.weights){
  #### this function simpy applies the weights found from the TUNEing PRS to the TeSTing PRS
  # TST - testing dataset
  # WORK.DIR - Working directory
  # prs.weights - weights to apply to the different prs

  WPNT.DIR=paste0(WORK.DIR,"/",TST,".TST/WPNT-TESTING")
  dir.create(paste0(WPNT.DIR),recursive = T, showWarnings = F)
  
  # find the names of the TRN datasets
  rownames(pnt.weights)=pnt.weights$TUNE
  # change the column names for easy access
  colnames(pnt.weights)=gsub(".wght","",colnames(pnt.weights))
  # retain the columns with usefull information
  pnt.weights=pnt.weights[,2:(ncol(pnt.weights)-1)]
  
  training.sets=colnames(pnt.weights)

  
  # read the PRS sets
  PNT=list()
  for(trn in training.sets){
    PNT[[trn]]=fread(paste0(WORK.DIR,"/",TST,".TST/PNT-TESTING/",trn,".TRN/",trn,".TRN-",TST,".TST.sscore"),header=T,data.table=F)
  }
  
  # determine the weighted PRS
  # general information
  WPNT=PNT[[1]]
  rownames(WPNT)=WPNT$IID
  # collect the scoresum
  WPNT$SCORE1_SUM=0
  for(TRN in training.sets){
    WPNT[PNT[[TRN]]$IID,"SCORE1_SUM"]=WPNT[PNT[[TRN]]$IID,"SCORE1_SUM"]+pnt.weights[TST,TRN]*scale(PNT[[TRN]]$SCORE1_SUM)
  }
  # determine the auc
  re.auc=auc(WPNT$PHENO1-1~WPNT$SCORE1_SUM)[1]
  # collect pertinent information
  WPNT=WPNT[,c("#FID","IID","PHENO1","SCORE1_SUM")]
  fwrite(WPNT,paste0(WPNT.DIR,"/WPNT-",TST,".TST.sscore"),row.names=F,col.names=T,quote=F,sep="\t")
  
  return(data.frame(TST,CEU.wght=pnt.weights[TST,"CEU"],YRI.wght=pnt.weights[TST,"YRI"],auc=re.auc))
  
}



#### Used in step 4a. We take the results from step 3b, tuning of the hyper parameters to get the P&T PRS with the largest AIC
#### The TUNEing results are used to obtain the optimal weights
wpntTuning <- function(TUNE,WORK.DIR,tuning.parameters){
  ### this function find the optiomal weights bewteen the two TRaiNing sets to obtain the weighted P&T
  # TUNE - TUNEing dataset 
  # WORK.DIR - working directory # WORK.DIR=paste0(work.dir,"/",sim)
  # tuning.parameters - the parameters to use foir each combination of TRaiNing and TUNEing data sets
  
  WPNT.DIR=paste0(WORK.DIR,"/",TUNE,".TUNE/WPNT-TUNING")
  dir.create(paste0(WPNT.DIR),recursive = T, showWarnings = F)
  
  # obtain the the two sets of PRS for the selected hyper parameters
  PRS=list()
  for(TRN in tuning.parameters[tuning.parameters$TUNE == TUNE,"TRN"]){
    TUNING.DIR=paste0(WORK.DIR,"/",TUNE,".TUNE/PNT-TUNING/",TRN,".TRN/CHR")
    r2=tuning.parameters[tuning.parameters$TRN == TRN & tuning.parameters$TUNE == TUNE,"r2"]
    pvalue=tuning.parameters[tuning.parameters$TRN == TRN & tuning.parameters$TUNE == TUNE,"pvalue"]
    sscore=NULL
    for(chr in 1:22){
      if(file.exists(paste0(TUNING.DIR,"/",TRN,"-chr",chr,sprintf("-r2_%2.1f",r2),"-p.",pvalue,".sscore"))){
        # guards against a non-existing file
        if(is.null(sscore)){
          sscore=fread(paste0(TUNING.DIR,"/",TRN,"-chr",chr,sprintf("-r2_%2.1f",r2),"-p.",pvalue,".sscore"),header=T,data.table=F)
          sscore$SCORE1_SUM=sscore$SCORE1_AVG*sscore$ALLELE_CT  # this is used to geta  score across the chromosomes
          rownames(sscore)=sscore$IID
        }else{
          tmp=fread(paste0(TUNING.DIR,"/",TRN,"-chr",chr,sprintf("-r2_%2.1f",r2),"-p.",pvalue,".sscore"),header=T,data.table=F)
          if(tmp$ALLELE_CT[1] > 0){
            # guards against non-contributing chromosomes
            tmp$SCORE1_SUM=tmp$SCORE1_AVG*tmp$ALLELE_CT
            # add this contribution to the score
            sscore[tmp$IID,"ALLELE_CT"]=sscore[tmp$IID,"ALLELE_CT"]+tmp$ALLELE_CT
            sscore[tmp$IID,"SCORE1_SUM"]=sscore[tmp$IID,"SCORE1_SUM"]+tmp$SCORE1_SUM
          }
        }
      }
    }
    
    PRS[[TRN]]=sscore
  }
  
  # find the optimal weight
  re.WGHT=NULL
  for(wght in seq(0,1,.01)){
    (auc=auc(PRS$CEU$PHENO1-1~eval(wght*scale(PRS$CEU$SCORE1_SUM)+(1-wght)*scale(PRS$YRI$SCORE1_SUM)))[1])
    re.WGHT=rbind.data.frame(re.WGHT,data.frame(CEU.wght=wght,YRI.wght=1-wght,auc))
  }
  
  fwrite(re.WGHT,paste0(WPNT.DIR,"/WPNT-TUNING.txt"),row.names=F,col.names=T,quote=F,sep="\t")
  return(data.frame(TUNE,re.WGHT[which.max(re.WGHT$auc),]))
}

##### Used in step 3B of the pipeline, this part calculates the PRS given the selected parameters
##### it used the GWAS from the TRainNing set, the parameters from the TUNEing set and the genotypes 
##### from the TeSTing set
pntTesting <- function(TRN,TST,r2,pvalue,WORK.DIR=WORK.DIR,plink2 = plink2){
  ### function to calculate the PRS given a set of r2 and pvalue cut-off
  # TRN - training set
  # TST - testing set
  # r2 - r2 to use (0.2, 0.5, 0.8)
  # p-value - p-value to use
  # WORK.DIR - working directory  # WORK.DIR=paste0(work.dir,"/",sim)
  
  # create a directory for the results of the PRS calculation
  PRS.DIR=paste0(WORK.DIR,"/",TST,".TST/PNT-TESTING/",TRN,".TRN")
  dir.create(paste0(PRS.DIR,"/CHR"),recursive = T, showWarnings = F)
  
  
  prsByCHR<-function(chr,TRN,TST,PRS.DIR,r2,pvalue,plink2){
    ##### calculates the score by chromosome
    # chr - chromosome to work on 
    # TRN - training set from which to use the GWAS
    # TST - testing set from which to use the genotypes
    # PRS.DIR - location in which to store the results
    # r2 - r2 value to use
    # pvalue - pvalue to use
    # plink2 - plink version to use
    
    # read the clumps by chromosome use the tuning set that matches the testing set
    index.snp=fread(paste0(WORK.DIR,"/",TST,".TUNE/PNT-PRUNING/",TRN,".TRN/CHR/",TRN,"-chr",chr,sprintf("-r2_%2.1f",r2),".clumped"),header=T,data.table=F)$SNP
    
    # read the gwas for the training set and select the fields of interst
    trn.gwas=fread(paste0(WORK.DIR,"/",TRN,".TRN/GWAS/CHR/",TRN,"-chr",chr,".PHENO1.glm.logistic.hybrid"),header=T,data.table=F)[,c("ID","A1","OR","P")]
    # select the snp to use
    trn.gwas=trn.gwas[trn.gwas$P <= pvalue,]
    trn.gwas$B=log(trn.gwas$OR)
    # select the index.snp
    trn.gwas=trn.gwas[trn.gwas$ID %in% index.snp,]
    fwrite(trn.gwas,paste0(PRS.DIR,"/CHR/chr",chr,"-selected",".snp"),row.names=F,col.names=T,quote=F,sep="\t")
      
    # use plink to calcualte the PRS score
    plink.command=paste(plink2,"--allow-no-sex","--nonfounders","--threads",8,"--memory",25000,
                          "--bfile",paste0(WORK.DIR,"/",TST,".TST/GNT/CHR/",TST,".TST-chr",chr),
                          "--score",paste0(PRS.DIR,"/CHR/chr",chr,"-selected",".snp"),1,2,5,"header",
                          "--out",paste0(PRS.DIR,"/CHR/",TRN,"-chr",chr),
                          sep=" ")
    system(plink.command)
      
    return(chr)
  }

  ### calculate the scores for each of the chromsomes
  system.time(re.prsBychr<- mclapply(1:22,prsByCHR,TRN=TRN,TST=TST,PRS.DIR=PRS.DIR,r2=r2,pvalue=pvalue,plink2=plink2,
                                        mc.cores=22,mc.preschedule = F) )
  
  
  # calculate the AUC
  # collect the scores for each chromosome
  sscore=NULL
  for(chr in 1:22){
    if(file.exists(paste0(PRS.DIR,"/CHR/",TRN,"-chr",chr,".sscore"))){
      # guards a gaainst a non existing file
      if(is.null(sscore)){
        sscore=fread(paste0(PRS.DIR,"/CHR/",TRN,"-chr",chr,".sscore"),header=T,data.table=F)
        sscore$SCORE1_SUM=sscore$SCORE1_AVG*sscore$ALLELE_CT  # this is used to geta  score across the chromosomes
        rownames(sscore)=sscore$IID
      }else{
       tmp=fread(paste0(PRS.DIR,"/CHR/",TRN,"-chr",chr,".sscore"),header=T,data.table=F)
       if(tmp$ALLELE_CT[1] > 0){
         # guards against a non contributing chromosome
         tmp$SCORE1_SUM=tmp$SCORE1_AVG*tmp$ALLELE_CT
         # add this contribution to the score
         sscore[tmp$IID,"ALLELE_CT"]=sscore[tmp$IID,"ALLELE_CT"]+tmp$ALLELE_CT
         sscore[tmp$IID,"SCORE1_SUM"]=sscore[tmp$IID,"SCORE1_SUM"]+tmp$SCORE1_SUM
       }
      }
    }
  }
  ### calculate the average
  sscore$SCORE1_AVG=sscore$SCORE1_SUM/sscore$ALLELE_CT
  
  # save the results
  fwrite(sscore,paste0(PRS.DIR,"/",TRN,".TRN-",TST,".TST",".sscore"),row.names=F,col.names=T,quote=F,sep="\t")
  
  # calculate the AUC
  re.auc=auc(sscore$PHENO1-1~sscore$SCORE1_SUM)[1]
  
  return(data.frame(TRN,TST,r2,pvalue,auc=re.auc))

  
}






####  used in step 3b of the pipeline , this part selects the best combination of r2 and pvalue for the P&T
#### for this it used the GWAS from TRN, the pruning results from step 3a and the genoptype from TUNE
#### it does this by processing the genotypes by chromosome
pntTuning <- function(TRN,TUNE,WORK.DIR,r2,pval.range.file,plink2){
  # function to select the hyper parameters for P&T, here the hyper parameters is a 
  # combination of r2 and pvalue cutoff
  # this function requires a a file with the range of pvalues to consider
  #
  # TRN - training set
  # TUNE - tuning set
  # WORK.DIR - general work directory                               # WORK.DIR=paste0(work.dir,"/",sim)
  # pval.range.file - file with the p-value cut-offs to consider    # pval.range.file=paste0(data.dir,"/pval-range-list.txt")
  #
  # Output
  #  hyper.r2 - hyper parameter for r2
  #  hyper.pcut - hyper parameter for the pvalue cut-off
  
  # create a directory for the results of the hyper parameter calculation
  TUNING.DIR=paste0(WORK.DIR,"/",TUNE,".TUNE/PNT-TUNING/",TRN,".TRN/CHR")
  dir.create(TUNING.DIR,recursive = T, showWarnings = F)
  
  
  RE.PNT.AUC=NULL
  # for each of the r2 in r2.list
  sscoreByCHR<-function(chr,TRN,TUNE,TUNING.DIR,r2,pval.range.file,plink2){
    ##### calculates the score by chromosome
    # chr - chromosome to work on 
    # TRN - training set from which to use the GWAS
    # TUNE - tuning set from which to use the genotypes
    # TUNING.DIR - location in which to store the results
    # r2 - r2 value to use
    # pval.range.file - loaction and filename with the p-values to use
    # plink2 - plink version to use
    
    # read the clumps by chromosome
    index.snp=fread(paste0(WORK.DIR,"/",TUNE,".TUNE/PNT-PRUNING/",TRN,".TRN/CHR/",TRN,"-chr",chr,sprintf("-r2_%2.1f",r2),".clumped"),header=T,data.table=F)$SNP
    
    # read the gwas for the training set and select the fields of interst
    trn.gwas=fread(paste0(WORK.DIR,"/",TRN,".TRN/GWAS/CHR/",TRN,"-chr",chr,".PHENO1.glm.logistic.hybrid"),header=T,data.table=F)[,c("ID","A1","OR","P")]
    trn.gwas$B=log(trn.gwas$OR)
    # select the index.snp
    trn.gwas=trn.gwas[trn.gwas$ID %in% index.snp,]
    fwrite(trn.gwas,paste0(TUNING.DIR,"/chr",chr,"-selected",sprintf("-r2_%2.1f",r2),".snp"),row.names=F,col.names=T,quote=F,sep="\t")
    
    # use plink to calcualte the PRS score
    plink.command=paste(plink2,"--allow-no-sex","--nonfounders","--threads",8,"--memory",25000,
                        "--bfile",paste0(WORK.DIR,"/",TUNE,".TUNE/GNT/CHR/",TUNE,".TUNE-chr",chr),
                        "--score",paste0(TUNING.DIR,"/chr",chr,"-selected",sprintf("-r2_%2.1f",r2),".snp"),1,2,5,"header",
                        "--q-score-range",pval.range.file,paste0(TUNING.DIR,"/chr",chr,"-selected",sprintf("-r2_%2.1f",r2),".snp"),1,4,"header",
                        "--out",paste0(TUNING.DIR,"/",TRN,"-chr",chr,sprintf("-r2_%2.1f",r2),"-p"),
                        sep=" ")
    system(plink.command)
    
    return(chr)
    
    
  }
  
  ### calculate the scores for each of the chromsomes
  system.time(re.sscoreBychr<- mclapply(1:22,sscoreByCHR,TRN=TRN,TUNE=TUNE,TUNING.DIR=TUNING.DIR,r2=r2,pval.range.file=pval.range.file,plink2=plink2,
                           mc.cores=22,mc.preschedule = F) )
  
  
  # process the scores
  files=list.files(pattern = sprintf("-r2_%2.1f",r2),path = TUNING.DIR,full.names = T)
  files=files[-grep(paste(c(".log",".snp"),collapse="|"),files)]

 
  calculateAUC <- function(pvalue,TRN,TUNE,TUNING.DIR,r2){
    # calculate the AUC fort each of the pvalues
    sscore=NULL
    for(chr in 1:22){
      if(file.exists(paste0(TUNING.DIR,"/",TRN,"-chr",chr,sprintf("-r2_%2.1f",r2),"-p.",pvalue,".sscore"))){
        # guards against a non existing file
        if(is.null(sscore)){
          sscore=fread(paste0(TUNING.DIR,"/",TRN,"-chr",chr,sprintf("-r2_%2.1f",r2),"-p.",pvalue,".sscore"),header=T,data.table=F)
          sscore$SCORE1_SUM=sscore$SCORE1_AVG*sscore$ALLELE_CT  # this is used to geta  score across the chromosomes
          rownames(sscore)=sscore$IID
        }else{
          tmp=fread(paste0(TUNING.DIR,"/",TRN,"-chr",chr,sprintf("-r2_%2.1f",r2),"-p.",pvalue,".sscore"),header=T,data.table=F)
          if(tmp$ALLELE_CT[1] > 0){
            # make sure this chromosome contributes
            tmp$SCORE1_SUM=tmp$SCORE1_AVG*tmp$ALLELE_CT
            # add this contribution to the score
            sscore[tmp$IID,"ALLELE_CT"]=sscore[tmp$IID,"ALLELE_CT"]+tmp$ALLELE_CT
            sscore[tmp$IID,"SCORE1_SUM"]=sscore[tmp$IID,"SCORE1_SUM"]+tmp$SCORE1_SUM
          }
        }
      }
    }
    
    # calculate the AUC
    re.auc=auc(sscore$PHENO1-1~sscore$SCORE1_SUM)[1]
    return(data.frame(TRN,TUNE,r2,pvalue,auc=re.auc))
  }
  
  # find the pvalues to process
  pvalues=unique(gsub(".sscore","",unlist(lapply(strsplit(files,"-p."),`[[`,2))))
  
  # calculate the auc for each of the p-values
  re.AUC<-rbindlist(lapply(pvalues,calculateAUC,TRN = TRN, TUNE = TUNE, TUNING.DIR = TUNING.DIR, r2 = r2))
  
  return(re.AUC)
}


####  used in step 3b of the pipeline , this part selects the best combination of r2 and pvalue for the P&T
#### for this it used the GWAS from TRN, the pruning results from step 3a and the genoptype from TUNE
pntSelectHyperParameters.v1 <- function(TRN,TUNE,WORK.DIR,r2.list,pval.range.file,plink2){
  # function to select the hyper parameters for P&T, here the hyper parameters is a 
  # combination of r2 and pvalue cutoff
  # this function requires a a file with the range of pvalues to consider
  #
  # TRN - training set
  # TUNE - tuning set
  # WORK.DIR - general work directory                               # WORK.DIR=paste0(work.dir,"/",sim)
  # pval.range.file - file with the p-value cut-offs to consider    # pval.range.file=paste0(data.dir,"/pval-range-list.txt")
  # r2.list - list of r2 values for which clumps were determined    # r2.list=c(0.2,0.5,0.8)
  #
  # Output
  #  hyper.r2 - hyper parameter for r2
  #  hyper.pcut - hyper parameter for the pvalue cut-off
  
  # create a directory for the results of the hyper parameter calculation
  TUNING.DIR=paste0(WORK.DIR,"/",TUNE,".TUNE/PNT-TUNING/",TRN,".TRN")
  dir.create(TUNING.DIR,recursive = T, showWarnings = F)
  
  
  RE.PNT.AUC=NULL
  # for each of the r2 in r2.list
  for(r2 in r2.list){
    
    # read the index SNP for each clump
    readClumpsByChr <- function(chr,TRN,TUNE,WORK.DIR,r2){
      # read the clumps by chromosome
      index.snp=fread(paste0(WORK.DIR,"/",TUNE,".TUNE/PNT-PRUNING/",TRN,".TRN/CHR/",TRN,"-chr",chr,sprintf("-r2_%2.1f",r2),".clumped"),header=T,data.table=F)$SNP
      return(index.snp)
    }
    index.snp=unlist(lapply(1:22,readClumpsByChr,TRN=TRN,TUNE=TUNE,WORK.DIR=WORK.DIR,r2=r2))
    
    # read the gwas for the training set and select the fields of interst
    trn.gwas=fread(paste0(WORK.DIR,"/",TRN,".TRN/GWAS/",TRN,".PHENO1.glm.logistic.hybrid"),header=T,data.table=F)[,c("ID","A1","OR","P")]
    trn.gwas$B=log(trn.gwas$OR)
    # select the index.snp
    trn.gwas=trn.gwas[trn.gwas$ID %in% index.snp,]
    fwrite(trn.gwas,paste0(TUNING.DIR,"/selected",sprintf("-r2_%2.1f",r2),".snp"),row.names=F,col.names=T,quote=F,sep="\t")
    
    # use plink to calcualte the PRS score
    plink.command=paste(plink2,"--allow-no-sex","--nonfounders",
                        "--bfile",paste0(WORK.DIR,"/",TUNE,".TUNE/GNT/",TUNE,".TUNE"),
                        "--score",paste0(TUNING.DIR,"/selected",sprintf("-r2_%2.1f",r2),".snp"),1,2,5,"header",
                        "--q-score-range",pval.range.file,paste0(TUNING.DIR,"/selected",sprintf("-r2_%2.1f",r2),".snp"),1,4,"header",
                        "--out",paste0(TUNING.DIR,"/",TRN,sprintf("-r2_%2.1f",r2),"-p"),
                        sep=" ")
    system(plink.command)
    
    # determine the AUC for each of the p-value cutoffs
    sscore.files=list.files(path=TUNING.DIR,pattern = paste0(TRN,sprintf("-r2_%2.1f",r2)),full.names = T)
    sscore.files=sscore.files[-grep("log",sscore.files)]
    for(file in sscore.files){
      prs=fread(file,header=T,data.table=F)
      AUC=auc(prs$PHENO1-1~prs$SCORE1_AVG)[1]
      RE.PNT.AUC=rbind.data.frame(RE.PNT.AUC,data.frame(r2=r2,p=as.numeric(gsub(".sscore","",unlist(strsplit(file,"-p."))[2])),auc=AUC))
    }
    
  }
  fwrite(RE.PNT.AUC,paste0(TUNING.DIR,"/",TRN,"-PNT-hyper.txt"),row.names=F,col.names=T,quote=F,sep="\t")
         
  return(data.frame(TRN,TUNE,RE.PNT.AUC[which.max(RE.PNT.AUC$auc),]))      
  
}







#### used in step 3a of the pipeline, pruning of the SNP based on the r2 of the SNP ain the tuning data and the p-values from the GWAS
pntPruning <- function(TRN,TUNE,WORK.DIR,r2,plink){
  # function to perform the pruning that is needed in the for the P&T genetic risk score. Prune is based on the p-value from the GWAS of the TRaiNing set
  # and the genotypes from the TUNEing set.
  # the appraoch is to perform this by chromosome  
  
  # TRN - training set (CEU or YRI)
  # TUNE - tuning set (CEU or YRI)
  # WORK.DIR - general work directory WORK.DIR=paste0(work.dir,"/",sim)
  # r2 - r2 to use
  
  # Results are written to the TUNING directory
  
  # create a direcory to write the results
  RES.DIR=paste0(WORK.DIR,"/",TUNE,".TUNE/PNT-PRUNING/",TRN,".TRN/CHR")
  dir.create(RES.DIR,recursive = T,showWarnings = F)
  
  clumpingByChr<-function(chr,r2,TRN,TUNE,RES.DIR,plink){
    # function using plink to perform the clumping by chromosome
    # run clumping
    plink.command=paste(plink,"--nonfounders","--allow-no-sex",
                        "--threads",1,
                        "--memory",25000,
                        "--bfile",paste0(WORK.DIR,"/",TUNE,".TUNE/GNT/CHR/",TUNE,".TUNE-chr",chr),
                        "--chr",chr,
                        "--clump-p1",1,
                        "--clump-p2",1,
                        "--clump-r2",r2,
                        "--clump-kb",500,
                        "--clump",paste0(WORK.DIR,"/",TRN,".TRN/GWAS/CHR/",TRN,"-chr",chr,".PHENO1.glm.logistic.hybrid",sep=""),
                        "--clump-snp-field","ID",
                        "--clump-field","P",
                        "--out",paste0(WORK.DIR,"/",TUNE,".TUNE/PNT-PRUNING/",TRN,".TRN/CHR/",TRN,"-chr",chr,sprintf("-r2_%2.1f",r2)),
                        sep=" ")
    system(plink.command)
    return(chr)
  }
  
  re=mclapply(1:22,clumpingByChr,r2=r2,TRN=TRN,TUNE=TUNE,RES.DIR=RES.DIR,plink=plink,mc.cores=22,mc.preschedule = F)
  
  return(re)
}



#### used in step 2 of the pipeline, calculating GWAS for the training set

runGWAS <-function(ANC,WORK.DIR,plink2){
  # function performs the GWAS on the TRN data 
  # ANC - ancestry to use
  # WORK.DIR - working directory #            WORK.DIR=paste0(work.dir,"/",sim,"/",anc,".","TRN")
  
  # create a directory to store the GWAS results
  dir.create(paste0(WORK.DIR,"/GWAS/CHR"),recursive = T,showWarnings = F)
  
  # complete data
  plink2.command=paste(plink2,"--nonfounders","--allow-no-sex","--threads",32,"--memory",25000,
                       "--bfile",paste(WORK.DIR,"/GNT/",ANC,".TRN",sep=""),
                       "--glm","allow-no-covars","omit-ref",
                       "--out",paste(WORK.DIR,"/GWAS/",ANC,sep=""),
                       sep=" ")
  system(plink2.command)
  
  # run the GWAS by chromosome
  gwasByCHR <- function(chr,ANC,WORK.DIR,plink2){
    plink2.command=paste(plink2,"--nonfounders","--allow-no-sex","--threads",16,"--memory",25000,
                         "--bfile",paste(WORK.DIR,"/GNT/",ANC,".TRN",sep=""),
                         "--glm","allow-no-covars","omit-ref",
                         "--chr",chr,
                         "--out",paste(WORK.DIR,"/GWAS/CHR/",ANC,"-chr",chr,sep=""),
                         sep=" ")
    system(plink2.command)
    return(chr) 
  }
  
  # run for each of the chromosomes
  re<-mclapply(1:22,gwasByCHR,ANC=ANC,WORK.DIR=WORK.DIR,plink2=plink2,mc.cores=8,mc.preschedule = F)
  
  
  return(ANC)
}



############## function for step 1


prepareGenotypes<-function(GNT.SET,ANC,ROLE,WORK.DIR,plink2){
  # for a given genotype set this function will make a copy of the data that can be processed by
  # both plink and plink. It will also split the data by chromosome
  
  # 
  # GNT.SET - genotype set to use     # GNT.SET=gsub("XXX",sim,TRN.SETS[["CEU"]])
  # ANC     - Ancestry CEU or YRI
  # ROLE    - TRN, TUNE, or TST
  # WORK.DIR - location where the results should be written to # WORK.DIR=paste0(work.dir,"/",sim,"/",ANC,".",ROLE)
  
  # create the WORKDIR if it does not exist, include one for the chromosome specific data
  dir.create(paste0(WORK.DIR,"/GNT/CHR"),recursive = T,showWarnings = F)
  
  # create a bed copy of the full data
  plink2.command=paste(plink2,"--threads",32,"--memory",25000,
                       "--pfile",GNT.SET,
                       "--make-bed",
                       "--out",paste0(WORK.DIR,"/GNT/",ANC,".",ROLE),
                       sep=" ")
  system(plink2.command)
  
  # and by chromosome
  byCHR <- function(chr,GNT.SET,ANC,ROLE,WORK.DIR,plink2){
    # split by chromosome
    plink2.command=paste(plink2,"--threads",16,"--memory",25000,
                         "--pfile",GNT.SET,
                         "--make-bed",
                         "--chr",chr,
                         "--out",paste0(WORK.DIR,"/GNT/CHR/","",ANC,".",ROLE,"-chr",chr),
                         sep=" ")
    system(plink2.command)
    return(chr=chr)
  }
  
  # run for each of the chromosomes
  re<-mclapply(1:22,byCHR,GNT.SET=GNT.SET,ANC=ANC,ROLE=ROLE,WORK.DIR=WORK.DIR,plink2=plink2,mc.cores=8,mc.preschedule = F)
  
  return(data.frame(ANC,ROLE))
}

