#### this contains a list of functions that are used in the LASSO (LS) part of the  pipeline

#### joint lassosum - selection of SNP and betas
jlsSelectSnp <- function(TRN1,TRN2, REF1, REF2, WORK.DIR,REF.DIR,ld.populations,N,shrink=0.8,lambda=exp(seq(log(0.01), log(0.02), length.out=10)), gamma=seq(0,1,.1),chrs=1:22){
  # this function performs the selection of SNP and their weights to be used in the joint lassosum
  # the joint lassom sum takes variables for shrinkage, lambda, and gamma. The mylassosum function can only handle a single value shrinkage and gamma
  # while multiple values can be suplied for lambda. To deal with this I create combinations of chr, shrink, and gamma and process these across multiple
  # processors
  
  # TRN1 - First TRaiNing dataset that provides the GWAS 
  # TRN2 - Second TRaiNing dataset that provides the GWAS 
  # REF1 - First REFerence dataset provides genotypes that are used in lassosum for LD patterns
  # REF2 - Second REFerence dataset provides genotypes that are used in lassosum for LD patterns
  # WORK.DIR - directory to which results are written WORK.DIR=paste0(work.dir,"/",sim)
  # REF.DIR - directory with the reference population genotypes
  # ld.populations - are used by lassosum to find ld blocks
  # N - contains number of cases and controls for each population used in the GWAS
  # shrink - shrinkage parameter (suggested shrink=0.9)
  # lambda - tuning parameter (suggested lambda = sort(0.021-exp(seq(log(0.001), log(0.02), length.out=20))))
  # gamma - tuning parameter (suggested gamma = seq(0.72,0.88,0.04))
  # chrs - chromosomes to analyse (mostly for debugging)
  
  
  # create a directory to store the results
  JLS.DIR=paste0(WORK.DIR,"/",TRN1,"-",TRN2,".TRN/JLS-SELECT-SNP/",REF1,"-",REF2,".REF")
  dir.create(paste0(JLS.DIR,"/CHR"),recursive = T,showWarnings = F)
  
  
  # this function selects the lasso SNP by ldblocks
  jlsByGammaChr <- function(i.combn,combn,TRN1,TRN2,REF1,REF2,ld.populations,REF.DIR,N,shrink,lambda){
    chr=combn[i.combn,"chr"]
    shrink=combn[i.combn,"shrink"]
    gamma=combn[i.combn,"gamma"]
    
    LDblocks=list()
    # obtain the LD blocks
    LDblocks[[TRN1]]=fread(system.file(paste0("data/Berisa.", 
                                              ld.populations[TRN1], ".bed"), 
                                       package="lassosum"), header=T,data.table=F)
    LDblocks[[TRN1]]=LDblocks[[TRN1]][LDblocks[[TRN1]][,"chr"] == paste0("chr",chr),]
    
    LDblocks[[TRN2]]=fread(system.file(paste0("data/Berisa.", 
                                              ld.populations[TRN2], ".bed"), 
                                       package="lassosum"), header=T,data.table=F)
    LDblocks[[TRN2]]=LDblocks[[TRN2]][LDblocks[[TRN2]][,"chr"] == paste0("chr",chr),]
    
    
    
    # read the training GWAS for both populations
    trn.gwas=list()
    trn.gwas[[TRN1]]=fread(paste0(WORK.DIR,"/",TRN1,".TRN/GWAS/CHR/",TRN1,"-chr",chr,".PHENO1.glm.logistic.hybrid"),header=T,data.table=F)
    rownames(trn.gwas[[TRN1]])=trn.gwas[[TRN1]][,"ID"]
    
    trn.gwas[[TRN2]]=fread(paste0(WORK.DIR,"/",TRN2,".TRN/GWAS/CHR/",TRN2,"-chr",chr,".PHENO1.glm.logistic.hybrid"),header=T,data.table=F)
    rownames(trn.gwas[[TRN2]])=trn.gwas[[TRN2]][,"ID"]
    
    
    
    # calculate the correlation for both populations
    r=list()
    r[[TRN1]] <- p2cor(p = trn.gwas[[TRN1]][,"P"], n = N[[TRN1]]$n.case+N[[TRN1]]$n.control, sign=log(trn.gwas[[TRN1]][,"OR"]))
    names(r[[TRN1]])=rownames(trn.gwas[[TRN1]])
    
    r[[TRN2]] <- p2cor(p = trn.gwas[[TRN2]][,"P"], n = N[[TRN2]]$n.case+N[[TRN2]]$n.control, sign=log(trn.gwas[[TRN2]][,"OR"]))
    names(r[[TRN2]])=rownames(trn.gwas[[TRN2]])
    
    # read the fam file fro the REFERENCE data so that a maximum of 5000 samples can be used
    fam=list()
    fam[[REF1]]=fread(paste0(REF.DIR,"/",REF1,".REF/CHR/",REF1,".REF-chr",chr,".fam"),header=F,data.table=F)
    fam[[REF1]]=fam[[REF1]][sample(min(nrow(fam[[REF1]]),5000),replace = F),1:2]
    
    fam[[REF2]]=fread(paste0(REF.DIR,"/",REF2,".REF/CHR/",REF2,".REF-chr",chr,".fam"),header=F,data.table=F)
    fam[[REF2]]=fam[[REF2]][sample(min(nrow(fam[[REF2]]),5000),replace = F),1:2]
    
    
    # run the joint lassosum
    re.jls<-mylassosum(cor1=r[[TRN1]], cor2=r[[TRN2]],
                       fileName1 = paste0(REF.DIR,"/",REF1,".REF/CHR/",REF1,".REF-chr",chr), 
                       fileName2 = paste0(REF.DIR,"/",REF2,".REF/CHR/",REF2,".REF-chr",chr),  
                       gamma = gamma, lambda = lambda, shrink=shrink,
                       chunk=TRUE, mem.limit=10e9,
                       trace=-1,
                       keep1=fam[[REF1]][,"V2"],keep2=fam[[REF2]][,"V2"],
                       LDblocks1=LDblocks[[TRN1]], LDblocks2=LDblocks[[TRN2]])
    
    # give the output row and column names
    colnames(re.jls$beta)=sprintf("s%.2f-g%.2f-l%.5f",re.jls$shrink,re.jls$gamma,re.jls$lambda)
    rownames(re.jls$beta)=rownames(trn.gwas[[TRN1]])
    
    return(re.jls$beta)
    
  }
  
  # determine the combinations of chromosomes and gammas
  # the program can only handle 1 gamma per run
  combn=data.frame( expand.grid(gamma,chrs,shrink)[,c(2,3,1)])
  colnames(combn)=c("chr","shrink","gamma")
  
  # run all chr by gamma combinations
  system.time(RE.jls <- mclapply(1:nrow(combn),jlsByGammaChr,combn=combn,TRN1 = TRN1, TRN2 = TRN2, REF1 = REF1, REF2 = REF2,
                                 ld.populations=ld.populations,REF.DIR=REF.DIR,N=N,shrink = shrink, lambda = lambda,mc.cores = 150, mc.preschedule = F))
  #  save(RE.jls,file="temp.RData")
  # now combine the results
  re.JLS=NULL
  for(chr in unique(combn$chr)){
    #    print(chr);flush.console()
    re.CHR=NULL
    i.chr= which( combn$chr == chr)
    for(i in i.chr){
      if(is.null(re.CHR)){
        re.CHR=RE.jls[[i]]
      }else{
        re.CHR=cbind.data.frame(re.CHR,RE.jls[[i]])
      }
    }
    # write the chromosome specific results
    # round betas to 6 digits
    re.CHR=round(re.CHR,6)
    re.CHR=cbind.data.frame(CHR = unlist(lapply(strsplit(rownames(re.CHR),":"),`[[`,1)),ID = rownames(re.CHR) ,A1=unlist(lapply(strsplit(rownames(re.CHR),":"),`[[`,4)),re.CHR)
    # save the chromosome specific betas
    fwrite(re.CHR,paste0(JLS.DIR,"/CHR/",REF1,"-",REF2,"-chr",chr,"-selected-snp.wghts"),row.names=F,col.names=T,quote=F,sep="\t")
    re.JLS=rbind.data.frame(re.JLS,re.CHR)
  }
  # save all the betas
  fwrite(re.JLS,paste0(JLS.DIR,"/",REF1,"-",REF2,"-selected-snp.wghts"),row.names=F,col.names=T,quote=F,sep="\t")
  
  return("done")
}


# use the TUNEing set to find the optimal combination of hyper parameters

jlsTuning <- function(TRN1, TRN2, REF1, REF2, TUNE, WORK.DIR, plink2){
  # this function takes the betas from the joint lassosum snp selection and determine the PRS based on these using a TUNEing 
  # datasest. It then determines the combination of s and lambda that provides the maximum AUC
  # TRN1 - first TRaiNing dataset for GWAS
  # TRN2 - second TRaiNing dataset for GWAS
  # REF1 - first REFerence dataset used to select the SNP
  # REF2 - second REFerence dataset used to select the SNP
  # TUNE - TUNEing directory
  # WORK.DIR - location of the working directory to which results are written # WORK.DIR=paste0(work.dir,"/",sim)
  # plink2 - location of the plink 2.0 program
  
  # create the directory to which the results will be written
  JLS.DIR=paste0(WORK.DIR,"/",TUNE,".TUNE/JLS-TUNING/",TRN1,"-",TRN2,".TRN/",REF1,"-",REF2,".REF")
  dir.create(paste0(JLS.DIR,"/CHR"),recursive = T,showWarnings = F)
  
  # function
  tuningByCHR <- function(chr,TRN1,TRN2,REF1,REF2,TUNE,WORK.DIR,JLS.DIR,plink2){
    # performs the calculations on a per chromosome basis allowing for the use of multiple processors
    # chr - chromosome to process
    # TRN1 - first TRaiNing dataset
    # TRN2 - second TRaiNing dataset
    # REF1 - first REFerence dataset used to select the SNP
    # REF2 - second REFerence dataset used to select the SNP
    # TUNE - TUNEing directory
    # WORK.DIR - location of the working directory to which results are written # WORK.DIR=paste0(work.dir,"/",sim)
    # JLS.DIR - directory to which the joint lasso results are written4
    # plink2 - location of the plink 2.0 program
    
    
    score.file=paste0(WORK.DIR,"/",TRN1,"-",TRN2,".TRN/JLS-SELECT-SNP/",REF1,"-",REF2,".REF/CHR/",REF1,"-",REF2,"-chr",chr,"-selected-snp.wghts")
    tune.file=paste0(WORK.DIR,"/",TUNE,".TUNE/GNT/CHR/",TUNE,".TUNE-chr",chr)
    tmp=fread(score.file,header=T,data.table=F,nrow=1)
    
    plink2.command = paste(plink2,"--nonfounders","--allow-no-sex","--threads",8,"--memory",25000,
                           "--bfile",tune.file,
                           "--score",score.file,"header-read",2,3,
                           "--score-col-nums",paste0("4","-",ncol(tmp)),
                           "--out",paste0(JLS.DIR,"/CHR/",TRN1,"-",TRN2,".TRN-",REF1,"-",REF2,".REF-chr",chr,"-TUNING"),
                           sep=" ")
    system(plink2.command)
    
    # read the results back in and return them to to the main function
    re.sscore=fread(paste0(JLS.DIR,"/CHR/",TRN1,"-",TRN2,".TRN-",REF1,"-",REF2,".REF-chr",chr,"-TUNING.sscore"),header=T,data.table=F)
    
    # replace the average scores by the sum
    i.avg=grep("AVG",colnames(re.sscore))
    re.sscore[,i.avg]=re.sscore[,i.avg]*re.sscore$ALLELE_CT
    colnames(re.sscore)=gsub("AVG","SUM",colnames(re.sscore))
    
    return(re.sscore)
    
  }
  
  # get teh scores for each chromosome
  system.time(re.SSCORE<-mclapply(1:22,tuningByCHR,TRN1 = TRN1, TRN2 = TRN2,REF1 = REF1,REF2 = REF2,TUNE = TUNE,
                                  WORK.DIR = WORK.DIR,JLS.DIR = JLS.DIR,plink2 = plink2,mc.cores = 22,mc.preschedule = F))
  
  
  # get the total score
  sscore=re.SSCORE[[1]]
  rownames(sscore)=sscore[,"IID"]
  scores=colnames(sscore)[6:ncol(sscore)]
  for(i in 2:length(re.SSCORE)){
    sscore[re.SSCORE[[i]]$IID,colnames(re.SSCORE[[i]])[-c(1:3)]]=sscore[re.SSCORE[[i]]$IID,colnames(re.SSCORE[[i]])[-c(1:3)]]+re.SSCORE[[i]][,colnames(re.SSCORE[[i]])[-c(1:3)]]
  }
  # save the total score
  fwrite(sscore,paste0(JLS.DIR,"/",TRN1,"-",TRN2,".TRN-",REF1,"-",REF2,".REF","-TUNING.sscore"),row.names=F,col.names=T,quote=F,sep="\t")
  
  # calculate the auc for each of the combination of parameters
  auc.tuning=NULL
  scores=colnames(sscore)[grep("SUM",colnames(sscore))][-1]
  for(score in scores){
    auc.tuning=rbind.data.frame(auc.tuning,data.frame(score=gsub("_SUM","",score),auc=auc(sscore$PHENO1-1~sscore[,score])[1]))
  }
  # save these results
  fwrite(auc.tuning,paste0(JLS.DIR,"/",TRN1,"-",TRN2,".TRN-",REF1,"-",REF2,".REF","-TUNING.txt"),row.names=F,col.names=T,quote=F,sep="\t")
  
  
  # determine the maximum AUC
  i.max=which.max(auc.tuning$auc)
  
  return(data.frame(TRN1,TRN2,REF1,REF2,TUNE,auc=auc.tuning$auc[i.max],shrink=as.numeric(gsub("s","",unlist(strsplit(auc.tuning$score[i.max],"-"))[1])),
                    gamma=as.numeric(gsub("g","",unlist(strsplit(auc.tuning$score[i.max],"-"))[2])),
                    lambda=as.numeric(gsub("l","",unlist(strsplit(auc.tuning$score[i.max],"-"))[3])),
                    combn=gsub("_SUM","",auc.tuning$score[i.max])))
  
  
}


##### this part calculates the PRS given the selected lassosum parameters
##### it used the GWAS from the two TRainNing sets, the parameters from the two REFerence sets and TUNEing set combination and the genotypes 
##### from the TeSTing set
jlsTesting <- function(TRN1, TRN2, REF1, REF2, TST,combn,WORK.DIR=WORK.DIR,plink2 = plink2){
  ### function to calculate the LS given a set of hyper-parameters shrink, gamma, and lambda (these are coded in combn)
  # TRN1 - first training set
  # TRN2 - second training set
  # REF1 - first reference set
  # REF2 - second reference set
  # TST - testing set
  # combn - combination of shrinkage, gamma and lambda parameters
  # WORK.DIR - working directory  # WORK.DIR=paste0(work.dir,"/",sim)
  
  # create a directory for the results of the PRS calculation
  JLS.DIR=paste0(WORK.DIR,"/",TST,".TST/JLS-TESTING/",TRN1,"-",TRN2,".TRN/",REF1,"-",REF2,".REF")
  dir.create(paste0(JLS.DIR,"/CHR"),recursive = T, showWarnings = F)
  
  
  jlsByCHR<-function(chr,TRN1,TRN2,REF1,REF2,TST,JLS.DIR,combn,plink2){
    ##### calculates the score by chromosome
    # chr - chromosome to work on 
    # TRN1 - first training set
    # TRN2 - second training set
    # REF1 - first reference set
    # REF2 - second reference set
    # TST - testing set from which to use the genotypes
    # LS.DIR - location in which to store the results
    # comb - combination of shrinkage and tuning parameter lambda s.ssss_l.llll
    # plink2 - plink version to use
    
    score.file=paste0(WORK.DIR,"/",TRN1,"-",TRN2,".TRN/JLS-SELECT-SNP/",REF1,"-",REF2,".REF/CHR/",REF1,"-",REF2,"-chr",chr,"-selected-snp.wghts")
    tst.file=paste0(WORK.DIR,"/",TST,".TST/GNT/CHR/",TST,".TST-chr",chr)
    tmp=fread(score.file,header=T,data.table=F,nrow=1)
    i.col=which(colnames(tmp) == combn)
    
    plink2.command = paste(plink2,"--nonfounders","--allow-no-sex","--threads",8,"--memory",25000,
                           "--bfile",tst.file,
                           "--score",score.file,"header-read",2,3,
                           "--score-col-nums",i.col,
                           "--out",paste0(JLS.DIR,"/CHR/",TRN1,"-",TRN2,".TRN-",REF1,"-",REF2,".REF-chr",chr,"-TESTING"),
                           sep=" ")
    system(plink2.command)
    
    # read the results back in and return them to to the main function
    re.sscore=fread(paste0(JLS.DIR,"/CHR/",TRN1,"-",TRN2,".TRN-",REF1,"-",REF2,".REF-chr",chr,"-TESTING.sscore"),header=T,data.table=F)
    
    # replace the average scores by the sum
    i.avg=grep("AVG",colnames(re.sscore))
    re.sscore[,i.avg]=re.sscore[,i.avg]*re.sscore$ALLELE_CT
    colnames(re.sscore)=gsub("AVG","SUM",colnames(re.sscore))
    
    return(re.sscore)
    
  }
  
  ### calculate the scores for each of the chromsomes
  system.time(re.SSCORE<- mclapply(1:22,jlsByCHR,TRN1=TRN1,TRN2=TRN2,REF1=REF1,REF2=REF2,TST=TST,JLS.DIR=JLS.DIR,combn=combn,plink2=plink2,
                                   mc.cores=22,mc.preschedule = F) )
  
  
  # get the total score
  sscore=re.SSCORE[[1]]
  rownames(sscore)=sscore[,"IID"]
  scores=colnames(sscore)[6:ncol(sscore)]
  for(i in 2:length(re.SSCORE)){
    sscore[re.SSCORE[[i]]$IID,colnames(re.SSCORE[[i]])[-c(1:3)]]=sscore[re.SSCORE[[i]]$IID,colnames(re.SSCORE[[i]])[-c(1:3)]]+re.SSCORE[[i]][,colnames(re.SSCORE[[i]])[-c(1:3)]]
  }
  # save the total score
  # fwrite(sscore,paste0(JLS.DIR,"/",TRN1,"-",TRN2,".TRN-",REF1,"-",REF2,".REF","-TESTING.sscore"),row.names=F,col.names=T,quote=F,sep="\t")
  
  # calculate the auc for each of the combination of parameters
  auc.testing=NULL
  scores=colnames(sscore)[grep("SUM",colnames(sscore))][-1]
  for(score in scores){
    auc.testing=rbind.data.frame(auc.testing,data.frame(score=score,auc=auc(sscore$PHENO1-1~sscore[,score])[1]))
  }
  
  # save the results
  fwrite(sscore,paste0(JLS.DIR,"/",TRN1,"-",TRN2,".TRN-",REF1,"-",REF2,".REF-",TST,".TST",".sscore"),row.names=F,col.names=T,quote=F,sep="\t")
  
  # prepare results
  i.max=which.max(auc.testing$auc)
  return(data.frame(TRN1,TRN2,REF1,REF2,TST,auc=auc.testing$auc[i.max],shrink=as.numeric(gsub("s","",unlist(strsplit(gsub("_SUM","",auc.testing$score[i.max]),"-"))[1])),
                    gamma=as.numeric(gsub("g","",unlist(strsplit(gsub("_SUM","",auc.testing$score[i.max]),"-"))[2])),
                    lambda=as.numeric(gsub("l","",unlist(strsplit(gsub("_SUM","",auc.testing$score[i.max]),"-"))[3])),
                    combn=gsub("_SUM","",auc.testing$score[i.max])))  
  
}




### step 1 - Determining the betas from the TRaNing data and the LD pattern in the REFerence population
lsSelectSNP <- function(TRN, REF, WORK.DIR,REF.DIR,ld.populations,N,s=seq(0.8,1.0,.05),lambda=exp(seq(log(0.001), log(0.02), length.out=20)), chrs=1:22){
  # this function performs the selection of SNP and their weights to be used in lassosum
  # the computing strategy is to do this by ld-block as defined by Berisa and Pickrell (2015)
  
  # TRN - TRaiNing dataset that provides the GWAS 
  # REF - REFerence dataset provides genotypes that are used in lassosum for LD calculations
  # WORK.DIR - directory to which results are written WORK.DIR=paste0(work.dir,"/",sim)
  # REF.DIR - directory with the reference population genotypes
  # ld.pops - are used by lassosum to find ld blocks
  # N - contains number of cases and controls for each population
  # s - shrinkage parameter (suggested s=seq(0.8,1.0,.05))
  # lambda - tuning parameter (suggested lambda=exp(seq(log(0.001), log(0.02), length.out=20))
  # chrs - chromosomes to analyse (mostly for debugging)
  
  # create a directory to store the results
  LASSO.DIR=paste0(WORK.DIR,"/",TRN,".TRN/LS-SELECT-SNP/",REF,".REF")
  dir.create(paste0(LASSO.DIR,"/CHR"),recursive = T,showWarnings = F)
  
  

  #### run by ld.block
  ld.map=fread(paste0(REF.DIR,"/LD/ld.block",".map"),header=T,data.table=F)

  # this function selects the lasso SNP by ldblocks
  byLDBlock <- function(ld.block,TRN,REF,REF.DIR,N,s,lambda){
    # read the ld blocks
    chr=as.numeric(gsub("chr","",unlist(strsplit(ld.block,"-"))[1]))
    # read all the SNP on the chromosome and their ld.block indformation
    ld.snp=fread(paste0(REF.DIR,"/LD/CHR/ld.block-chr",chr,".map"),header=T,data.table=F)
    # select the ones that are in the ld.block
    ld.snp=ld.snp[ld.snp[,paste0(REF,".blk")] == ld.block,"ID"]

    # read the training GWAS
    trn.gwas=fread(paste0(WORK.DIR,"/",TRN,".TRN/GWAS/CHR/",TRN,"-chr",chr,".PHENO1.glm.logistic.hybrid"),header=T,data.table=F)
    rownames(trn.gwas)=trn.gwas$ID
    trn.gwas=trn.gwas[ld.snp,]
    
    # calculate the correlation
    r <- p2cor(p = trn.gwas$P, n = N[[TRN]]$n.case+N[[TRN]]$n.control, sign=log(trn.gwas$OR))
    names(r)=rownames(trn.gwas)
  
    # run the lassosum code
    system.time(re.chr <- lassosum.pipeline(cor=r, chr=trn.gwas$`#CHROM`, snp=trn.gwas$ID, pos=trn.gwas$POS, 
                                A1=trn.gwas$A1, A2=trn.gwas$REF, # A2 is not required but advised
                                ref.bfile=paste0(REF.DIR,"/",REF,".REF/CHR/",REF,".REF-chr",chr), 
                                test.bfile=NULL, #paste0(WORK.DIR,"/",TUNE,".TUNE/GNT/CHR/",TUNE,".TUNE-chr",chr),
                                s=s,lambda=lambda,
                                cluster = NULL,
                                sample=min(N[[REF]]$n.case+N[[REF]]$n.control,5000),
                                LDblocks = ld.populations[REF]))
    
    # set up the output
    shrink.names=names(re.chr$beta)
    re.out=data.frame(CHR=re.chr$sumstats$chr,ID=re.chr$sumstats$snp,A1=re.chr$sumstats$A1)
    for(shrink in shrink.names){
      re.out=cbind.data.frame(re.out,re.chr$beta[[shrink]])
    }
    column.names=expand.grid(round(lambda,5),s)
    column.names=c("CHR","ID","A1",paste(column.names[,2],column.names[,1],sep="_"))
    colnames(re.out)=column.names
                   
    # send results back               
    return(re.out)
  }
  
  ### this allows you to select certain chromosomes
 
  ld.blocks=unique(ld.map[ld.map$CHROM %in% paste0("chr",chr),paste0(REF,".blk")])
  
  re.ls.select<-mclapply(ld.blocks,byLDBlock,TRN=TRN,REF=REF,REF.DIR,N=N,s=s,lambda=lambda,mc.cores=min(length(ld.blocks),150),mc.preschedule = F, mc.silent=T)
  re.ls.select<-rbindlist(re.ls.select)
  
  # round to 6 digits
  re.ls.select[,4:ncol(re.ls.select)]=round(re.ls.select[,4:ncol(re.ls.select)],6)
  # write the result to a file
  fwrite(re.ls.select,paste0(LASSO.DIR,"/",REF,"-selected-snp.wghts"),row.names=F,col.names=T,quote=F,sep="\t")
  
  # and by chromosome
  for(chr in unique(re.ls.select$CHR)){
    fwrite(re.ls.select[re.ls.select$CHR == chr,],paste0(LASSO.DIR,"/CHR/",REF,"-chr",chr,"-selected-snp.wghts"),row.names=F,col.names=T,quote=F,sep="\t")
  }
  
  rm(re.ls.select); gc()
  
  return(c(TRN,REF))
}


lsTuning <- function(TRN, REF, TUNE, WORK.DIR, plink2){
  # this function takes the betas from the lassosum snp selection and determine the PRS based on these using a TUNEing 
  # datasest. It then determines the combination of s and lambda that provides the maximum AUC
  # TRN - TRaiNing dataset
  # REF - REFerence dataset used to select the SNP
  # TUNE - TUNEing directory
  # WORK.DIR - location of the working directory to which results are written # WORK.DIR=paste0(work.dir,"/",sim)
  # plink2 - location of the plink 2.0 program
  
  # create the directory to which the results will be written
  LASSO.DIR=paste0(WORK.DIR,"/",TUNE,".TUNE/LS-TUNING/",TRN,".TRN/",REF,".REF")
  dir.create(paste0(LASSO.DIR,"/CHR"),recursive = T,showWarnings = F)
  
  # function
  tuningByCHR <- function(chr,TRN,REF,TUNE,WORK.DIR,LASSO.DIR,plink2){
    # performs the calculations on a per chromosome basis allowing for the use of multiple processors
    # chr - chromosome to process
    # TRN - TRaiNing dataset
    # REF - REFerence dataset used to select the SNP
    # TUNE - TUNEing directory
    # WORK.DIR - location of the working directory to which results are written # WORK.DIR=paste0(work.dir,"/",sim)
    # LASSO.DIR - directory to which the lasso results are written4
    # plink2 - location of the plink 2.0 program
    
  
    score.file=paste0(WORK.DIR,"/",TRN,".TRN/LS-SELECT-SNP/",REF,".REF/CHR/",REF,"-chr",chr,"-selected-snp.wghts")
    tune.file=paste0(WORK.DIR,"/",TUNE,".TUNE/GNT/CHR/",TUNE,".TUNE-chr",chr)
    tmp=fread(score.file,header=T,data.table=F,nrow=1)
    
    plink2.command = paste(plink2,"--nonfounders","--allow-no-sex","--threads",8,"--memory",25000,
                          "--bfile",tune.file,
                          "--score",score.file,"header-read",2,3,
                          "--score-col-nums",paste0("4","-",ncol(tmp)),
                          "--out",paste0(LASSO.DIR,"/CHR/",TRN,".TRN-",REF,".REF-chr",chr,"-TUNING"),
                          sep=" ")
    system(plink2.command)
    
    # read the results back in and return them to to the main function
    re.sscore=fread(paste0(LASSO.DIR,"/CHR/",TRN,".TRN-",REF,".REF-chr",chr,"-TUNING.sscore"),header=T,data.table=F)
    
    # replace the average scores by the sum
    i.avg=grep("AVG",colnames(re.sscore))
    re.sscore[,i.avg]=re.sscore[,i.avg]*re.sscore$ALLELE_CT
    colnames(re.sscore)=gsub("AVG","SUM",colnames(re.sscore))
    
    return(re.sscore)
  
  }
  
  # get teh scores for each chromosome
  re.SSCORE=mclapply(1:22,tuningByCHR,TRN = TRN,REF = REF,TUNE = TUNE,WORK.DIR = WORK.DIR,LASSO.DIR = LASSO.DIR,plink2 = plink2,mc.cores = 22,mc.preschedule = F)
  
  
  # get the total score
  sscore=re.SSCORE[[1]]
  rownames(sscore)=sscore[,"IID"]
  scores=colnames(sscore)[6:ncol(sscore)]
  for(i in 2:length(re.SSCORE)){
    sscore[re.SSCORE[[i]]$IID,colnames(re.SSCORE[[i]])[-c(1:3)]]=sscore[re.SSCORE[[i]]$IID,colnames(re.SSCORE[[i]])[-c(1:3)]]+re.SSCORE[[i]][,colnames(re.SSCORE[[i]])[-c(1:3)]]
  }
  # save the total score
  fwrite(sscore,paste0(LASSO.DIR,"/",TRN,".TRN-",REF,".REF","-TUNING.sscore"),row.names=F,col.names=T,quote=F,sep="\t")
  
  # calculate the auc for each of the combination of parameters
  auc.tuning=NULL
  scores=colnames(sscore)[grep("SUM",colnames(sscore))][-1]
  for(score in scores){
    auc.tuning=rbind.data.frame(auc.tuning,data.frame(score=score,auc=auc(sscore$PHENO1-1~sscore[,score])[1]))
  }
  # save these results
  fwrite(auc.tuning,paste0(LASSO.DIR,"/",TRN,".TRN-",REF,".REF","-TUNING.txt"),row.names=F,col.names=T,quote=F,sep="\t")
  
  
  # determine the maximum AUC
  i.max=which.max(auc.tuning$auc)
  
  return(data.frame(TRN,REF,TUNE,auc=auc.tuning$auc[i.max],s=unlist(strsplit(auc.tuning$score[i.max],"_"))[1],
                    lambda=unlist(strsplit(auc.tuning$score[i.max],"_"))[2],combn=gsub("_SUM","",auc.tuning$score[i.max])))
  
  
}
 

##### this part calculates the PRS given the selected lassosum parameters
##### it used the GWAS from the TRainNing set, the parameters from the REFerence and TUNEing set combination and the genotypes 
##### from the TeSTing set
lsTesting <- function(TRN,REF,TST,combn,WORK.DIR=WORK.DIR,plink2 = plink2){
  ### function to calculate the LS given a set of hyper-parameters s and lambda (these are coded in combn)
  # TRN - training set
  # REF - reference set
  # TST - testing set
  # combn - combination of shrinkage and tuning parameter lambda s.ssss_l.llll
  # WORK.DIR - working directory  # WORK.DIR=paste0(work.dir,"/",sim)
  
  # create a directory for the results of the PRS calculation
  LS.DIR=paste0(WORK.DIR,"/",TST,".TST/LS-TESTING/",TRN,".TRN/",REF,".REF")
  dir.create(paste0(LS.DIR,"/CHR"),recursive = T, showWarnings = F)
  
  
  lsByCHR<-function(chr,TRN,REF,TST,LS.DIR,combn,plink2){
    ##### calculates the score by chromosome
    # chr - chromosome to work on 
    # TRN - training set from which to use the GWAS
    # REF - referesence set
    # TST - testing set from which to use the genotypes
    # LS.DIR - location in which to store the results
    # comb - combination of shrinkage and tuning parameter lambda s.ssss_l.llll
    # plink2 - plink version to use
    
    score.file=paste0(WORK.DIR,"/",TRN,".TRN/LS-SELECT-SNP/",REF,".REF/CHR/",REF,"-chr",chr,"-selected-snp.wghts")
    tst.file=paste0(WORK.DIR,"/",TST,".TST/GNT/CHR/",TST,".TST-chr",chr)
    tmp=fread(score.file,header=T,data.table=F,nrow=1)
    i.col=which(colnames(tmp) == combn)
    
    plink2.command = paste(plink2,"--nonfounders","--allow-no-sex","--threads",8,"--memory",25000,
                           "--bfile",tst.file,
                           "--score",score.file,"header-read",2,3,
                           "--score-col-nums",i.col,
                           "--out",paste0(LS.DIR,"/CHR/",TRN,".TRN-",REF,".REF-chr",chr,"-TESTING"),
                           sep=" ")
    system(plink2.command)
    
    # read the results back in and return them to to the main function
    re.sscore=fread(paste0(LS.DIR,"/CHR/",TRN,".TRN-",REF,".REF-chr",chr,"-TESTING.sscore"),header=T,data.table=F)
    
    # replace the average scores by the sum
    i.avg=grep("AVG",colnames(re.sscore))
    re.sscore[,i.avg]=re.sscore[,i.avg]*re.sscore$ALLELE_CT
    colnames(re.sscore)=gsub("AVG","SUM",colnames(re.sscore))
    
    return(re.sscore)
    
  }
  
  ### calculate the scores for each of the chromsomes
  system.time(re.SSCORE<- mclapply(1:22,lsByCHR,TRN=TRN,REF=REF,TST=TST,LS.DIR=LS.DIR,combn=combn,plink2=plink2,
                                     mc.cores=22,mc.preschedule = F) )
  
  
  # get the total score
  sscore=re.SSCORE[[1]]
  rownames(sscore)=sscore[,"IID"]
  scores=colnames(sscore)[6:ncol(sscore)]
  for(i in 2:length(re.SSCORE)){
    sscore[re.SSCORE[[i]]$IID,colnames(re.SSCORE[[i]])[-c(1:3)]]=sscore[re.SSCORE[[i]]$IID,colnames(re.SSCORE[[i]])[-c(1:3)]]+re.SSCORE[[i]][,colnames(re.SSCORE[[i]])[-c(1:3)]]
  }
  # save the total score
  fwrite(sscore,paste0(LS.DIR,"/",TRN,".TRN-",REF,".REF","-TESTING.sscore"),row.names=F,col.names=T,quote=F,sep="\t")
  
  # calculate the auc for each of the combination of parameters
  auc.testing=NULL
  scores=colnames(sscore)[grep("SUM",colnames(sscore))][-1]
  for(score in scores){
    auc.testing=rbind.data.frame(auc.testing,data.frame(score=score,auc=auc(sscore$PHENO1-1~sscore[,score])[1]))
  }
  
  # save the results
  fwrite(sscore,paste0(LS.DIR,"/",TRN,".TRN-",REF,".REF-",TST,".TST",".sscore"),row.names=F,col.names=T,quote=F,sep="\t")
  
  # prepare results
  i.max=which.max(auc.testing$auc)
  return(data.frame(TRN,REF,TST,auc=auc.testing$auc[i.max],s=unlist(strsplit(auc.testing$score[i.max],"_"))[1],
                    lambda=unlist(strsplit(auc.testing$score[i.max],"_"))[2],combn=gsub("_SUM","",auc.testing$score[i.max])))  
  
}


#### We take the results tuning of the hyper parameters to get the P&T PRS with the largest AIC
#### The TUNEing results are used to obtain the optimal weights
wlsTuning <- function(TUNE,WORK.DIR,tuning.parameters){
  ### this function find the optiomal weights bewteen the two TRaiNing sets to obtain the weighted P&T
  # TUNE - TUNEing dataset 
  # WORK.DIR - working directory # WORK.DIR=paste0(work.dir,"/",sim)
  # tuning.parameters - the parameters to use for the hyper parameters
  
  WLS.DIR=paste0(WORK.DIR,"/",TUNE,".TUNE/WLS-TUNING")
  dir.create(paste0(WLS.DIR),recursive = T, showWarnings = F)
  
  # obtain the the two sets of PRS for the selected hyper parameters
  PRS=list()
  for(TRN in tuning.parameters[ tuning.parameters$TUNE == TUNE,"TRN"]){
    TUNING.DIR=paste0(WORK.DIR,"/",TUNE,".TUNE/LS-TUNING/",TRN,".TRN/",TRN,".REF")
    
    sscore=fread(paste0(TUNING.DIR,"/",TRN,".TRN-",TRN,".REF","-TUNING.sscore"),header=T,data.table=F)
    combn=tuning.parameters[tuning.parameters$TRN == TRN & tuning.parameters$TUNE == TUNE,"combn"]
    i.col=which(colnames(sscore) == paste0(combn,"_SUM"))
    PRS[[TRN]]=sscore[,c(1:3,i.col)]
  }
  
  # find the optimal weight
  re.WGHT=NULL
  for(wght in seq(0,1,.01)){
    (auc=auc(PRS$CEU$PHENO1-1~eval(wght*scale(PRS$CEU[,4])+(1-wght)*scale(PRS$YRI[,4])))[1])
    re.WGHT=rbind.data.frame(re.WGHT,data.frame(CEU.wght=wght,YRI.wght=1-wght,auc))
  }
  
  fwrite(re.WGHT,paste0(WLS.DIR,"/WLS-TUNING.txt"),row.names=F,col.names=T,quote=F,sep="\t")
  return(data.frame(TUNE,re.WGHT[which.max(re.WGHT$auc),]))
}


##### Testing the weighted lassosum. Here the weights from the previous step are applied to our TeSTing set
wlsTesting <- function(TST,WORK.DIR,ls.weights){
  #### this function simpy applies the weights found from the TUNEing PRS to the TeSTing PRS
  # TST - testing dataset
  # WORK.DIR - Working directory
  # ls.weights - weights to apply to the different lassosum PRS
  
  WLS.DIR=paste0(WORK.DIR,"/",TST,".TST/WLS-TESTING")
  dir.create(paste0(WLS.DIR),recursive = T, showWarnings = F)
  
  # find the names of the TRN datasets
  rownames(ls.weights)=ls.weights$TUNE
  # change the column names for easy access
  colnames(ls.weights)=gsub(".wght","",colnames(ls.weights))
  # retain the columns with usefull information
  ls.weights=ls.weights[,2:(ncol(ls.weights)-1)]
  
  training.sets=colnames(ls.weights)
  
  
  # read the PRS sets
  LS=list()
  for(trn in training.sets){
    LS[[trn]]=fread(paste0(WORK.DIR,"/",TST,".TST/LS-TESTING/",trn,".TRN/",trn,".REF/",trn,".TRN-",trn,".REF-",TST,".TST.sscore"),header=T,data.table=F)
  }
  
  # determine the weighted PRS
  # general information
  WLS=LS[[1]]
  colnames(WLS)[ncol(WLS)]="SCORE1_SUM"
  rownames(WLS)=WLS$IID
  # collect the scoresum
  WLS$SCORE1_SUM=0
  for(TRN in training.sets){
    WLS[LS[[TRN]]$IID,"SCORE1_SUM"]=WLS[LS[[TRN]]$IID,"SCORE1_SUM"]+ls.weights[TST,TRN]*scale(LS[[TRN]][,ncol(LS[[TRN]])])
  }
  # determine the auc
  re.auc=auc(WLS$PHENO1-1~WLS$SCORE1_SUM)[1]
  # collect pertinent information
  WLS=WLS[,c("#FID","IID","PHENO1","SCORE1_SUM")]
  fwrite(WLS,paste0(WLS.DIR,"/WLS-",TST,".TST.sscore"),row.names=F,col.names=T,quote=F,sep="\t")
  
  return(data.frame(TST,CEU.wght=ls.weights[TST,"CEU"],YRI.wght=ls.weights[TST,"YRI"],auc=re.auc))
  
}
