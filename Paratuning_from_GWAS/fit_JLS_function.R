
#### Functions to use
#Joint Lassosum Function Peng created
mylassosumFunction <- function(chr,gamma,lambda,shrink,
                               mem.limit,gwasANC,COR,
                               referenceFiles,LDblocks){
  re=mylassosum(cor1=COR[COR$CHR == chr,gwasANC[1]], 
                cor2=COR[COR$CHR == chr,gwasANC[2]],
                fileName1 = referenceFiles[[chr]][[gwasANC[1]]], 
                fileName2 = referenceFiles[[chr]][[gwasANC[2]]],  
                gamma = gamma, lambda = lambda, shrink=shrink,
                chunk=TRUE, mem.limit=mem.limit,
                trace=2,
                LDblocks1=as.data.frame(LDblocks[[gwasANC[1]]][LDblocks[[gwasANC[1]]]$chr == paste0("chr",chr),]), 
                LDblocks2=as.data.frame(LDblocks[[gwasANC[2]]][LDblocks[[gwasANC[2]]]$chr == paste0("chr",chr),]))
  return(re)
}

wrapperFunction <- function(i.combn, input.df, gwasANC, lambda, shrink, main.dir, work.dir, CHR=1:22, 
                            mem.limit, mymc.cores = 8){
  
  # collect gamma
  gamma=input.df$gamma[i.combn]  
  
  # collectthe number of samples for the two populations
  gwasN=list()
  gwasN[["CEU"]]=input.df$N1[i.combn]
  gwasN[["YRI"]]=input.df$N2[i.combn]
  
  # > gwasN
  # $CEU
  # [1] 20000
  # 
  # $YRI
  # [1] 20000
  
  
  #### determine correlations
  ### Summary statistics ###
  # read a map for the first ancestry to set up a dataframe
  anc=gwasANC[1]
  map <- fread(paste0(work.dir, 'TST/', anc, '.TST.pvar'), header=T, data.table=F)
  # map=fread(paste0(work.dir,anc,".GWAS/",anc,".GWAS-",gwasN[[anc]],".pvar"),header=T,data.table=F)
  # > head(map)
  # #CHROM     POS            ID REF ALT
  # 1      1 1962845 1:1962845:T:C   T   C
  # 2      1 1962899 1:1962899:A:C   A   C
  # 3      1 1963406 1:1963406:G:A   G   A
  # 4      1 1963538 1:1963538:T:C   T   C
  # 5      1 1963738 1:1963738:C:T   C   T
  # 6      1 1964101 1:1964101:A:G   A   G
  
  
  COR=data.frame(CHR=map$`#CHROM`,ID=map$ID)
  
  # ##> head(COR)
  # CHR            ID
  # 1   1 1:1962845:T:C
  # 2   1 1:1962899:A:C
  # 3   1 1:1963406:G:A
  # 4   1 1:1963538:T:C
  # 5   1 1:1963738:C:T
  # 6   1 1:1964101:A:G
  
  for(anc in gwasANC){
    # glm=fread(paste0(work.dir,anc,".GWAS/Assoc/",anc,".GWAS-",gwasN[[anc]],".PHENO1.glm.logistic.hybrid"),header=T,data.table=F)
    # glm <- fread(paste0(work.dir, anc, '.TRN.PHENO1.glm.logistic.hybrid'), header=T, data.table=F)
    glm <- fread(paste0(work.dir, 'TRN/', anc, '.TRN.PHENO1.glm.logistic.hybrid'),
                 header=T,data.table=F)
    #'/raid6/Ron/prs/data/bert_sample/CEU.TRN.PHENO1.glm.logistic.hybrid'
    COR[,anc]=p2cor(p = glm$P, n = gwasN[[anc]], sign=log(glm$OR))
  }
  rownames(COR)=COR$ID
  
  # head(COR)
  # 
  # > head(COR)
  # CHR            ID          CEU          YRI
  # 1:1962845:T:C   1 1:1962845:T:C -0.008296110 0.0073444745
  # 1:1962899:A:C   1 1:1962899:A:C  0.013275142 0.0086428250
  # 1:1963406:G:A   1 1:1963406:G:A -0.006433539 0.0007753638
  # 1:1963538:T:C   1 1:1963538:T:C  0.004405373 0.0137580552
  # 1:1963738:C:T   1 1:1963738:C:T -0.006433539 0.0007753638
  # 1:1964101:A:G   1 1:1964101:A:G -0.006433539 0.0007753638
  
  # load the LD block information
  LDblocks=list()
  for(anc in gwasANC){
    if(anc == "CEU"){
      ld.anc="EUR"
    }else if(anc == "YRI"){
      ld.anc="AFR"
    }
    LDblocks[[anc]]=read.table2(system.file(paste0("data/Berisa.", 
                                                   paste0(ld.anc,".hg38"), ".bed"), 
                                            package="lassosum"), header=T)
    #  LDblocks[[anc]]$chr=gsub("chr","",LDblocks[[anc]]$chr)
  }
  
  # > head(LDblocks)
  # $CEU
  # chr     start      stop
  # 1     chr1   1961168   3666172
  # 2     chr1   3666172   4320751
  # 3     chr1   4320751   5853833
  # 4     chr1   5853833   7187275
  # 5     chr1   7187275   9305140
  # 6     chr1   9305140  10746927
  # 7     chr1  10746927  11717784
  # 8     chr1  11717784  12719464
  # 9     chr1  12719464  14565015
  
  
  # reference files by chromosome
  referenceFiles=list()
  for(chr in 1:22){
    referenceFiles[[chr]]=list()
    for(anc in gwasANC){
      # referenceFiles[[chr]][[anc]]=paste0("/data3/DownLoadedData/GWAS-Populations-SimulationInput/",anc,"_reference_LDblocks/CHR/",anc,"-chr",chr)
      referenceFiles[[chr]][[anc]] <- paste0(main.dir,"Data/Reference-LDblocks/", anc,"/CHR/", anc,"-chr",chr)
      # paste0(work.dir, "GWAS-Populations-SimulationInput/", anc, "_reference_LDblocks/CHR/", anc,"-chr",chr)
    }
  }
  
  # > referenceFiles
  # [[1]]
  # [[1]]$CEU
  # [1] "/raid6/Ron/prs/data/bert_sample/GWAS-Populations-SimulationInput/CEU_reference_LDblocks/CHR/CEU-chr1"
  # 
  # [[1]]$YRI
  # [1] "/raid6/Ron/prs/data/bert_sample/GWAS-Populations-SimulationInput/YRI_reference_LDblocks/CHR/YRI-chr1"
  # 
  # 
  # [[2]]
  # [[2]]$CEU
  # [1] "/raid6/Ron/prs/data/bert_sample/GWAS-Populations-SimulationInput/CEU_reference_LDblocks/CHR/CEU-chr2"
  # 
  # run lassosum
  print(CHR); flush.console()
  
  system.time(
    re.chr<-mclapply(CHR, mylassosumFunction,
                     gamma=gamma,lambda=lambda,
                     shrink=shrink,mem.limit = mem.limit,
                     gwasANC = gwasANC,COR=COR,
                     referenceFiles = referenceFiles,
                     LDblocks = LDblocks,
                     mc.cores = mymc.cores,
                     mc.preschedule = F)
  )
  # print('loss in combined lassosum')
  # print(re.chr$loss)
  # print(re.chr$trainerror1)
  # print(re.chr$trainerror2)
  
  # merge the results from the 22 chromosomes
  re.lasso = merge.mylassosum(re.chr)
  print('loss in combined lassosum, after merging')
  print(re.lasso$loss)
  print(re.lasso$trainerror1)
  print(re.lasso$trainerror2)
  # save ther results
  dir.create(paste0(work.dir,"JointLassoSum/"),showWarnings = F,recursive = T)
  save(re.lasso,file=paste(work.dir,"JointLassoSum/JointLassosum-",sprintf("-gamma-%.2f",gamma),".Rdata",sep=""))
  return(i.combn)
}
