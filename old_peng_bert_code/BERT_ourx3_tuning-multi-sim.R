#rm(list=ls()); 
gc()
options(stringsAsFactors = F)

#### software needed
plink="/data3/Software/Plink/plink"
plink2="/data3/Software/Plink2/plink2"

#### load the functions that are needed
source("R-code/simulation-functions.R")

### Load functions ###
source(paste0("R-code/LassosumFunctions/parseselect.R"))
source(paste0("R-code/LassosumFunctions/parseblocks.R"))
source(paste0("R-code/LassosumFunctions/ncol.bfile.R"))
source(paste0("R-code/LassosumFunctions/nrow.bfile.R"))
source(paste0("R-code/LassosumFunctions/read.table2.R"))
source(paste0("R-code/LassosumFunctions/selectregion.R"))
source(paste0("R-code/LassosumFunctions/parse.pheno.covar.R"))
source(paste0("R-code/LassosumFunctions/myelnet.R"))
source(paste0("R-code/LassosumFunctions/mylassosum.R"))
source(paste0("R-code/LassosumFunctions/splitgenome.R"))
source(paste0("R-code/LassosumFunctions/validation.R"))
source(paste0("R-code/LassosumFunctions/merge.mylassosum.R"))
Rcpp::sourceCpp(paste0("R-code/LassosumFunctions/myfunctions.cpp"))

library(data.table)
library(pryr) # check memory useage
library(lassosum) #transform p value to correlation
library(doParallel) # foreach
library(pROC) # for AUC 
library(R.utils)


#wrap pipeline into a function
predictPGS.blk <- function(i.result, shrink, results, chrom, anc1, anc2, cor1, cor2, N1, N2, snpCHROM, LDblocks1, LDblocks2,
                           root1, root2, tune.set, mem.limit = 4*e9){
  print(root2)
  print(tune.set)
  print(file.exists(paste0(root2,tune.set,"/",tune.set,".bed")))
  # set the parameters
  gamma=results$gamma[i.result]
  lambda=results$lambda[i.result]
  strt.time=Sys.time()
  
  extract.n=0
  df=0
  beta=NULL
  pgscore=NULL
  pred1=0
  pred2=0
  loss=0
  for(i in chrom){
    print(paste0("chr",i))
    #fit model for each chromosome seperately
    chr=paste0("chr",i)
    
    myfit <- mylassosum(cor1[snpCHROM==chr], cor2[snpCHROM==chr],
                        fileName1 = paste0(root1,"Data/Reference-LDblocks/",anc1,"/CHR/",anc1,"-",chr), 
                        fileName2 = paste0(root1,"Data/Reference-LDblocks/",anc2,"/CHR/",anc2,"-",chr),  
                        gamma = gamma, lambda = lambda, shrink=shrink,
                        chunk=TRUE, mem.limit=mem.limit,
                        trace=2,
                        LDblocks1=LDblocks1[LDblocks1$chr == chr,], LDblocks2=LDblocks2[LDblocks2$chr == chr,])
    
    # collect the results
    extract.n=extract.n + nrow(myfit$beta)
    df=df+ sum(myfit$beta!=0)
    loss = loss + myfit$loss
    beta = rbind(beta, myfit$beta)
    
    pred1 = pred1 + myfit$pred1
    pred2 = pred2 + myfit$pred2
    #validate for each chromosome seperately
    validate.result <- validate(myfit, paste0(root2,tune.set,"/CHR/",tune.set,"-",chr), extract=NULL, keep=NULL, distance=NULL, mem.limit=mem.limit)
    #validate.result.list[[i]]<-validate.result
    if(is.null(pgscore)) {
      pgscore = validate.result$pgscore
    } else {
      pgscore = pgscore + validate.result$pgscore
    }
  }
  print("Then calculate R2 and AUC")
  # ----- calculate r2 and AUC -----# 
  parsed<-parseselect(paste0(root2,tune.set,"/",tune.set), extract=NULL, keep=NULL, distance=NULL)
#  parsed<-parseselect(paste0(root2,tune.set,"/CHR/",tune.set,"-",chr), extract=NULL, keep=NULL, distance=NULL)
  phcovar <- parse.pheno.covar(pheno=NULL, covar=NULL, parsed=parsed, trace=2)
  pheno=phcovar$pheno
  r2=cor(as.numeric(pgscore), pheno, use="complete.obs") #remove missing values
  (AUC=as.numeric(ci.auc(roc(pheno, as.numeric(pgscore))))[2])
  AUC.up=as.numeric(ci.auc(roc(pheno, as.numeric(pgscore))))[3]
  AUC.dn=as.numeric(ci.auc(roc(pheno, as.numeric(pgscore))))[1]
  
  print("Then calculate AIC BIC")
  x11=sum(pred1^2)
  x12=t(beta) %*% cor1
  x21=sum(pred2^2)
  x22=t(beta) %*% cor2
  sig2_1=1/(N1-1)
  sig2_2=1/(N2-1)
  l1=1 + 10*sum(pred1^2) - 2 * t(beta) %*% cor1
  l2=1 + 10*sum(pred2^2) - 2 * t(beta) %*% cor2
  AIC = gamma*(l1/(sig2_1) + 2/N1*df) + (1-gamma)*(l2/(sig2_2) + 2/N2*df)
  BIC = gamma*(l1/(sig2_1) + (log(N1)/N1)*df) + (1-gamma)*(l2/(sig2_2) + (log(N2)/N2)*df)
  l11=1 + sum(pred1^2) - 2 * t(beta) %*% cor1
  l22=1 + sum(pred2^2) - 2 * t(beta) %*% cor2
  AICc = gamma*(l11/(sig2_1) + 2/N1*df) + (1-gamma)*(l22/(sig2_2) + 2/N2*df)
  BICc = gamma*(l11/(sig2_1) + (log(N1)/N1)*df) + (1-gamma)*(l22/(sig2_2) + (log(N2)/N2)*df)
  end.time=Sys.time()
  run.time=end.time-strt.time
  stats= data.frame(gamma = gamma, lambda = lambda, extract.n=extract.n, df=df, r2=r2, AUC=AUC, AUC.up=AUC.up, AUC.dn=AUC.dn, 
                    AIC=AIC, BIC=BIC, AICc=AICc, BICc=BICc, l1=l1, l2=l2, x11=x11, x12 =x12, x21=x21, x22=x22, run.time=run.time)
  oneres = list(beta=beta, pheno=pheno,pgs=pgscore,stats=stats
                #,fit=myfit.list,validate=validate.result.list
  )
  return(oneres)
}

wrapperFunction<- function(i.sim, max.cores=30){
  source("R-code/simulation-functions.R")
  
  ### Load functions ###
  source(paste0("R-code/LassosumFunctions/parseselect.R"))
  source(paste0("R-code/LassosumFunctions/parseblocks.R"))
  source(paste0("R-code/LassosumFunctions/ncol.bfile.R"))
  source(paste0("R-code/LassosumFunctions/nrow.bfile.R"))
  source(paste0("R-code/LassosumFunctions/read.table2.R"))
  source(paste0("R-code/LassosumFunctions/selectregion.R"))
  source(paste0("R-code/LassosumFunctions/parse.pheno.covar.R"))
  source(paste0("R-code/LassosumFunctions/myelnet.R"))
  source(paste0("R-code/LassosumFunctions/mylassosum.R"))
  source(paste0("R-code/LassosumFunctions/splitgenome.R"))
  source(paste0("R-code/LassosumFunctions/validation.R"))
  source(paste0("R-code/LassosumFunctions/merge.mylassosum.R"))
  Rcpp::sourceCpp(paste0("R-code/LassosumFunctions/myfunctions.cpp"))
  
  
  # load the parameters for this simulation
  load(paste0("Work/Sim-",i.sim,"/simulation-params.RData"))
  
  # xet directories for this simulation
  root1=params$run.info$main.dir
  root2=params$run.info$work.dir
  

  ### Parameter panel
  chr.n=1:22
  anc1="CEU"
  anc2="YRI"
  ### Tuning parameters ###
  shrink=0.9
  gamma= seq(1,0,by=-0.1) #; gamma=gamma[c(8,9)]
  lambda = seq(0.0001,0.1, length.out=10) #; lambda=lambda[c(8,9)]
  results = data.frame(gamma=rep(gamma, each = length(lambda)),  lambda=rep(lambda, length(gamma)))
  N1=params[[paste0(anc1,".TRN")]]$n.case+params[[paste0(anc1,".TRN")]]$n.control
  N2=params[[paste0(anc2,".TRN")]]$n.case+params[[paste0(anc2,".TRN")]]$n.control
  print(paste0("N1:",N1,", N2:",N2))
  
  ### Read information for SNPs ###
  snp=fread(paste0(root1,"Data/chr1-22-qc-frq-ld.block.map"),header=T,data.table=F)
  snpCHROM=snp$CHROM
  
  ### Read LD region file ###
  read.table2 <- function(file, header=F, data.table=F, check.names=TRUE, ...) {
    return(data.table::fread(file, header=header, data.table=data.table, 
                             check.names=check.names, ...))
  }
  
  LDblocks <- "EUR.hg38" # This will use LD regions as defined in Berisa and Pickrell (2015) for the European population and the hg19 genome.
  LDblocks_EUR <- read.table2(system.file(paste0("data/Berisa.", 
                                                 LDblocks, ".bed"), 
                                          package="lassosum"), header=T)
  
  LDblocks<- "AFR.hg38" # This will use LD regions as defined in Berisa and Pickrell (2015) for the European population and the hg19 genome.
  LDblocks_AFR <- read.table2(system.file(paste0("data/Berisa.", 
                                                 LDblocks, ".bed"), 
                                          package="lassosum"), header=T)
  
  ### Summary statistics ###
  glm1=fread(paste(root2,anc1,".TRN/Assoc/",anc1,".TRN.PHENO1.glm.logistic.hybrid",sep=""),header=T,data.table=F)
  glm2=fread(paste(root2,anc2,".TRN/Assoc/",anc2,".TRN.PHENO1.glm.logistic.hybrid",sep=""),header=T,data.table=F)
  
  # Transform p value to correlation 
  cor1 <- p2cor(p = glm1$P, n = N1, sign=log(glm1$OR)) # n: sample size in total?
  cor2 <- p2cor(p = glm2$P, n = N2, sign=log(glm2$OR))
  rm(glm1,glm2, snp); gc()
  
#  cor1=cor1[snpCHROM == "chr22"]
#  cor2=cor2[snpCHROM == "chr22"]
#  snpCHROM=snpCHROM[snpCHROM == "chr22"]
  ### Let's go! ###
#  results=results[105:110,]
  for(tune.set in c("CEU.TUNE","YRI.TUNE")){
  
    system.time(results.all.blk <- mclapply(1:nrow(results),predictPGS.blk, shrink = shrink, results = results, chrom = 1:22, anc1 = anc1, anc2 = anc2, cor1 = cor1, cor2 = cor2, N1 = N1, N2 = N2,
                                            snpCHROM = snpCHROM, LDblocks1 = LDblocks_EUR, LDblocks2 = LDblocks_AFR, 
                   root1 = root1, root2 = root2, tune.set = tune.set,mem.limit=12e9, mc.cores = min(nrow(results),max.cores), mc.preschedule = F, mc.silent = T))
  
    dir.create(paste0(root2,tune.set,"/Lassosum/",tune.set),showWarnings=F,recursive=T)
    save(results.all.blk, file=paste0(root2,tune.set,"/Lassosum/MULTI.TRN-lassosum-",tune.set,".test.RData"))
    
    # find the best combination of gamma and lambda
    RE=NULL
    for(i in 1:length(results.all.blk)){
      RE=rbind.data.frame(RE,results.all.blk[[i]]$stats)
    }
    sel.lasso.results=results.all.blk[[which.max(RE$AUC)]]
    save(sel.lasso.results,file=paste0(root2,tune.set,"/Lassosum/MULTI.TRN-lassosum-BEST-AUC-",tune.set,".test.RData"))
    
    rm(results.all.blk); gc()
  }
  return("done")
}

#wrapperFunction(600)


system.time(x<-mclapply(sims,wrapperFunction,max.cores=30,mc.cores=5,mc.preschedule = F,mc.silent=T))

# retain information that you want to carry over to another part of the simulation
rm.list=ls()
rm.list=rm.list[!(rm.list %in% c("i.sim","sims"))]
rm(list = rm.list); flush.console()
