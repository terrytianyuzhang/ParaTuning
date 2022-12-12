#rm(list=ls()); 
gc()
options(stringsAsFactors = F)

library(data.table)
library(doParallel)
library(snpStats)
library(R.utils)
library(lassosum)
library(pROC)

#### software needed
# plink="/data3/Software/Plink/plink"
# plink2="/data3/Software/Plink2/plink2"
if(!exists("plink") | !exists("plink2")){
  plink <- "/usr/local/bin/plink"
  plink2 <- "/usr/local/bin/plink2"
}

#### load the functions that are needed
source("simulation-functions.R")

# load the parameters for this simulation
load(paste0("Work/Sim-",i.sim,"/simulation-params.RData"))
main.dir=params$run.info$main.dir
work.dir=params$run.info$work.dir

train.sets=names(params)[grep(".TRN",names(params))]
tune.sets=names(params)[grep(".TUNE",names(params))]
test.sets=names(params)[grep(".TST",names(params))]

re.wght=NULL
par(mfrow=c(1,2))
for(tune.set in tune.sets){
  PRS=list()
  for(train.set in train.sets){
    train.dir=paste0(work.dir,train.set,"/")
    tune.dir=paste0(work.dir,tune.set,"/")
    PRS[[train.set]]=fread(paste0(tune.dir,"Lassosum/",tune.set,"/",train.set,"-",tune.set,"-lassosum.score"),header=T,data.table=F)
    PRS[[train.set]]$best.pgs=as.vector(scale(PRS[[train.set]]$best.pgs))
  }
  
  # obtain optimal weight between the two using the training set
  wghts=seq(0,1,0.05)
  
  AUC.weighted=NULL
  for(wght in wghts){
    AUC.weighted=rbind.data.frame(AUC.weighted,data.frame(wght=wght,AUC=auc(PRS[[train.sets[1]]]$pheno~eval(wght*PRS[[train.sets[1]]]$best.pgs+(1-wght)*PRS[[train.sets[2]]]$best.pgs))))
  }
  
  plot(AUC.weighted$wght,AUC.weighted$AUC,main=tune.set)
  
  (max.wght=AUC.weighted$wght[which.max(AUC.weighted$AUC)])
  abline(v=max.wght)  
  
  # now for the testing set
  test.set=paste0(substr(tune.set,1,3),".TST")
  PRS=list()
  for(train.set in train.sets){
    train.dir=paste0(work.dir,train.set,"/")
    test.dir=paste0(work.dir,test.set,"/")
    PRS[[train.set]]=fread(paste0(test.dir,"Lassosum/",test.set,"/",train.set,"-",test.set,"-lassosum.score"),header=T,data.table=F)
    PRS[[train.set]]$best.pgs=as.vector(scale(PRS[[train.set]]$best.pgs))
  }
  
  wght.PGS=cbind.data.frame(PRS[[train.sets[1]]][,1:3],PRS[[train.sets[1]]][,"best.pgs"],PRS[[train.sets[[2]]]][,"best.pgs"])
  colnames(wght.PGS)[(ncol(wght.PGS)-1):ncol(wght.PGS)]=paste0(train.sets,".scaled")
                            
  # calculate the weighted PGS
  wght.PGS=cbind.data.frame(wght.PGS,data.frame(WEIGHTED.scaled=scale(max.wght*wght.PGS[,paste0(train.sets[1],".scaled")]+
                                                                        (1-max.wght)*wght.PGS[,paste0(train.sets[2],".scaled")])))

  (auc(wght.PGS$pheno~wght.PGS$CEU.TRN.scaled))
  (auc(wght.PGS$pheno~wght.PGS$YRI.TRN.scaled))
  (auc(wght.PGS$pheno~wght.PGS$WEIGHTED.scaled))
  
  fwrite(wght.PGS,paste0(test.dir,"Lassosum/",test.set,"/","WEIGHTED","-",test.set,"-lassosum.score"),row.names=F,col.names=T,quote=F,sep="\t")
  
  re.wght=rbind.data.frame(re.wght,data.frame(anc=substr(tune.set,1,3),wght=max.wght))
  
}

fwrite(re.wght,paste0(work.dir,"weighted-lassosum.wghts"),row.names=F,col.names=T,quote=F,sep="\t")
print(re.wght)

print('finished PGS-lassosum-weighted')

rm.list=ls()
rm.list=rm.list[!(rm.list %in% c("i.sim","sims"))]
rm(list = rm.list); flush.console()
