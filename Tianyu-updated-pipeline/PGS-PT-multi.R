#rm(list=ls()); gc()
options(stringsAsFactors = F)

library(data.table)
library(doParallel)
library(snpStats)
library(R.utils)
library(pROC)

#### software needed
plink="/data3/Software/Plink/plink"
plink2="/data3/Software/Plink2/plink2"


#### load the functions that are needed
source("simulation-functions.R")

#### load the parameters for the simulation
load(paste0("Work/Sim-",i.sim,"/simulation-params.RData"))
main.dir=params$run.info$main.dir
work.dir=params$run.info$work.dir

train.sets=names(params)[grep(".TRN",names(params))]
tune.sets=names(params)[grep(".TUNE",names(params))]
test.sets=names(params)[grep(".TST",names(params))]

# read the p-value cutoffs
p.cutoff=fread(paste0(main.dir,"Data/pval-range-list.txt"),header=F,data.table=F,colClasses = "character")
r2.list=c(0.1,0.2,0.5,0.8)

# Determine the optimal r2 and p-cutoff combination for each training set
re.AUC=NULL
for(train.set in train.sets){
  for(tune.set in tune.sets){
    train.dir=paste0(work.dir,train.set,"/")
    tune.dir=paste0(work.dir,tune.set,"/")
    
    for(r2 in r2.list){
      for(p in p.cutoff[,1]){
        prs=fread(paste0(tune.dir,"PRS/",train.set,"/",train.set,"-",tune.set,sprintf("-r2_%2.1f",r2),"-p.",p,".sscore"),header=T,data.table=F)
        if(!is.na(prs$SCORE1_AVG[1])){
          re.AUC=rbind.data.frame(re.AUC,data.frame(tune.set,train.set,r2,p,n.snp=prs$ALLELE_CT[1]/2,AUC=auc(prs$PHENO1-1~prs$SCORE1_AVG)[1]))
        }
        
      }
    }
  }
}

# find the max AUC for the training sets
df.train=NULL
par(mfrow=c(2,2))
for(tune.set in tune.sets){
  for(train.set in train.sets){
    auc.train=re.AUC[re.AUC$tune.set  == tune.set & re.AUC$train.set == train.set,]
    i.max=which.max(auc.train$AUC)
    df.train=rbind.data.frame(df.train,data.frame(tune.set,train.set,r2=auc.train$r2[i.max],p=auc.train$p[i.max],AUC=auc.train$AUC[i.max]))
    plot(auc.train$AUC)
  }
}
fwrite(df.train,paste0(work.dir,"multi-PT.r2_p"),row.names=F,col.names=T,quote=F,sep="\t")

# load the PRS for the best r2 and p combination for each dataset
par(mfrow=c(1,2))
re.wght=NULL
for(tune.set in tune.sets){
  PRS=list()
  for(train.set in train.sets){
    train.dir=paste0(work.dir,train.set,"/")
    tune.dir=paste0(work.dir,tune.set,"/")
    i.select=which(df.train$tune.set == tune.set & df.train$train.set == train.set)
    PRS[[train.set]]=fread(paste0(tune.dir,"PRS/",train.set,"/",train.set,"-",tune.set,sprintf("-r2_%2.1f",df.train[i.select,"r2"]),"-p.",df.train[i.select,"p"],".sscore"),header=T,data.table=F)
  }
    
  # obtain optimal weight between the two using the training set
  wghts=seq(0,1,0.05)
  
  AUC.multi=NULL
  for(wght in wghts){
    AUC.multi=rbind.data.frame(AUC.multi,data.frame(wght=wght,AUC=auc(PRS[[train.sets[1]]]$PHENO1-1~eval(wght*PRS[[train.sets[1]]]$SCORE1_AVG+(1-wght)*PRS[[train.sets[2]]]$SCORE1_AVG))))
  }
  
  plot(AUC.multi$wght,AUC.multi$AUC,main=tune.set)
  
  (max.wght=AUC.multi$wght[which.max(AUC.multi$AUC)])
  abline(v=max.wght)  
  re.wght=rbind.data.frame(re.wght,data.frame(anc=substr(tune.set,1,3),wght=max.wght))
  

  
  
}
fwrite(re.wght,paste0(work.dir,"multi-PT.wghts"),row.names=F,col.names=T,quote=F,sep="\t")



##################################################################################
#### now for the validation set calculate the PRS based on df.train and max.wght
re.p_r2=fread(paste0(work.dir,"multi-PT.r2_p"),header=T,data.table=F)
re.wght=fread(paste0(work.dir,"multi-PT.wghts"),header=T,data.table=F)

for(test.set in test.sets){
  test.dir=paste0(work.dir,test.set,"/")
  dir.create(paste0(test.dir,"PRS"),showWarnings = F)

  for(train.set in train.sets){
    dir.create(paste0(test.dir,"PRS/",train.set,"/"),showWarnings = F,recursive = T)
    train.dir=paste0(work.dir,train.set,"/")
    
    # testing genotype file
    test.pfile=paste0(test.dir,test.set)
    # training association data file
    assoc.file=paste(train.dir,"Assoc/",train.set,".PHENO1.glm.logistic.hybrid",sep="")
    
    ### read the training association results 
    assoc=fread(assoc.file,header=T,data.table=F)
    # determine the regression coefficient from the odds ratios
    assoc$B=log(assoc$OR)
    # write a short version with the essential information
    fwrite(assoc[,c(3,6,15)],paste0(test.dir,"PRS/",train.set,"/",train.set,".SNP.values"),row.names=F,col.names=T,quote=F,sep="\t")
    
    # read the clumping file
    i.select=which(re.p_r2$tune.set == gsub("TST","TUNE",test.set) & re.p_r2$train.set == train.set)
    clump.file=paste0(train.dir,"Clumping/",train.set,sprintf("-r2_%2.1f",df.train[i.select,"r2"]),".clumped")
    clumps=fread(clump.file,header=T,data.table=F)
    # select the set of SNP with the p-value requirement
    clumps=clumps[clumps$P < as.numeric(df.train[i.select,"p"]),]
    # write a file with the inde SNP id for each of the clumps
    write.table(as.data.frame(clumps$SNP),paste0(test.dir,"PRS/",train.set,"/",train.set,sprintf("-r2_%2.1f",df.train[i.select,"r2"]),".clumped.valid.snp"),row.names=F,col.names=F,quote=F,sep="\t")
    
    # use plink to calculate the PRS score
    plink.command=paste(plink2,"--allow-no-sex","--nonfounders",
                        "--pfile",test.pfile,
                        "--score",paste0(test.dir,"PRS/",train.set,"/",train.set,".SNP.values"),"header",
                        "--extract",paste0(test.dir,"PRS/",train.set,"/",train.set,sprintf("-r2_%2.1f",df.train[i.select,"r2"]),".clumped.valid.snp"),
                        "--out",paste0(test.dir,"PRS/",train.set,"/",test.set),
                        sep=" ")
    system(plink.command)
  
  
  }

  ## now combine the results using the weight
  # read the PRS
  test.dir=paste0(work.dir,test.set,"/")
  
  PRS=list()
  for(train.set in train.sets){
    PRS[[train.set]]=fread(paste0(test.dir,"PRS/",train.set,"/",test.set,".sscore"),header=T,data.table=F)
  }

  MULTI.PRS=cbind(PRS$CEU.TRN[,c(1,2,3,6)],PRS$YRI.TRN[,6])
  colnames(MULTI.PRS)[(ncol(MULTI.PRS)-1):ncol(MULTI.PRS)]=train.sets
  wght=re.wght[re.wght$anc == substr(test.set,1,3),"wght"]
  MULTI.PRS$MULTI=wght*PRS[[train.sets[1]]]$SCORE1_AVG+(1-wght)*PRS[[train.sets[2]]]$SCORE1_AVG
  fwrite(MULTI.PRS,paste0(test.dir,"PRS/",test.set,"-multi-PGS.txt"),row.names=F,col.names=T,quote=F,sep="\t")
  (auc(MULTI.PRS$PHENO1-1~MULTI.PRS$CEU.TRN)[1])
  (auc(MULTI.PRS$PHENO1-1~MULTI.PRS$YRI.TRN)[1])
  (auc(MULTI.PRS$PHENO1-1~MULTI.PRS$MULTI)[1])
  #  write.table(AUC,paste0(test.dir,"PRS/validation-AUC.txt"),row.names=F,col.names=T,quote=F,sep="\t")

}

# retain information that you want to carry over to another part of the simulation
rm.list=ls()
rm.list=rm.list[!(rm.list %in% c("i.sim"))]
rm(list = rm.list); flush.console()

  