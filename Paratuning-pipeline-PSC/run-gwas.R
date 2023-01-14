#rm(list=ls()); gc()
options(stringsAsFactors = F)

library(data.table)
library(doParallel)
library(snpStats)
library(R.utils)

# #### software needed
# plink="/data3/Software/Plink/plink"
# plink2="/data3/Software/Plink2/plink2"
# directory changes
if(!exists("plink") | !exists("plink2")){
  plink <- "/usr/local/bin/plink"
  plink2 <- "/usr/local/bin/plink2"
}

#### load the functions that are needed
source("simulation-functions.R")

#### load the parameters for the simulation
load(paste0("Work/Sim-",i.sim,"/simulation-params.RData"))
main.dir=params$run.info$main.dir
work.dir=params$run.info$work.dir

train.sets=names(params)[grep(".TRN",names(params))]
for(set in train.sets){

  set.dir=paste0(work.dir,set,"/")
  dir.create(paste0(set.dir,"Assoc/"),showWarnings = F)
  
  plink2.command=paste(plink2,"--nonfounders","--allow-no-sex",
                       "--pfile",paste(set.dir,set,sep=""),
                       "--glm","allow-no-covars","omit-ref",
                       "--out",paste(set.dir,"Assoc/",set,sep=""),
                       sep=" ")
  system(plink2.command)
  
}
  
# retain information that you want to carry over to another part of the simulation
rm.list=ls()
rm.list=rm.list[!(rm.list %in% c("i.sim","sims"))]
rm(list = rm.list); flush.console()


#gwas=fread(paste0("Work/Sim-",i.sim,"/CEU.TRN/Assoc/CEU.TRN.PHENO1.glm.logistic.hybrid"),header=T,data.table=F)
#i=match(params$run.info$sel.snp$ID,gwas$ID)
#plot(params$run.info$sel.snp$CEU.beta,log(gwas[i,"OR"]))
#abline(a=0,b=1)

#sel.snp=params$run.info$sel.snp
#sel.snp[,c("beta","p","a1")]=gwas[i,c("OR","P","A1")]
#sel.snp$beta=log(sel.snp$beta)
#sel.snp$beta[sel.snp$riskAllele == "REF"]=-sel.snp$beta[sel.snp$riskAllele == "REF"]
#summary(lm(sel.snp$beta~sel.snp$CEU.beta))#

#n.snp=table(sel.snp$CEU.blk)
#sel.snp$n.snp=n.snp[match(sel.snp$CEU.blk,names(n.snp))]
