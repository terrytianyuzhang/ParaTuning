options(stringsAsFactors = F)

library(data.table)
library(doParallel)
library(snpStats)
library(R.utils)

#### software needed
plink="/data3/Software/Plink/plink"
plink2="/data3/Software/Plink2/plink2"



#### load the functions that are needed
source("R-code/simulation-functions.R")

#### load the parameters for the simulation
i.sim="NC8000"
load(paste0("Work/Sim-C0.80-Y0.60-rep",i.sim,"/simulation-params.RData"))


# create different subsets of samples
work.dir=params$run.info$work.dir
sets=names(params)[c(grep(".GWAS",names(params)),grep(".TRN",names(params)))]
# subset size
N=seq(20000,20000,4000)
for(set in sets){
  fam=fread(paste0(work.dir,set,"/",set,".psam"),header=T,data.table=F)
  rownames(fam)=fam$IID
  for(n in N){
    cases=sample(rownames(fam)[fam$PHENO == 2],n/2)
    controls=sample(rownames(fam)[fam$PHENO == 1],n/2)
    fwrite(fam[c(cases,controls),1:2],paste0(work.dir,set,"/",set,"-",n,".ind"),row.names=F,col.names=F,quote=F,sep="\t")
    
    plink.command=paste(plink2,
                        "--pfile",paste0(work.dir,set,"/",set),
                        "--keep",paste0(work.dir,set,"/",set,"-",n,".ind"),
                        "--make-pgen",
                        "--out",paste0(work.dir,set,"/",set,"-",n),
                        sep=" ")
    system(plink.command)
  }
}

# run GWAS on each of the subsets
gwas.sets=names(params)[grep(".GWAS",names(params))]
for(set in gwas.sets){
  for(n in N){
    set.dir=paste0(work.dir,set,"/")
    dir.create(paste0(set.dir,"Assoc/"),showWarnings = F)
    
    plink2.command=paste(plink2,"--nonfounders","--allow-no-sex",
                         "--pfile",paste0(set.dir,set,"-",n),
                         "--glm","allow-no-covars","omit-ref",
                         "--out",paste0(set.dir,"Assoc/",set,"-",n),
                         sep=" ")
    system(plink2.command)
    
  }
  
}
  
# translate pfile into bfile so it can be used by lassosum
for(set in sets){
  for(n in N){
    plink.command=paste(plink2,
                        "--pfile",paste0(work.dir,set,"/",set,"-",n),
                        "--make-bed",
                        "--out",paste0(work.dir,set,"/",set,"-",n),
                        sep=" ")
    system(plink.command)
  }
}

### the training sets need to be split by chr
trn.sets=names(params)[c(grep(".TRN",names(params)))]

for(set in trn.sets){
  for(n in N){
    for(chr in 1:22){
      plink.command=paste(plink2,
                          "--pfile",paste0(work.dir,set,"/",set,"-",n),
                          "--chr",chr,
                          "--make-bed",
                          "--out",paste0(work.dir,set,"/",set,"-",n,"-chr",chr),
                          sep=" ")
      system(plink.command)
    }
  }
}




