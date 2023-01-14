#rm(list=ls()); gc()
#options(stringsAsFactors = F)

library(data.table)
library(doParallel)
library(snpStats)
library(R.utils)
library(lassosum)

# #### software needed
# plink="/data3/Software/Plink/plink"
# plink2="/data3/Software/Plink2/plink2"
# directory changes
plink <- "/usr/local/bin/plink"
plink2 <- "/usr/local/bin/plink2"

#### load the parameters for the simulation
load(paste0("Work/Sim-",i.sim,"/simulation-params.RData"))
main.dir=params$run.info$main.dir
work.dir=params$run.info$work.dir


### split plink datasets by chromosome and create bfile formatfor CEU and YRI reference
#for(set in c("CEU","YRI")){
#  pdir=paste0("Data/Reference-LDblocks/",set,"/")
#  dir.create(paste0(pdir,"CHR"),showWarnings = F)
#  
#  for(chr in 1:22){
#    plink.command=paste(plink2,
#                        "--pfile",paste0(pdir,set),
#                        "--chr",chr,
#                        "--make-bed",
#                        "--out",paste0(pdir,"CHR","/",set,"-chr",chr))
#    system(plink.command)
#  }
#}

# tuning set
tune.sets=names(params)[grep(".TUNE",names(params))]
for(set in tune.sets){
  pdir=paste0(params$run.info$work.dir,set,"/")
  dir.create(paste0(pdir,"CHR"),showWarnings = F)
  
  for(chr in 1:22){
    plink.command=paste(plink2,
                        "--pfile",paste0(pdir,set),
                        "--chr",chr,
                        "--make-bed",
                        "--out",paste0(pdir,"CHR","/",set,"-chr",chr))
    system(plink.command)
  }
  plink.command=paste(plink2,
                      "--pfile",paste0(pdir,set),
                      "--make-bed",
                      "--out",paste0(pdir,set))
  system(plink.command)
}

# testing set
test.sets=names(params)[grep(".TST",names(params))]
for(set in test.sets){
pdir=paste0(params$run.info$work.dir,set,"/")
dir.create(paste0(pdir,"CHR"),showWarnings = F)

for(chr in 1:22){
  plink.command=paste(plink2,
                      "--pfile",paste0(pdir,set),
                      "--chr",chr,
                      "--make-bed",
                      "--out",paste0(pdir,"CHR","/",set,"-chr",chr))
  system(plink.command)
}
plink.command=paste(plink2,
                    "--pfile",paste0(pdir,set),
                    "--make-bed",
                    "--out",paste0(pdir,set))
system(plink.command)

}
# retain information that you want to carry over to another part of the simulation
rm.list=ls()
rm.list=rm.list[!(rm.list %in% c("i.sim","sims", "plink", "plink2"))]
rm(list = rm.list); flush.console()

