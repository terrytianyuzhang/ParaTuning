rm(list=ls()); gc()
options(stringsAsFactors = F)

#### software needed
plink <- "/usr/local/bin/plink"
plink2 <- "/usr/local/bin/plink2"
setwd('/raid6/Tianyu/PRS/bert_sample/ReferencePopulation-Package')

# translate the reference datasets from plink 2.0 format to plink 1.9 format
sets=c("CEU-20K","YRI-4K")
for(set in sets){
  pdir=paste0(set,"/")
  dir.create(paste0(pdir,"CHR"),showWarnings = F)
  # complete genotype set from plink 2.0 format to plink 1.9
  # plink.command=paste(plink2,
  #                     "--pfile",paste0(pdir,set),
  #                     "--make-bed",
  #                     "--out",paste0(pdir,set))
  # system(plink.command)
  
  # now by chromosomes
  for(chr in 1:22){
    plink.command=paste(plink2,
                        "--pfile",paste0(pdir,set),
                        "--chr",chr,
                        "--make-bed",
                        "--out",paste0(pdir,"CHR","/",set,"-chr",chr))
    system(plink.command)
  }
}
