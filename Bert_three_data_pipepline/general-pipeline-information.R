# this piece of R-code contains essential information that is being used across the analyses

#### location of the main and work directories
main.dir="/data3/Bert/PengLiu/Tianyu"          # location where the code is run from
data.dir=paste0(main.dir,"/Data")              # this is the root for the location where the original genotype files are located
work.dir=paste0(main.dir,"/Work")              # root of location where all the work files are written to
ref.dir=paste0(data.dir,"/REF.GNT")            # directoiry with information on the reference populations

# location of plink executables
plink <- "/data3/Software/Plink/plink"
plink2 <- "/data3/Software/Plink2/plink2"

# libraries that are used by the pipeline
library(data.table)
library(doParallel)
library(snpStats)
library(R.utils)
library(lassosum)
library(pROC)


#### GENOTYPE DATA SETS and their locations
GNT.SETS=list()
# training sets and their location 
# in the simulation setting we replace XXX with the simulation number/name
# these provide the GWAS information
GNT.SETS[["TRN"]]=list()
GNT.SETS[["TRN"]][["CEU"]]=paste0(data.dir,"/XXX/CEU.TRN")
GNT.SETS[["TRN"]][["YRI"]]=paste0(data.dir,"/XXX/YRI.TRN")

# tuning sets
# these are used when (hyper)-parameters need to be chosen
GNT.SETS[["TUNE"]]=list()
GNT.SETS[["TUNE"]][["CEU"]]=paste0(data.dir,"/XXX/CEU.TUNE")
GNT.SETS[["TUNE"]][["YRI"]]=paste0(data.dir,"/XXX/YRI.TUNE")

# testing sets
# these are used when we test 
GNT.SETS[["TST"]]=list()
GNT.SETS[["TST"]][["CEU"]]=paste0(data.dir,"/XXX/CEU.TST")
GNT.SETS[["TST"]][["YRI"]]=paste0(data.dir,"/XXX/YRI.TST")

# reference datasets
GNT.SETS[["REF"]]=list()
GNT.SETS[["REF"]][["CEU"]]="/data3/DownLoadedData/GWAS-Populations-SimulationInput/Reference-C20K-Y4K/CEU.REF"
GNT.SETS[["REF"]][["YRI"]]="/data3/DownLoadedData/GWAS-Populations-SimulationInput/Reference-C20K-Y4K/YRI.REF"

# define the LD informationm to use for the two populations
ld.populations=c("EUR.hg38","AFR.hg38"); names(ld.populations)=c("CEU","YRI")

# information on the number of samples
N=list()
N[["CEU"]]=data.frame(n.case=10000,n.control=10000)
N[["YRI"]]=data.frame(n.case=2000,n.control=2000)

# lassosum parameters
s=seq(0.8,1.0,.05)        # shrinkage parameter
lambda=abs(exp(seq(log(.02), log(0.001), length.out=20))-.021) # tuning parameter, this can be chosen


