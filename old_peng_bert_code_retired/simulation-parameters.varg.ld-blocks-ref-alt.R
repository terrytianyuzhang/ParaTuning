#rm(list=ls()); gc()
options(stringsAsFactors = F)

# this version of the simulation sets the heritability and uses this to determine the allele effect

library(data.table)
library(doParallel)
library(snpStats)
library(R.utils)

# simulation and seed for this run; seed needs to change from simulation to simulation
# seed=918s

# number of causal SNP 
n.causal=4000 ###tianyu don't have to specify it when just creating people
i.sim=paste0("NC",n.causal)

# set the total SNP heritability
h2=data.frame(CEU = 0.8,YRI = .6)
vare=1 # needs to be 1 for logistic regression, environment variance
varg=(h2/(1-h2))*vare

# set the random number generator using seed
seed=sample(1e8,1)

# store the information
params=list()
params[["run.info"]]$sim=i.sim
params[["run.info"]]$seed=seed
params[["run.info"]]$n.causal=n.causal
params[["run.info"]]$varg=varg
params[["run.info"]]$vare=vare

# two versions of plink that might be needed
plink="/data3/Software/Plink/plink"
plink2="/data3/Software/Plink2/plink2"

# directories used, create the workdir if needed
main.dir="/data3/Bert/PengLiu/SimulationMay2022/"
work.dir=paste0(main.dir,"Work/Sim-",sprintf("C%.2f-Y%.2f",h2$CEU,h2$YRI),"-rep",i.sim,"/")
dir.create(work.dir,showWarnings = F, recursive = T)
params$run.info$main.dir=main.dir
params$run.info$work.dir=work.dir

# read the map for our analysis, tianyu may not need this
snp.map=fread("/data3/DownLoadedData/GWAS-Populations-SimulationInput/chr1-22-qc-frq-ld.block.map",header=T,data.table=F)
snp.map$CHROM=as.numeric(gsub("chr","",snp.map$CHROM))
rownames(snp.map)=snp.map$ID

# use equal spacing of causal snp
#find the length of each chromosome
df=aggregate(snp.map$POS~snp.map$CHROM,FUN="range")
df=cbind.data.frame(df$`snp.map$CHROM`,df$`snp.map$POS`)
colnames(df)=c("chr","start","stop")
df$length=df$stop-df$start
total.length=sum(df$length)

sel.chr=sample(1:22,n.causal+200,replace=T,prob=df$length/sum(df$length))
df$n.causal=as.vector(table(sel.chr))

re=NULL
selected.snp=NULL
for(chr in 1:22){
  print(chr)
  locs=seq(df$start[chr]+runif(1,min=100e3,max=250e3),df$stop[chr],ceiling(df$length[chr]/df$n.causal[chr]))
  if(max(locs) > df$stop[chr] | length(locs) != df$n.causal[chr]){
    print(chr); flush.console()
    stop
  }
  locs=c(locs,df$stop[chr])
  tmp=as.numeric(cut(snp.map[snp.map$CHROM == chr,"POS"],breaks=locs))
  i.snp=which(!duplicated(tmp))
  selected.snp=c(selected.snp,snp.map[snp.map$CHROM == chr,"ID"][i.snp[-1]])
#  print(c(df$n.causal[chr],length(i.snp)-1))
  re=rbind(re,data.frame(exp=df$n.causal[chr],obs=length(i.snp)-1))
}

# select n.causal of the selected SNP
i=sort(sample(length(selected.snp),n.causal))
selected.snp=selected.snp[i]

snp=snp.map[selected.snp,]

# determine the effect sizes for the two population
rownames(snp)=snp$ID
# randomly assing reference or alternative allele as the risk allele
snp$riskAllele=sample(c("REF","ALT"),nrow(snp),replace=T)
# let the frequency reflect the frequency of the risk allele
snp$CEU.frq[snp$riskAllele == "REF"]=1-snp$CEU.frq[snp$riskAllele == "REF"]
snp$YRI.frq[snp$riskAllele == "REF"]=1-snp$YRI.frq[snp$riskAllele == "REF"]

# equal weights
snp$CEU.wght=1
snp$YRI.wght=1

### assign the weights based on beta distribution
#snp$CEU.wght=rbeta(nrow(snp),1,10)
#snp$YRI.wght=rbeta(nrow(snp),1,10)

# turn off a random 25% of the SNP
#snp$CEU.wght[sample(n.causal,1000)]=0
#snp$YRI.wght[sample(n.causal,1000)]=0

# redistribute the weights
snp$CEU.wght=snp$CEU.wght/sum(snp$CEU.wght)
snp$YRI.wght=snp$YRI.wght/sum(snp$YRI.wght)



# estimate the effect sizes from the genetic variance
snp$CEU.beta=sqrt((params$run.info$varg$CEU*snp$CEU.wght)/(2*snp$CEU.frq*(1-snp$CEU.frq)))
snp$YRI.beta=sqrt((params$run.info$varg$YRI*snp$YRI.wght)/(2*snp$YRI.frq*(1-snp$YRI.frq)))


# randomly assign positive and negative betas
snp$CEU.beta=sample(c(-1,1),nrow(snp),replace=T)*snp$CEU.beta
snp$YRI.beta=sign(snp$CEU.beta)*snp$YRI.beta

# make sure the variance is maintained
sum(2*snp$CEU.frq*snp$CEU.beta)
sum(2*snp$CEU.frq*(1-snp$CEU.frq)*snp$CEU.beta^2)
sum(2*snp$YRI.frq*snp$YRI.beta)
sum(2*snp$YRI.frq*(1-snp$YRI.frq)*snp$YRI.beta^2)

### determine the correlation
(cor.CY=cor(2*snp$CEU.frq*snp$CEU.beta,2*snp$YRI.frq*snp$YRI.beta))


# save the selected parameters
params$run.info$sel.snp=snp #####START SELECTING SNP
#########SELECTING SNP ENDS HERE####
# determine the mean genetic contribution
params$run.info$gnt.mn=data.frame(CEU=sum(2*snp$CEU.beta*snp$CEU.frq), YRI=sum(2*snp$YRI.beta*snp$YRI.frq))


# set up the different subsets of data to generate with their parameters

anc="CEU.GWAS" ####CEU.REFERENCE1
N=20000     # training
params[[anc]]$ancestry=unlist(strsplit(anc,"\\."))[1]
params[[anc]]$seed=sample(100000,1)
params[[anc]]$chr=1:22
params[[anc]]$blocks=unique(snp.map$CEU.blk)
params[[anc]]$n.haplo=358
params[[anc]]$n.causal=nrow(snp) #set to NA
params[[anc]]$prev=0.015
params[[anc]]$n.case=N/2 #need to be NA
params[[anc]]$n.control=N/2 #need to be NA
params[[anc]]$n.sample=min(1000000,N*50) # = N ####duplicate the above for several times

anc="CEU.TRN"
N=20000     # training
params[[anc]]$ancestry=unlist(strsplit(anc,"\\."))[1]
params[[anc]]$seed=sample(100000,1)
params[[anc]]$chr=1:22
params[[anc]]$blocks=unique(snp.map$CEU.blk)
params[[anc]]$n.haplo=358
params[[anc]]$n.causal=nrow(snp)
params[[anc]]$prev=0.015
params[[anc]]$n.case=N/2
params[[anc]]$n.control=N/2
params[[anc]]$n.sample=min(1000000,N*50)


anc="YRI.GWAS"
N=20000    # training, testing, validating 
params[[anc]]$ancestry=unlist(strsplit(anc,"\\."))[1]
params[[anc]]$seed=sample(100000,1)
params[[anc]]$chr=1:22
params[[anc]]$blocks=unique(snp.map$YRI.blk)
params[[anc]]$n.haplo=356
params[[anc]]$n.causal=nrow(snp)
params[[anc]]$prev=0.015
params[[anc]]$n.case=N/2
params[[anc]]$n.control=N/2
params[[anc]]$n.sample=min(1000000,N*50)

anc="YRI.TRN"
params[[anc]]$ancestry=unlist(strsplit(anc,"\\."))[1]
N=20000    # training, testing, validating 
params[[anc]]$seed=sample(100000,1)
params[[anc]]$chr=1:22
params[[anc]]$blocks=unique(snp.map$YRI.blk)
params[[anc]]$n.haplo=356
params[[anc]]$n.causal=nrow(snp)
params[[anc]]$prev=0.015
params[[anc]]$n.case=N/2
params[[anc]]$n.control=N/2
params[[anc]]$n.sample=min(1000000,N*50)



# save the simulation parameters for this simulation
save(params,file=paste0(work.dir,"simulation-params.RData"))

# retain information that you want to carry over to another part of the simulation
rm.list=ls()
rm.list=rm.list[!(rm.list %in% c("i.sim","sims"))]
rm(list = rm.list); flush.console()

