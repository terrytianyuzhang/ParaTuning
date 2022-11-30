#rm(list=ls()); gc()
options(stringsAsFactors = F)

library(data.table)
library(doParallel)
library(snpStats)
library(R.utils)

# this version uses the heritability and the underlying normal to simulate the outcome

# #### software needed
# plink="/data3/Software/Plink/plink"
# plink2="/data3/Software/Plink2/plink2"
# directory changes
if(!exists("plink") | !exists("plink2")){
  plink <- "/usr/local/bin/plink"
  plink2 <- "/usr/local/bin/plink2"
}

# Bert has 32 cores
# nodes <- 32
# Try less on Kathryn's computer
nodes <- 24 

#### load the functions that are needed
source("simulation-functions.R")

#### load the parameters for the simulation
load(paste0("Work/Sim-",i.sim,"/simulation-params.RData"))
## to generate the reference population use
#load("Data/Reference-LDblocks/simulation-params.RData")


###############specify the population to use ##########################
for(set in names(params)[-1]){

  #######################################################################
  # generate causal genotypes for a large set of samples of cases and controls
  set.seed(params[[set]]$seed)
  selected.cc.haplotypes=NULL
  if(!is.na(params[[set]]$n.case) & !is.na(params[[set]]$n.control)){
    n.cycle=0
    n.case=0; n.control=0
    while(n.case < params[[set]]$n.case | n.control < params[[set]]$n.control){
      print(c(n.case,n.control)); flush.console()
      cc.haplotypes=data.frame(IID=paste(params[[set]]$ancestry,as.integer((1:params[[set]]$n.sample)+n.cycle*params[[set]]$n.sample),sep="."),DX=0,GNT=0,UNL=0)
      rownames(cc.haplotypes)=cc.haplotypes$IID
      # for each of the chromosomes randomly selected the two haplotypes

      system.time(tmp<-mclapply(params[[set]]$blocks,selectHaplotypes,n.sample=params[[set]]$n.sample,n.haplo=params[[set]]$n.haplo,mc.cores=50,mc.preschedule = T,mc.silent=T))
      tmp=cbind.data.frame(tmp)
      
      ###
      print("finished the 1st mclapply")
      
      ####
      cc.haplotypes=cbind.data.frame(cc.haplotypes,tmp)
      rm(tmp); gc()
      
    
      #### determine the binary outcome
      # for each block, calculate the genetic values
      blocks=unique(params$run.info$sel.snp[,paste0(params[[set]]$ancestry,".blk")])
      block.12=which(unlist(lapply(strsplit(gsub("H:","",colnames(cc.haplotypes)),"\\."),`[[`,1)) %in% blocks)
#      rm(gnt.contrib)
      # all in memory, tends to be slow.
#      system.time(gnt.contrib<-mclapply(blocks,determineGeneticValues.ld.blocks,anc=params[[set]]$ancestry,sel.snp=params$run.info$sel.snp,cc.haplotypes=cc.haplotypes[,block.12],
#                                     haplo.dir=paste0(params$run.info$main.dir,"Tmp/",params[[set]]$ancestry,"_haplotypes/"),
#                                     mc.cores=50,mc.preschedule = T,mc.silent=T))
#      gnt.contrib=rbindlist(gnt.contrib)
#      # sum across the chromosomes
#      gnt.contrib=colSums(gnt.contrib)
      
      # alternative using a write and read, much faster
      tmp.dir=paste0(params$run.info$work.dir,set,"/GNT.contrib")
      dir.create(tmp.dir,recursive = T,showWarnings = F)
      
      print("haplo.dir ir")
      print(paste0(params$run.info$main.dir,"Tmp/",params[[set]]$ancestry,"_haplotypes/"))
      
      system.time(re.block<-mclapply(blocks,determineGeneticValues.ld.blocks.with.write.ref.alt,anc=params[[set]]$ancestry,sel.snp=params$run.info$sel.snp,cc.haplotypes=cc.haplotypes[,block.12],
                                     haplo.dir=paste0(params$run.info$main.dir,"Tmp/",params[[set]]$ancestry,"_haplotypes/"),tmp.dir=tmp.dir,
                                     mc.cores= nodes,mc.preschedule = T,mc.silent=T))
      
      ###
      print("finished the 2nd mclapply")
      
      # collect the genetic contributions
      gnt.contrib=rep(0,nrow(cc.haplotypes))
      for(block in blocks){
        gnt.contrib=gnt.contrib+fread(paste0(tmp.dir,"/gnt.contrib-",block,".txt"),header=F,data.table=F,nThread= nodes)
      }
      gnt.contrib=unlist(gnt.contrib)
      names(gnt.contrib)=cc.haplotypes$IID
      
      # remove the temporary files
      system.command=paste("rm",paste0(tmp.dir,"/gnt.contrib-","*.txt"))
      system(system.command)

      # adjust for the mean of the population
      gnt.contrib=gnt.contrib-as.numeric(params$run.info$gnt.mn[params[[set]]$ancestry])
      cc.haplotypes$GNT=gnt.contrib
      # add a random residual
      underlying.normal=gnt.contrib+rnorm(length(gnt.contrib),mean = 0,sd=sqrt(params$run.info$vare))
      cc.haplotypes$UNL=underlying.normal
      
      
#      # turn this into a case control 
      cc.haplotypes$DX=.bincode(underlying.normal,breaks = c(-Inf,qnorm(params[[set]]$prev,mean=0,sd=as.numeric(sqrt(params$run.info$varg[params[[set]]$ancestry]+params$run.info$vare)),lower.tail = F),Inf))
#      cc.haplotypes$DX=rbinom(params[[set]]$n.sample,1,boot::inv.logit(log(params[[set]]$prev)+underlying.normal))+1      
      
      
      #### select samples as cases and control
      # cases 
      selected.cc.haplotypes=rbind(selected.cc.haplotypes,cc.haplotypes[sample(rownames(cc.haplotypes)[cc.haplotypes$DX == 2],min(sum(cc.haplotypes$DX == 2),params[[set]]$n.case-n.case)),])  
      # controls
      selected.cc.haplotypes=rbind(selected.cc.haplotypes,cc.haplotypes[sample(rownames(cc.haplotypes)[cc.haplotypes$DX == 1],min(sum(cc.haplotypes$DX == 1),params[[set]]$n.control-n.control)),])  
    
      
      # count the cases and controls
      n.case=sum(selected.cc.haplotypes$DX == 2)
      n.control=sum(selected.cc.haplotypes$DX == 1)
      n.cycle=n.cycle+1
    }
    print(c(n.case,n.control)); flush.console()
    rm(cc.haplotypes,gnt.contrib); gc()
  }else{
    # create the reference population
    cc.haplotypes=data.frame(IID=paste(params[[set]]$ancestry,as.integer(1:params[[set]]$n.sample),"REF",sep="."),DX=1)
    rownames(cc.haplotypes)=cc.haplotypes$IID
    
    system.time(tmp<-mclapply(params[[set]]$blocks,selectHaplotypes,n.sample=params[[set]]$n.sample,n.haplo=params[[set]]$n.haplo,mc.cores=50,mc.preschedule = T,mc.silent=T))
    tmp=cbind.data.frame(tmp)
    
    ###
    print("finished the 3rd mclapply")
    
    ####
    cc.haplotypes=cbind.data.frame(cc.haplotypes,tmp)
    rm(tmp); gc()

    selected.cc.haplotypes=cc.haplotypes
    rm(cc.haplotypes,haplotypes);gc()
  }
  # save the selected haplotypes for this run
  set.dir=paste0(params$run.info$work.dir,set,"/")
  dir.create(path =set.dir,showWarnings = F)
  fwrite(selected.cc.haplotypes,paste0(params$run.info$work.dir,set,"/sample-haplotypes.txt"),row.names=F,col.names=T,quote=F,sep="\t")
      
  
  ############################################################################## 
  # create the the directory to use for this data
  
  # and create a path for temporarily storing intermediate files
  dir.create(path =paste0(set.dir,"Tmp/"),showWarnings = F)
  
  
  # Now for each of the ld block create the genotypes to use in plink using the samples and haplotypes that were selected
  
  blocks=gsub("H:","",unique(unlist(lapply(strsplit(colnames(selected.cc.haplotypes),"\\."),`[[`,1))))
  blocks=blocks[!(blocks %in% c("IID","DX","GNT","UNL"))]
  
  haplo.dir=paste0(params$run.info$main.dir,"Tmp/")
  gnt.tmp.dir=paste0(params$run.info$work.dir,set,"/GNT.block/")
  dir.create(gnt.tmp.dir,showWarnings = F, recursive = T)
  
  
  nodes= nodes # seems to work best when writing intermediate files
  system.time(x<-mclapply(blocks,multiNodeGenotypeSimulation.ld.blocks,set=set,anc=params[[set]]$ancestry,selected.cc.haplotypes = selected.cc.haplotypes,
                          main.dir=params$run.info$main.dir,haplo.dir=haplo.dir,gnt.dir=gnt.tmp.dir,mc.cores=nodes,mc.preschedule=T,mc.silent=F))
  
  ###
  print("finished the 4th mclapply")
  
  ## combine all the pieces into one plink data file
  df=data.frame(files=list.files(paste(gnt.tmp.dir,sep=""),pattern = ".bed",full.names = T))
  df$files=gsub(".bed","",df$files)
  # order the files by chromosome and location, this can be done based on the subset number
  df$chr=unlist(lapply(strsplit(df$file,"-chr"),`[[`,2))
  df$n=as.numeric(unlist(lapply(strsplit(df$chr,"-b"),`[[`,2)))
  df$chr=as.numeric(unlist(lapply(strsplit(df$chr,"-b"),`[[`,1)))
  df=df[order(df$chr,df$n),]
  
  write.table(df$files,paste(set.dir,"/bed.files",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
  
  plink2.command=paste(plink2,"--nonfounders","--allow-no-sex",
                       "--pmerge-list",paste(set.dir,"/bed.files",sep=""),"bfile",
                       "--out",paste(set.dir,"/",set,sep=""),
                       sep=" ")
  system(plink2.command)
  
  # remove the temporary files
  system.command=paste("rm",paste0(gnt.tmp.dir,"*chr*"))
  system(system.command)
}

# retain information that you want to carry over to another part of the simulation
rm.list=ls()
rm.list=rm.list[!(rm.list %in% c("i.sim","sims"))]
rm(list = rm.list); flush.console()

