selectHaplotypes <- function(block,n.sample,n.haplo){
  haplotypes=matrix(sample(n.haplo,n.sample*2,replace=T),nrow=n.sample,ncol=2)
  same=which(haplotypes[,1]== haplotypes[,2])
  while(length(same) > 0){
    haplotypes[same,2]=sample(params[[set]]$n.haplo,length(same),replace=T)
    same=which(haplotypes[,1]== haplotypes[,2])
  }
  # just give this data some row and column names for easier access
  colnames(haplotypes)=c(paste("H:",block,".1",sep=""),paste("H:",block,".2",sep=""))
  return(haplotypes)
}


determineGeneticValues <- function(chr,anc,sel.snp,cc.haplotypes,main.dir){
  if(sum(sel.snp$CHROM == chr) > 0){
    # read the haplotypes for this chromosome
    haplos=fread(paste(main.dir,"Data/",anc,"_haplotypes/",anc,"-chr",chr,"-qc.haplo.gz",sep=""),header=T,data.table=F)
    rownames(haplos)=haplos$ID; haplos=as.matrix(haplos[,-1])
    # select the markers that are in the causal snp
    if(sum(rownames(haplos) %in% sel.snp$ID) > 1){
      haplos=haplos[rownames(haplos) %in% sel.snp$ID,]
    }else{
      col.name=colnames(haplos)
      row.name=sel.snp$ID[sel.snp$CHROM == chr]
      haplos=matrix(haplos[rownames(haplos) %in% sel.snp$ID,],nrow=length(row.name),ncol=length(col.name),byrow=T)
      colnames(haplos)=col.name
      rownames(haplos)=row.name
    }
    # create the genetic value attributable to this chromosome
    gnt.contrib=as.vector(sel.snp[rownames(haplos),"beta"]%*%(haplos[,cc.haplotypes[,paste("H",chr,".1",sep="")]]+haplos[,cc.haplotypes[,paste("H",chr,".2",sep="")]]))
    names(gnt.contrib)=cc.haplotypes$IID
  }else{
    gnt.contrib=rep(0,nrow(cc.haplotypes))
    names(gnt.contrib)=cc.haplotypes$IID
  }
  return(data.frame(t(gnt.contrib)))
}

determineGeneticValues.H2 <- function(chr,anc,sel.snp,cc.haplotypes,main.dir){
  if(sum(sel.snp$CHROM == chr) > 0){
    # read the haplotypes for this chromosome
    haplos=fread(paste(main.dir,"Data/",anc,"_haplotypes/",anc,"-chr",chr,"-qc.haplo.gz",sep=""),header=T,data.table=F)
    rownames(haplos)=haplos$ID; haplos=as.matrix(haplos[,-1])
    # select the markers that are in the causal snp
    if(sum(rownames(haplos) %in% sel.snp$ID) > 1){
      haplos=haplos[rownames(haplos) %in% sel.snp$ID,]
    }else{
      col.name=colnames(haplos)
      row.name=sel.snp$ID[sel.snp$CHROM == chr]
      haplos=matrix(haplos[rownames(haplos) %in% sel.snp$ID,],nrow=length(row.name),ncol=length(col.name),byrow=T)
      colnames(haplos)=col.name
      rownames(haplos)=row.name
    }
    # create the genetic value attributable to this chromosome
    gnt.contrib=as.vector(sel.snp[rownames(haplos),paste0(anc,".beta")]%*%(haplos[,cc.haplotypes[,paste("H",chr,".1",sep="")]]+haplos[,cc.haplotypes[,paste("H",chr,".2",sep="")]]))
    names(gnt.contrib)=cc.haplotypes$IID
  }else{
    gnt.contrib=rep(0,nrow(cc.haplotypes))
    names(gnt.contrib)=cc.haplotypes$IID
  }
  return(data.frame(t(gnt.contrib)))
}

determineGeneticValues.ld.blocks <- function(block,anc,sel.snp,cc.haplotypes,haplo.dir){
  if(sum(sel.snp[paste0(anc,".blk")] == block) > 0){
    # read the haplotypes for this chromosome
    haplos=fread(paste(haplo.dir,anc,"-",block,"-qc.haplo.gz",sep=""),header=T,data.table=F)
    rownames(haplos)=haplos$ID; haplos=as.matrix(haplos[,-1])
    # select the markers that are in the causal snp
    if(sum(rownames(haplos) %in% sel.snp$ID) > 1){
      haplos=haplos[rownames(haplos) %in% sel.snp$ID,]
    }else{
      col.name=colnames(haplos)
      row.name=sel.snp$ID[sel.snp[paste0(anc,".blk")] == block]
      haplos=matrix(haplos[rownames(haplos) %in% sel.snp$ID,],nrow=length(row.name),ncol=length(col.name),byrow=T)
      colnames(haplos)=col.name
      rownames(haplos)=row.name
    }
    # create the genetic value attributable to this chromosome
    gnt.contrib=as.vector(sel.snp[rownames(haplos),paste0(anc,".beta")]%*%(haplos[,cc.haplotypes[,paste("H:",block,".1",sep="")]]+haplos[,cc.haplotypes[,paste("H:",block,".2",sep="")]]))
    names(gnt.contrib)=cc.haplotypes$IID
  }else{
    gnt.contrib=rep(0,nrow(cc.haplotypes))
    names(gnt.contrib)=cc.haplotypes$IID
  }
  return(data.frame(t(gnt.contrib)))
}


determineGeneticValues.ld.blocks.with.write <- function(block,anc,sel.snp,cc.haplotypes,haplo.dir,tmp.dir){
  if(sum(sel.snp[paste0(anc,".blk")] == block) > 0){
    # read the haplotypes for this chromosome
    haplos=fread(paste(haplo.dir,anc,"-",block,"-qc.haplo.gz",sep=""),header=T,data.table=F)
    rownames(haplos)=haplos$ID; haplos=as.matrix(haplos[,-1])
    # select the markers that are in the causal snp
    if(sum(rownames(haplos) %in% sel.snp$ID) > 1){
      haplos=haplos[rownames(haplos) %in% sel.snp$ID,]
    }else{
      col.name=colnames(haplos)
      row.name=sel.snp$ID[sel.snp[paste0(anc,".blk")] == block]
      haplos=matrix(haplos[rownames(haplos) %in% sel.snp$ID,],nrow=length(row.name),ncol=length(col.name),byrow=T)
      colnames(haplos)=col.name
      rownames(haplos)=row.name
    }
    # create the genetic value attributable to this chromosome
    gnt.contrib=as.vector(sel.snp[rownames(haplos),paste0(anc,".beta")]%*%(haplos[,cc.haplotypes[,paste("H:",block,".1",sep="")]]+haplos[,cc.haplotypes[,paste("H:",block,".2",sep="")]]))
    names(gnt.contrib)=cc.haplotypes$IID
  }else{
    gnt.contrib=rep(0,nrow(cc.haplotypes))
    names(gnt.contrib)=cc.haplotypes$IID
  }
  fwrite(as.data.frame(gnt.contrib),paste0(tmp.dir,"/gnt.contrib-",block,".txt"),row.names=F,col.names=F,quote=F,sep="\t",na=NA)
  return(block)
}


determineGeneticValues.ld.blocks.with.write.ref.alt <- function(block,anc,sel.snp,cc.haplotypes,haplo.dir,tmp.dir){
  if(sum(sel.snp[paste0(anc,".blk")] == block) > 0){
    # read the haplotypes for this chromosome
    haplos=fread(paste(haplo.dir,anc,"-",block,"-qc.haplo.gz",sep=""),header=T,data.table=F)
    rownames(haplos)=haplos$ID; haplos=as.matrix(haplos[,-1])
    # select the markers that are in the causal snp
    if(sum(rownames(haplos) %in% sel.snp$ID) > 1){
      haplos=haplos[rownames(haplos) %in% sel.snp$ID,]
    }else{
      col.name=colnames(haplos)
      row.name=sel.snp$ID[sel.snp[paste0(anc,".blk")] == block]
      haplos=matrix(haplos[rownames(haplos) %in% sel.snp$ID,],nrow=length(row.name),ncol=length(col.name),byrow=T)
      colnames(haplos)=col.name
      rownames(haplos)=row.name
    }
    # check to see if the risk allele count needs to be reversed
    i.ref=which(sel.snp[rownames(haplos),"riskAllele"]=="REF")
    haplos[i.ref,]=1-haplos[i.ref,]
    # create the genetic value attributable to this chromosome
    gnt.contrib=as.vector(sel.snp[rownames(haplos),paste0(anc,".beta")]%*%(haplos[,cc.haplotypes[,paste("H:",block,".1",sep="")]]+haplos[,cc.haplotypes[,paste("H:",block,".2",sep="")]]))
    names(gnt.contrib)=cc.haplotypes$IID
  }else{
    gnt.contrib=rep(0,nrow(cc.haplotypes))
    names(gnt.contrib)=cc.haplotypes$IID
  }
  fwrite(as.data.frame(gnt.contrib),paste0(tmp.dir,"/gnt.contrib-",block,".txt"),row.names=F,col.names=F,quote=F,sep="\t",na=NA)
  return(block)
}



###### CREATE THE COMPLETE SET OF GENOTYPES USING MULTIPLE NODES
multiNodeGenotypeSimulation <- function(i.subset,subsets,anc,selected.cc.haplotypes,main.dir,set.dir){
  #### read the haplotypes
  chr=subsets[i.subset,"chr"]
  haplos=fread(paste(main.dir,"Data/",anc,"_haplotypes/",anc,"-chr",chr,"-qc.haplo.gz",sep=""),header=T,data.table=F)
  rownames(haplos)=haplos$ID
  # select the subset of haplotypes needed
  haplos=as.matrix(t(haplos[subsets[i.subset,"start"]:subsets[i.subset,"end"],-1]))  # transpose the haplo matrix for quicker processing
  # We count the number of ALTernative alleles, due to a quirk in write.plink in the writing step it counts the number of REFerence alleles
  # there we need to recode the haplotypes as 1-haplotypes. We can then be sure that once the results come out we still have the alternative allele
  # as the number of counted alleles (it is confusing, but trist me.)
  haplos=1-haplos
  # create the genotypes for this subset
  system.time(gnt<-haplos[selected.cc.haplotypes[,paste("H",chr,".1",sep="")],]+
                haplos[selected.cc.haplotypes[,paste("H",chr,".2",sep="")],])
  rownames(gnt)=selected.cc.haplotypes$IID
  rm(haplos)
  
  
  # create a map data.frame
  map=data.frame(CHR=chr,SNP=colnames(gnt),CM=0,
                 POS=as.numeric(unlist(lapply(strsplit(colnames(gnt),":"),`[[`,2))),
                 A1=unlist(lapply(strsplit(colnames(gnt),":"),`[[`,4)),
                 A2=unlist(lapply(strsplit(colnames(gnt),":"),`[[`,3)))
  map$POS[grep("000000:",map$SNP)]=map$POS[grep("000000:",map$SNP)]+1 # this is a quirk in write.plink where any location that ends in 000000: is put in scientific notation
  
  # put genotypes in SnpMatrix format, this is the time  intensive part
  system.time(GNT<-as(gnt,"SnpMatrix"))
  rm(gnt)
  
  # write the genotypes to a plink binary format file
  system.time(write.plink(paste(set.dir,"Tmp/",anc,"-chr",chr,"-s",i.subset,sep=""),
                          snp.major=T,
                          snps=GNT,
                          pedigree=selected.cc.haplotypes$IID,
                          id=selected.cc.haplotypes$IID,
                          phenotype=selected.cc.haplotypes$DX,
                          chromosome=map$CHR,
                          genetic.distance=map$CM,
                          position=map$POS,
                          allele.1=map$A1,
                          allele.2=map$A2,
                          na.code=0))
  
  return(i.subset)
  
}


multiNodeGenotypeSimulation.ld.blocks <- function(block,set,anc,selected.cc.haplotypes,main.dir,haplo.dir,gnt.dir){
  #### read the haplotypes
  haplos=fread(paste(haplo.dir,anc,"_haplotypes/",anc,"-",block,"-qc.haplo.gz",sep=""),header=T,data.table=F)
  rownames(haplos)=haplos$ID
  # select the subset of haplotypes needed
  haplos=as.matrix(t(haplos[,-1]))  # transpose the haplo matrix for quicker processing
  # We count the number of ALTernative alleles, due to a quirk in write.plink in the writing step it counts the number of REFerence alleles
  # there we need to recode the haplotypes as 1-haplotypes. We can then be sure that once the results come out we still have the alternative allele
  # as the number of counted alleles (it is confusing, but trust me.)
  haplos=1-haplos
  # create the genotypes for this subset
  system.time(gnt<-haplos[selected.cc.haplotypes[,paste("H:",block,".1",sep="")],]+
                haplos[selected.cc.haplotypes[,paste("H:",block,".2",sep="")],])
  rownames(gnt)=selected.cc.haplotypes$IID
  rm(haplos)
  
  
  # create a map data.frame
  map=data.frame(CHR=as.numeric(gsub("chr","",unlist(strsplit(block,"-"))[1])),SNP=colnames(gnt),CM=0,
                 POS=as.numeric(unlist(lapply(strsplit(colnames(gnt),":"),`[[`,2))),
                 A1=unlist(lapply(strsplit(colnames(gnt),":"),`[[`,4)),
                 A2=unlist(lapply(strsplit(colnames(gnt),":"),`[[`,3)))
  map$POS[grep("000000:",map$SNP)]=map$POS[grep("000000:",map$SNP)]+1 # this is a quirk in write.plink where any location that ends in 000000: is put in scientific notation
  
  # put genotypes in SnpMatrix format, this is the time  intensive part
  system.time(GNT<-as(gnt,"SnpMatrix"))
  rm(gnt)
  
  # write the genotypes to a plink binary format file
  system.time(write.plink(paste(gnt.dir,anc,"-",block,sep=""),
                          snp.major=T,
                          snps=GNT,
                          pedigree=selected.cc.haplotypes$IID,
                          id=selected.cc.haplotypes$IID,
                          phenotype=selected.cc.haplotypes$DX,
                          chromosome=map$CHR,
                          genetic.distance=map$CM,
                          position=map$POS,
                          allele.1=map$A1,
                          allele.2=map$A2,
                          na.code=0))
  
  return(block)
  
}

### clump the markers by chromsomes
clumpingFunction<-function(chr,r2,set,set.dir,plink=plink){
  # create the directory
  dir.create(paste0(set.dir,"Clumping/"),showWarnings = F)
  # run clumping
  plink.command=paste(plink,"--nonfounders","--allow-no-sex","--threads",32,
                      "--bfile",paste(set.dir,set,sep=""),
                      "--chr",chr,
                      "--clump-p1",1,
                      "--clump-p2",1,
                      "--clump-r2",r2,
                      "--clump-kb",500,
                      "--clump",paste(set.dir,"Assoc/",set,"-chr",chr,".PHENO1.glm.logistic.hybrid",sep=""),
                      "--clump-snp-field","ID",
                      "--clump-field","P",
                      "--out",paste0(set.dir,"Clumping/",set,"-chr",chr,sprintf("-r2_%2.1f",r2)),
                      sep=" ")
  
  system(plink.command)
  return(chr)
  
}
