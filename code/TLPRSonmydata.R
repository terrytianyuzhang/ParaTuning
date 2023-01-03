##plink-1.9 needs to be pre-installed.
# library(devtools)
# install_github("terrytianyuzhang/TLPRS")
#########Step 1. Run TL-PRS using example data#########
library(data.table)
library(lassosum)
library(TLPRS)
library(parallel)
library(pROC)

if(!exists("plink") | !exists("plink2")){
  plink <- "/usr/local/bin/plink"
  plink2 <- "/usr/local/bin/plink2"
}

if(!exists("i.sim")){
  i.sim <- 800
}

# load the parameters for this simulation
load(paste0("Work/Sim-",i.sim,"/simulation-params.RData"))
main.dir=params$run.info$main.dir #"/raid6/Tianyu/PRS/SimulationPipeline/"
work.dir=params$run.info$work.dir #"/raid6/Tianyu/PRS/SimulationPipeline/Work/Sim-800/"


ped_file <- paste0(work.dir, '/YRI.TUNE/YRI-TUNE-fam-forTL.txt')
Y_name <- "Y" ##what is Y name
Ytype <- "B" 

train_file <- paste0(work.dir, 'YRI.TRN/YRI.TRN')

validate_file <- paste0(work.dir, 'YRI.TUNE/YRI.TUNE')

sum_stats_file <- paste0(work.dir, 'CEU.TRN/CEU-TRN-Lassosum-forTL.txt')

target_sumstats_file <- paste0(work.dir, 'YRI.TRN/YRI-TRN-GWAS-forTL.txt')

LDblocks <- "AFR.hg38" ##it used to be "AFR.hg19"

outfile <- paste0(main.dir, 'firstTL')

###inside TL_PRS
test_file = validate_file
tempfile=outfile
Covar_name = NULL
out1=TLPRS:::PRStr_main_check(ped_file,Covar_name,Y_name, Ytype,train_file,test_file,sum_stats_file,LDblocks)
if (out1!=0){stop(out1)}

sum_stats=data.frame(fread(sum_stats_file)) #right now the beta is correct
if (ncol(sum_stats)==3){ 
  if (sum(colnames(sum_stats) %in% c("V1","V2","V3"))==3){
    colnames(sum_stats)=c("SNP","A1","Beta")
  }
} 

sum_stats=sum_stats[,c("SNP","A1","Beta")] 
sum_stats_file=paste0(tempfile,"_original_sum_stats.txt")
write.table(sum_stats, file=sum_stats_file,col.names=F,row.names=F,quote=F)

ped=data.frame(fread(ped_file,header=T))[,setdiff(c("FID","IID",Covar_name,Y_name),"")]

##obj=calculate_betaPRS(train_file,sum_stats_file,ped,Covar_name,Y_name,paste0(tempfile,"_step0") ) ##need to remove sum_stats_file and plink command later.

sum_stats_target=fread(target_sumstats_file)
sum_stats_target=merge(sum_stats,sum_stats_target,by="SNP",sort=F)
if (sum(sum_stats_target$p<=1E-320)>0){ sum_stats_target$p[sum_stats_target$p<=1E-320]=1E-320}

sum_stats_target$cor=lassosum::p2cor(p = sum_stats_target$p, 
                                     n = median(sum_stats_target$N,na.rm=T), 
                                     sign=sum_stats_target$beta)
# > head(sum_stats_target)
# SNP A1.x Beta A1.y         beta    N        p           cor
# 1 1:1962845:T:C    C    0    C  0.005862780 4000 0.913745  0.0017131611
# 2 1:1962899:A:C    C    0    C -0.003230212 4000 0.967944 -0.0006356117
# 3 1:1963406:G:A    A    0    A -0.019352049 4000 0.780846 -0.0044003466
# 4 1:1963538:T:C    C    0    C -0.012912002 4000 0.772044 -0.0045820035
# 5 1:1963738:C:T    T    0    T -0.019352049 4000 0.780846 -0.0044003466
# 6 1:1964101:A:G    G    0    G -0.019352049 4000 0.780846 -0.0044003466

flag=which(sum_stats_target$A1.x !=sum_stats_target$A1.y)
# > flag
# integer(0) #no miss-match
if (length(flag)>0){sum_stats_target$cor[flag]=-sum_stats_target$cor[flag]}
sum_stats_target=sum_stats_target[,c("SNP","A1.x","Beta","cor")];colnames(sum_stats_target)[2]="A1";
# > head(sum_stats_target)
# SNP A1 Beta           cor
# 1 1:1962845:T:C  C    0  0.0017131611
# 2 1:1962899:A:C  C    0 -0.0006356117
# 3 1:1963406:G:A  A    0 -0.0044003466
# 4 1:1963538:T:C  C    0 -0.0045820035
# 5 1:1963738:C:T  T    0 -0.0044003466
# 6 1:1964101:A:G  G    0 -0.0044003466
gc()

# next step is time-consuming
# beta_list=as.data.frame(PRStr_calculation2(sum_stats_target, train_file, sum_stats, LDblocks, cluster=cluster,temp.file=paste0(tempfile,"_step1"),
#                                            plink = plink))
####inside PRStr_calculation2
if(1){
    temp.file=paste0(tempfile,"_step1")
    possible.LDblocks <- c("EUR.hg19", "AFR.hg19", "ASN.hg19", 
                           "EUR.hg38", "AFR.hg38", "ASN.hg38") 
    if(!is.null(LDblocks)) { #T
      if(is.character(LDblocks) && length(LDblocks) == 1) { #T
        if(LDblocks %in% possible.LDblocks) { #T
          LDblocks <- data.table::fread(system.file(paste0("data/Berisa.",  LDblocks, ".bed"),  package="lassosum"), header=T)
          # > head(LDblocks)
          # chr   start    stop
          # 1: chr1 1641723 2433551
          # 2: chr1 2433551 4224127
          # 3: chr1 4224127 4794254
          # 4: chr1 4794254 5854827
          # 5: chr1 5854827 7110659
          # 6: chr1 7110659 8142348
        } else {
          stop(paste("I cannot recognize this LDblock. Specify one of", 
                     paste(possible.LDblocks, collapse=", ")))
        }
      }
      if(is.factor(LDblocks)) LDblocks <- as.integer(LDblocks)
      if(is.vector(LDblocks)) stopifnot(length(LDblocks) == length(cor)) else 
        if(is.data.frame(LDblocks) || is.data.table(LDblocks)) {
          LDblocks <- as.data.frame(LDblocks)
          stopifnot(ncol(LDblocks) == 3)
          stopifnot(all(LDblocks[,3] >= LDblocks[,2]))
          LDblocks[,1] <- as.character(sub("^chr", "", LDblocks[,1], ignore.case = T))
        }
    } else {
      stop(paste0("LDblocks must be specified. Specify one of ", 
                  paste(possible.LDblocks, collapse=", "), 
                  ". Alternatively, give an integer vector defining the blocks, ", 
                  "or a .bed file with three columns read as a data.frame."))
    }
    
    ref.bim <- fread(paste0(train_file, ".bim"))
    ref.bim$V1 <- as.character(sub("^chr", "", ref.bim$V1, ignore.case = T))
    ref.bim$order=1:nrow(ref.bim)
    bim_sum_stats=merge(ref.bim, sum_stats_target,by.x="V2",by.y="SNP",order=F)
    bim_sum_stats=bim_sum_stats[order(bim_sum_stats$order),]
    bim_sum_stats$Beta2=NA
    
    # > head(bim_sum_stats)
    # V2 V1 V3      V4 V5 V6 order A1 Beta           cor Beta2
    # 1: 1:1962845:T:C  1  0 1962845  C  T     1  C    0  0.0017131611    NA
    # 2: 1:1962899:A:C  1  0 1962899  C  A     2  C    0 -0.0006356117    NA
    # 3: 1:1963406:G:A  1  0 1963406  A  G     3  A    0 -0.0044003466    NA
    # 4: 1:1963538:T:C  1  0 1963538  C  T     4  C    0 -0.0045820035    NA
    # 5: 1:1963738:C:T  1  0 1963738  T  C     5  T    0 -0.0044003466    NA
    # 6: 1:1964101:A:G  1  0 1964101  G  A     6  G    0 -0.0044003466    NA
    
    flag1=which(bim_sum_stats$V5==bim_sum_stats$A1)
    if (length(flag1)>0){  bim_sum_stats$Beta2[flag1]=bim_sum_stats$Beta[flag1]}
    flag2=which(bim_sum_stats$V6==bim_sum_stats$A1)
    # > flag2
    # integer(0)
    if (length(flag2)>0){  bim_sum_stats$Beta2[flag2]=-bim_sum_stats$Beta[flag2];  bim_sum_stats$cor[flag2]=-bim_sum_stats$cor[flag2];}
    
    bim_sum_stats=bim_sum_stats[which(! is.na(bim_sum_stats$Beta2)),c("V2","V1","V4","V5","V6","order","Beta2","cor")]
    ##so fat the beta's are correct
    
    ref.extract <- rep(FALSE, nrow(ref.bim))
    ref.extract[bim_sum_stats$order] <- TRUE
    # > summary(ref.extract)
    # Mode   FALSE    TRUE 
    # logical  896272 4734473 
    
    if(!is.null(LDblocks)) {
      LDblocks2 <- TLPRS:::splitgenome2(CHR = ref.bim$V1[ ref.extract], 
                                POS = ref.bim$V4[ ref.extract],
                                ref.CHR = LDblocks[,1], 
                                ref.breaks = LDblocks[,3])
      # Assumes base 1 for the 3rd column of LDblocks (like normal bed files)
    } 
    
    # > head(LDblocks2[[1]])
    # [1] 1_[1,2.4336e+06] 1_[1,2.4336e+06] 1_[1,2.4336e+06] 1_[1,2.4336e+06]
    # [5] 1_[1,2.4336e+06] 1_[1,2.4336e+06]
    # 2583 Levels: 1_[1,2.4336e+06] ... 22_(5.08e+07,Inf]
    # > head(LDblocks2[[2]])
    # chr   start     end counts
    # 1   1       1 2433551    990
    # 2   1 2433551 4224127   3660
    # 3   1 4224127 4794254   1672
    # 4   1 4794254 5854827   2509
    # 5   1 5854827 7110659   2149
    # 6   1 7110659 8142348   1802
    i <- unique(LDblocks2[[1]])[[1]]
    
    # results.list <- lapply(unique(LDblocks2[[1]]), function(i) {
    #   block_calculation2(cor=bim_sum_stats[which(LDblocks2[[1]]==i),], 
    #                      num=which(i==unique(LDblocks2[[1]])),
    #                      train_file=train_file,
    #                      nsnp=nrow(bim_sum_stats),
    #                      temp.file,
    #                      plink)
    # })
    if(1){
        #####NOW we go into block_calculation2
        cor=bim_sum_stats[which(LDblocks2[[1]]==i),]
        # > cor
        # V2 V1      V4 V5 V6 order Beta2           cor
        # 1: 1:1962845:T:C  1 1962845  C  T     1     0  0.0017131611
        # 2: 1:1962899:A:C  1 1962899  C  A     2     0 -0.0006356117
        # 3: 1:1963406:G:A  1 1963406  A  G     3     0 -0.0044003466
        # 4: 1:1963538:T:C  1 1963538  C  T     4     0 -0.0045820035
        # 5: 1:1963738:C:T  1 1963738  T  C     5     0 -0.0044003466
        # ---                                                         
        #   986: 1:2430393:C:T  1 2430393  T  C  1122     0  0.0022518035
        # 987: 1:2430581:G:T  1 2430581  T  G  1123     0 -0.0281986335
        # 988: 1:2431388:C:T  1 2431388  T  C  1125     0 -0.0243703083
        # 989: 1:2431982:C:T  1 2431982  T  C  1127     0  0.0070223993
        # 990: 1:2433069:A:G  1 2433069  G  A  1128     0  0.0064016868
        
        # > cor[V2 == '1:1991822:G:A',]
        # V2 V1      V4 V5 V6 order      Beta2         cor
        # 1: 1:1991822:G:A  1 1991822  A  G   112 0.00641277 -0.01995747
        num=which(i==unique(LDblocks2[[1]])) ##num = 1
        nsnp=nrow(bim_sum_stats) #4734473
        
        temp_file=paste0(temp.file,"_block_",num)
        write.table(cor$V2,file=temp_file,col.names=F,row.names=F,quote=F)
        ##wrote the SNP names
        cmd = paste0(plink, " --bfile ",train_file," --extract ",temp_file,   " --recodeA  --out ", temp_file,"_Geno.txt")
        system(cmd)
        
        Gtemp=try(as.data.frame(fread(paste0(temp_file,"_Geno.txt.raw"),header=T)),silent=T)
        # > Gtemp[1:5,1:7]
        # FID     IID PAT MAT SEX PHENOTYPE 1:1962845:T:C_C
        # 1  YRI.70  YRI.70   0   0   0         1               0
        # 2  YRI.91  YRI.91   0   0   0         1               0
        # 3 YRI.165 YRI.165   0   0   0         2               1
        # 4 YRI.201 YRI.201   0   0   0         1               0
        # 5 YRI.240 YRI.240   0   0   0         2               1
        if (file.exists(temp_file)) {file.remove(temp_file)}
        if (file.exists(paste0(temp_file,"_Geno.txt.nosex"))) {file.remove(paste0(temp_file,"_Geno.txt.nosex"))}
        if (file.exists(paste0(temp_file,"_Geno.txt.log"))) {file.remove(paste0(temp_file,"_Geno.txt.log"))}
        if (file.exists(paste0(temp_file,"_Geno.txt.raw"))) {file.remove(paste0(temp_file,"_Geno.txt.raw"))}
        ##TTTT
        
        if (class(Gtemp)=="try-error"){
          return(NULL)
          #GG=diag(nrow(cor));colnames(GG)=paste0(cor$V2,"_",cor$V5) 
          #geno_info=as.data.frame(t(sapply(colnames(GG),   split_SNPandA1   ) )) ;colnames(geno_info)=c("SNP","A1")
          #geno_info$mean=NA; geno_info$maf=NA; geno_info$sd=NA
        }else{
          GG=cor(as.matrix(Gtemp[,7:ncol(Gtemp)]))
          geno_info=as.data.frame(t(sapply(colnames(Gtemp)[7:ncol(Gtemp)],   TLPRS:::split_SNPandA1   ) )) ;colnames(geno_info)=c("SNP","A1")
          geno_info$mean=colMeans(as.matrix(Gtemp[,7:ncol(Gtemp)]),na.rm=T); geno_info$maf=geno_info$mean/2; geno_info$sd=sqrt(2*geno_info$maf*(1-geno_info$maf))
          # > head(geno_info)
          # SNP A1    mean      maf        sd
          # 1:1962845:T:C_C 1:1962845:T:C  C 0.43650 0.218250 0.5841523
          # 1:1962899:A:C_C 1:1962899:A:C  C 0.17175 0.085875 0.3962335
          # 1:1963406:G:A_A 1:1963406:G:A  A 0.22650 0.113250 0.4481617
          # 1:1963538:T:C_C 1:1963538:T:C  C 0.86825 0.434125 0.7009429
          # 1:1963738:C:T_T 1:1963738:C:T  T 0.22650 0.113250 0.4481617
          # 1:1964101:A:G_G 1:1964101:A:G  G 0.22650 0.113250 0.4481617
        
          # > geno_info[geno_info$SNP == '1:1991822:G:A',]
          # SNP A1    mean      maf        sd
          # 1:1991822:G:A_G 1:1991822:G:A  G 0.12425 0.062125 0.3413663
        }
        
        list1=which(geno_info$sd==0)
        if (length(list1)>0){ #length(list1) = 0
          geno_info=geno_info[-list1,]
          GG=GG[-list1,-list1]
        }
        if (nrow(geno_info)==0){
          return(NULL)
        } else {
          geno_info$order=1:nrow(geno_info)
          geno_info2=merge(cor[,c("V2","V5","Beta2","cor")],geno_info, by.x="V2",by.y="SNP",sort=F)
          flag_nomatch=which(geno_info2$A1 != geno_info2$V5)
          if (length(flag_nomatch)>0){
            geno_info2$Beta2[flag_nomatch]=-geno_info2$Beta2[flag_nomatch]
            geno_info2$cor[flag_nomatch]=-geno_info2$cor[flag_nomatch]
          }
          
          # > geno_info2[V2 == '1:1991822:G:A',]
          # V2 V5       Beta2        cor A1    mean      maf        sd order
          # 1: 1:1991822:G:A  A -0.00641277 0.01995747  G 0.12425 0.062125 0.3413663    98
          # 
          GG2=as.matrix(GG[geno_info2$order,geno_info2$order])
          # > GG2[1:3,1:3]
          # 1:1962845:T:C_C 1:1962899:A:C_C 1:1963406:G:A_A
          # 1:1962845:T:C_C       1.0000000      -0.1609391       0.6879457
          # 1:1962899:A:C_C      -0.1609391       1.0000000      -0.1028879
          # 1:1963406:G:A_A       0.6879457      -0.1028879       1.0000000
          gy=geno_info2$cor
          betatemp=geno_info2$Beta2*geno_info2$sd
          # > betatemp[90:100]
          # [1]  0.000000000  0.000000000  0.000000000  0.000000000  0.000000000
          # [6]  0.000000000  0.000000000  0.000000000 -0.002189104  0.000000000
          # [11] -0.005240255
          
          u0=gy-GG2%*%betatemp #this is gradient
          
          # > u0[90:100]
          # [1] -0.0123208635 -0.0032165046  0.0003374752 -0.0160604999  0.0065016681
          # [6]  0.0016290792 -0.0223940175 -0.0111549832  0.0240602572 -0.0143946666
          # [11]  0.0094579082
          beta.all=cbind(u0, betatemp)
          for (factor1 in c(1,10,100,1000)){
            k=1
            betatemp=beta.all[,2]
            u0=beta.all[,1]
            while (k<=15){
              ##betanew=c()
              learningrate=1/nsnp*factor1
              if (learningrate>1){learningrate=1}
              ##print(learningrate)
              for (j in 1:length(betatemp)){
                beta_old=betatemp[j]
                betatemp[j]=(learningrate*u0[j]+beta_old)/ 1
                u0=u0-GG2[,j]*(betatemp[j]-beta_old)
              }
              beta.all=cbind(beta.all,betatemp)
              k=k+1
            } 
          }
          geno_info2=cbind(geno_info2,beta.all)
          
          # > geno_info2[V2 == '1:1991822:G:A',1:15]
          # V2 V5       Beta2        cor A1    mean      maf        sd order
          # 1: 1:1991822:G:A  A -0.00641277 0.01995747  G 0.12425 0.062125 0.3413663    98
          # V1     betatemp     betatemp     betatemp     betatemp     betatemp
          # 1: 0.02406026 -0.002189104 -0.002189099 -0.002189093 -0.002189088 -0.002189083
          # return(geno_info2)
        } #else
    }
    results.list<-do.call("rbind", results.list)
}
# cl <- makeCluster(24, type="FORK")

##########FINALLY
#######CAME BACK TO TLPRS
beta_list <- geno_info2
# beta_list=as.data.frame(beta_list[,-c(5,9)]) ###you got rid of the correct A1?!?!
# ##the fifth column is the A1 and the nineth is the order.

beta_list[,2] <- beta_list[,5] ###now the second column is the correct A1
beta_list <- beta_list[,-c(5,9)]
beta_list <- as.data.frame(beta_list)


# beta_list <- fread(file=paste0(tempfile,"_beta.candidates.txt"))
# beta_list=as.data.frame(beta_list)

# out1=PRStr_tuning(beta_list, ped,Covar_name, Y_name, Ytype,test_file)
####into PRStr_tuning
if(1){#PRStr_tuning
  Beta.all = beta_list
  
  beta.info=Beta.all[,1:2]
  for (i in 9:ncol(Beta.all)){
    sdtemp=sd(Beta.all[,i],na.rm=T)
    if (sdtemp>1){
      Beta.all[,i:ncol(Beta.all)]=1;##print("fdsf")
    }
  }
  
  beta.all=Beta.all[,9:ncol(Beta.all)]/Beta.all$sd
  
  # PRS.all=Calculate_PRS(test_file,beta.info,beta.all)
  
  if(1){
    
    test.bfile <- test_file
    B.beta.info <- beta.info
    B.beta.all <- beta.all
    
    test.bim=fread(paste0(test.bfile,".bim"),header=F)
    test.bim$V1 <- as.character(sub("^chr", "", test.bim$V1, ignore.case = T))
    test.bim$.index.ref <- 1:nrow(test.bim)
    
    # > test.bim
    # V1              V2 V3       V4 V5 V6 .index.ref
    # 1:  1   1:1962845:T:C  0  1962845  C  T          1
    # 2:  1   1:1962899:A:C  0  1962899  C  A          2
    # 3:  1   1:1963406:G:A  0  1963406  A  G          3
    # 4:  1   1:1963538:T:C  0  1963538  C  T          4
    # 5:  1   1:1963738:C:T  0  1963738  T  C          5
    # ---                                                
    #   5630741: 22 22:50797551:G:A  0 50797551  A  G    5630741
    # 5630742: 22 22:50798021:A:G  0 50798021  G  A    5630742
    # 5630743: 22 22:50798635:T:C  0 50798635  C  T    5630743
    # 5630744: 22 22:50801260:T:C  0 50801260  C  T    5630744
    # 5630745: 22 22:50802392:C:T  0 50802392  T  C    5630745
    
    B.beta.info$.index.tomatch <- 1:nrow(B.beta.info)
    
    merged <- merge(test.bim,   B.beta.info,all=F, 
                    by.x="V2", by.y="SNP",sort=F)
    
   
    # > merged[V2 == '1:1991822:G:A',]
    # V2 V1 V3      V4 V5 V6 .index.ref A1 .index.tomatch
    # 1: 1:1991822:G:A  1  0 1991822  A  G        112  A             98
    # > which.max(merged$V2 == '1:1991822:G:A')
    # [1] 98
    # > B.beta.all[98,1:5]
    # betatemp   betatemp.1  betatemp.2   betatemp.3  betatemp.4
    # 98 -0.00641277 -0.006412755 -0.00641274 -0.006412725 -0.00641271
    
    ###check the sign##
    # list1=which(merged$A2==merged$V5)
    list1=which(merged$A1!=merged$V5)
    # list1=which(merged$A1==merged$V5)
    if (length(list1)>0){
      B.beta.all[list1,]=-B.beta.all[list1,]
    }
    
    flag<-rep(FALSE, nrow(test.bim))
    flag[merged$.index.ref]<-TRUE
    
    BM=as.matrix(B.beta.all[merged$.index.tomatch,])
    BM0=list()
    BM0[[1]]=BM
    
    system.time({
      pgs <- lapply(BM0, function(x) lassosum:::pgs(bfile=test.bfile, weights = x, 
                                                    extract=flag, keep=NULL, 
                                                    cluster=NULL))
    }) ####is this correct for my case? TZ
    
    fam=read.table(paste0(test.bfile,".fam"))[,c(1,2)]
    colnames(fam)[1:2]=c("FID","IID")
    out.all=cbind(fam,pgs[[1]])
  }#Calculate_PRS
  
}##PRStr_tuning


system.time(out.beta <- TL_PRS(ped_file = ped_file,
                               Covar_name = NULL,
                               Y_name = Y_name, 
                               Ytype= Ytype, 
                               train_file = train_file,
                               test_file = validate_file,
                               sum_stats_file = sum_stats_file,
                               target_sumstats_file = target_sumstats_file,
                               LDblocks = LDblocks,
                               outfile = outfile,
                               cluster= cl,
                               plink = plink))
# stopCluster(cl) 

summary(out.beta)
save(out.beta, file = 'TL.first.out.RData')
write.table(out.beta, file='out.beta.txt',quote=F,row.names=F,col.names=F)


