##cor=bim_sum_stats[which(LDblocks2[[1]]==i),]; num=which(i==unique(LDblocks2[[1]]));nsnp=nrow(bim_sum_stats)
block_calculation2<-function(cor,num,train_file,nsnp,temp.file, plink){
  temp_file=paste0(temp.file,"_block_",num)
  write.table(cor,file=temp_file,col.names=F,row.names=F,quote=F)
  cmd = paste0(plink, " --bfile ",train_file," --extract ",temp_file," --keep-allele-order"," --recode A  --out ", temp_file,"_Geno.txt")
  system(cmd)
  
  Gtemp=try(as.data.frame(fread(paste0(temp_file,"_Geno.txt.raw"),header=T)),silent=T)
  if (file.exists(temp_file)) {file.remove(temp_file)}
  if (file.exists(paste0(temp_file,"_Geno.txt.nosex"))) {file.remove(paste0(temp_file,"_Geno.txt.nosex"))}
  if (file.exists(paste0(temp_file,"_Geno.txt.log"))) {file.remove(paste0(temp_file,"_Geno.txt.log"))}
  if (file.exists(paste0(temp_file,"_Geno.txt.raw"))) {file.remove(paste0(temp_file,"_Geno.txt.raw"))}
  
  if (class(Gtemp)=="try-error"){
    return(NULL)
    #GG=diag(nrow(cor));colnames(GG)=paste0(cor$V2,"_",cor$V5) 
    #geno_info=as.data.frame(t(sapply(colnames(GG),   split_SNPandA1   ) )) ;colnames(geno_info)=c("SNP","A1")
    #geno_info$mean=NA; geno_info$maf=NA; geno_info$sd=NA
  }else{
    GG=cor(as.matrix(Gtemp[,7:ncol(Gtemp)]))
    geno_info=as.data.frame(t(sapply(colnames(Gtemp)[7:ncol(Gtemp)],   split_SNPandA1   ) )) ;colnames(geno_info)=c("SNP","A1")
    geno_info$mean=colMeans(as.matrix(Gtemp[,7:ncol(Gtemp)]),na.rm=T); geno_info$maf=geno_info$mean/2; geno_info$sd=sqrt(2*geno_info$maf*(1-geno_info$maf))
    
  }
  
  list1=which(geno_info$sd==0)
  if (length(list1)>0){
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
    GG2=as.matrix(GG[geno_info2$order,geno_info2$order])
    gy=geno_info2$cor
    betatemp=geno_info2$Beta2*geno_info2$sd
    u0=gy-GG2%*%betatemp
    beta.all=cbind(u0, betatemp)
    # for (factor1 in c(1,10,100,1000)){
    for (factor1 in c(10,1e3,1e4)){
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
    return(geno_info2)
  } 
}##function end




##PRStr_calculation2(sum_stats_target, train_file, sum_stats, LDblocks, cluster=cluster,temp.file=paste0(tempfile,"_step1"))
##temp.file=paste0(tempfile,"_step1")
PRStr_calculation2<-function(sum_stats_target, train_file, sum_stats, LDblocks, cluster=NULL,temp.file,
                             plink){
  possible.LDblocks <- c("EUR.hg19", "AFR.hg19", "ASN.hg19", 
                         "EUR.hg38", "AFR.hg38", "ASN.hg38") 
  if(!is.null(LDblocks)) {
    if(is.character(LDblocks) && length(LDblocks) == 1) {
      if(LDblocks %in% possible.LDblocks) {
        LDblocks <- data.table::fread(system.file(paste0("data/Berisa.",  LDblocks, ".bed"),  package="lassosum"), header=T)
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
  flag1=which(bim_sum_stats$V5==bim_sum_stats$A1)
  if (length(flag1)>0){  bim_sum_stats$Beta2[flag1]=bim_sum_stats$Beta[flag1]}
  flag2=which(bim_sum_stats$V6==bim_sum_stats$A1)
  if (length(flag2)>0){  bim_sum_stats$Beta2[flag2]=-bim_sum_stats$Beta[flag2];  bim_sum_stats$cor[flag2]=-bim_sum_stats$cor[flag2];}
  
  bim_sum_stats=bim_sum_stats[which(! is.na(bim_sum_stats$Beta2)),c("V2","V1","V4","V5","V6","order","Beta2","cor")]
  
  
  ref.extract <- rep(FALSE, nrow(ref.bim))
  ref.extract[bim_sum_stats$order] <- TRUE
  
  
  if(!is.null(LDblocks)) {
    LDblocks2 <- splitgenome2(CHR = ref.bim$V1[ ref.extract], 
                              POS = ref.bim$V4[ ref.extract],
                              ref.CHR = LDblocks[,1], 
                              ref.breaks = LDblocks[,3])
    # Assumes base 1 for the 3rd column of LDblocks (like normal bed files)
  } 
  
  if(is.null(cluster)) {
    results.list <- lapply(unique(LDblocks2[[1]]), function(i) {
      block_calculation2(cor=bim_sum_stats[which(LDblocks2[[1]]==i),], num=which(i==unique(LDblocks2[[1]])),train_file=train_file,nsnp=nrow(bim_sum_stats),temp.file,
                         plink)
    })
  } else {
    results.list <-  parallel::parLapplyLB(cluster,unique(LDblocks2[[1]]), function(i) {
      block_calculation2(cor=bim_sum_stats[which(LDblocks2[[1]]==i),], num=which(i==unique(LDblocks2[[1]])),train_file=train_file,nsnp=nrow(bim_sum_stats),temp.file,
                         plink)
    })
  }
  
  results.list<-do.call("rbind", results.list)
  
  return(results.list)
}



######################The function used summary statistics for training############### 
##ped_file=ped.file;Covar_name="";Y_name=kword;Ytype="C"; train_file=train.bfile;test_file=test.bfile;sum_stats_file=beta.file;LDblocks="EUR.hg19"
##ped.file,"",kword, Ytype="C",train.bfile,test.bfile,beta.file,target_sumstats_file,LDblocks="EUR.hg19",tempfile
TL_PRS<-function(ped_file,Covar_name,Y_name, Ytype="C",train_file,test_file,sum_stats_file,target_sumstats_file, LDblocks="EUR.hg19",outfile,cluster=NULL,
                 plink = "plink-1.9"){
  tempfile=outfile
  out1=PRStr_main_check(ped_file,Covar_name,Y_name, Ytype,train_file,test_file,sum_stats_file,LDblocks)
  if (out1!=0){stop(out1)}
  
  sum_stats=data.frame(fread(sum_stats_file))
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
  
  sum_stats_target$cor=lassosum::p2cor(p = sum_stats_target$p, n = median(sum_stats_target$N,na.rm=T), sign=sum_stats_target$beta)
  flag=which(sum_stats_target$A1.x !=sum_stats_target$A1.y)
  if (length(flag)>0){sum_stats_target$cor[flag]=-sum_stats_target$cor[flag]}
  sum_stats_target=sum_stats_target[,c("SNP","A1.x","Beta","cor")];colnames(sum_stats_target)[2]="A1";
  gc()
  
  beta_list=as.data.frame(PRStr_calculation2(sum_stats_target, train_file, sum_stats, LDblocks, cluster=cluster,temp.file=paste0(tempfile,"_step1"),
                                             plink = plink))
  fwrite(beta_list,file=paste0(tempfile,"_beta.candidates.fresh.txt"),row.names=F,quote=F,col.names=T)
  
  return(1)
  # beta_list=as.data.frame(beta_list[,-c(5,9)]) ###you got rid of the correct A1?!?!
  # ##the fifth column is the A1 and the nineth is the order.
  
  # beta_list[,2] <- beta_list[,5] ###now the second column is the correct A1
  # beta_list <- as.data.frame(beta_list)
  # 
  # colnames(beta_list)[1:2]=c("SNP","A1")
  # write.table(beta_list,file=paste0(tempfile,"_beta.candidates.txt"),row.names=F,quote=F,col.names=T)
  
  
  # 	out1=PRStr_tuning(beta_list, ped,Covar_name, Y_name, Ytype,test_file)
  # 
  #   	if (file.exists(paste0(tempfile,"_original_sum_stats.txt"))) {file.remove(paste0(tempfile,"_original_sum_stats.txt"))}
  #   	if (file.exists(paste0(tempfile,"_step0.train.PRS.nosex"))) {file.remove(paste0(tempfile,"_step0.train.PRS.nosex"))}
  #   	if (file.exists(paste0(tempfile,"_step0.train.PRS.log"))) {file.remove(paste0(tempfile,"_step0.train.PRS.log"))}
  #   	if (file.exists(paste0(tempfile,"_step0.train.PRS.profile"))) {file.remove(paste0(tempfile,"_step0.train.PRS.profile"))}  
  # 	write.table(out1$best.beta,file=paste0(tempfile,"_best.beta.txt"),row.names=F,quote=F,col.names=T)
  # 	write.table(out1$best.PRS,file=paste0(tempfile,"_best.PRS.txt"),row.names=F,quote=F,col.names=T)
  # 
  # 	return(out1)
}

# 
# AUC.from.beta.candidates <- function(beta.can.file,
#                                      test_file, stats_file, ped_test_file,Y_name){
#   beta.all <- fread(beta.can.file) ##read in all the beta
# 
#   ##normalize and make sure the signs are correct
#   beta.all[,11:ncol(beta.all)] <- beta.all[,11:ncol(beta.all)] / beta.all$sd
#   not.same.A1 <- which(beta.all$V5 != beta.all$A1)
#   beta.all[not.same.A1,11:ncol(beta.all)] <- -beta.all[not.same.A1,11:ncol(beta.all)]
#   names(beta.all)[1:2] <- c("SNP","A1")
# 
#   Y_test=read.table(ped_test_file,header=T, sep = ' ')
# 
#   beta.all <- as.data.frame(beta.all)
#   aucs <- rep(0, length(12:ncol(beta.all)))
#   k <- 1
#     for(i in 12:ncol(beta.all)){
#     temp.beta <- beta.all[, c(1,2,i)]
#     fwrite(temp.beta,
#                 file = stats_file, sep = ' ',
#                 col.names = FALSE)
#     cmd <- paste0("/usr/local/bin/plink --bfile ",test_file,
#                "  --score ",stats_file, " sum",
#                " --out ",stats_file,".test.PRS",
#                " --allow-no-sex")
#     system(cmd)
#     print("0")
#     temp=read.table(paste0(stats_file,".test.PRS.profile"),header=T)
#     print("2")
#     merged=merge(temp,Y_test,by=c("FID","IID"))
#     print("3")
#     aucs[k] <- auc(merged[,Y_name],merged$SCORESUM)
#     k <- k+1
#     print(aucs)
#     }
#   
#   return(aucs)
#   # > aucs
#   # [1] 0.7124512 0.7124523 0.7124569 0.7124571 0.7124580 0.7124589 0.7124601
#   # [8] 0.7124647 0.7124666 0.7124677 0.7124701 0.7124725 0.7124764 0.7124796
#   # [15] 0.7124797 0.7124677 0.7124936 0.7125233 0.7125501 0.7125733 0.7125966
#   # [22] 0.7126212 0.7126488 0.7126644 0.7126870 0.7127111 0.7127335 0.7127590
#   # [29] 0.7127785 0.7128108 0.7126876 0.7129319 0.7131436 0.7133729 0.7136186
#   # [36] 0.7138386 0.7140674 0.7142802 0.7145073 0.7147295 0.7149423 0.7151689
#   # [43] 0.7153975 0.7155980 0.7158003 0.7147240 0.7168759 0.7187756 0.7206672
#   # [50] 0.7223919 0.7239199 0.7254027 0.7266540 0.7278350 0.7288401 0.7297561
#   # [57] 0.7305475 0.7312015 0.7318841 0.7323994
#   }
# 
# test_file <- paste0(work.dir, 'YRI.TUNE/YRI.TUNE')
# ped_test_file <- paste0(work.dir, '/YRI.TUNE/YRI-TUNE-fam-forTL.txt')
# Y_name <- "Y"
# stats_file <- paste0(main.dir,'firstTL_temp.beta.txt')
# beta.can.file <- paste0(main.dir,"firstTL_beta.candidates.fresh.txt")
# 
# tune.auc <- AUC.from.beta.candidates(beta.can.file,
#                                      test_file, stats_file, ped_test_file, Y_name)
# 
# AUC.best.beta <- function(beta.can.file,
#                           test_file, stats_file, ped_test_file,Y_name,
#                           best.beta.index){
#   beta.all <- fread(beta.can.file) ##read in all the beta
#   
#   ##normalize and make sure the signs are correct
#   beta.all[,11:ncol(beta.all)] <- beta.all[,11:ncol(beta.all)] / beta.all$sd
#   not.same.A1 <- which(beta.all$V5 != beta.all$A1)
#   beta.all[not.same.A1,11:ncol(beta.all)] <- -beta.all[not.same.A1,11:ncol(beta.all)]
#   names(beta.all)[1:2] <- c("SNP","A1")
#   
#   Y_test=read.table(ped_test_file,header=T, sep = ' ')
#   
#   beta.all <- as.data.frame(beta.all)
#   best.auc <- -1
#   
#   k <- 1
#   for(i in 12:ncol(beta.all)){
#     if(k == best.beta.index){
#     temp.beta <- beta.all[, c(1,2,i)]
#     fwrite(temp.beta,
#            file = stats_file, sep = ' ',
#            col.names = FALSE)
#     cmd <- paste0("/usr/local/bin/plink --bfile ",test_file,
#                   "  --score ",stats_file, " sum",
#                   " --out ",stats_file,".test.PRS",
#                   " --allow-no-sex")
#     system(cmd)
#     print("0")
#     temp=read.table(paste0(stats_file,".test.PRS.profile"),header=T)
#     print("2")
#     merged=merge(temp,Y_test,by=c("FID","IID"))
#     print("3")
#     best.auc <- auc(merged[,Y_name],merged$SCORESUM)
#     k <- k+1
#     print(best.auc)
#     }##if best beta
#   }
#   
#   return(aucs)
# }
# 
# test_file <- paste0(work.dir, 'YRI.TST/YRI.TST')
# ped_test_file <- paste0(work.dir, '/YRI.TST/YRI-TST-fam-forTL.txt')
# Y_name <- "Y"
# stats_file <- paste0(main.dir,'firstTL_temp.beta.txt')
# beta.can.file <- paste0(main.dir,"firstTL_beta.candidates.fresh.txt")
# 
# test.auc <- AUC.best.beta(beta.can.file,
#                           test_file, stats_file, ped_test_file, Y_name,
#                           which.max(tune.auc))
# > test.aucs
# [1] 0.7124512 0.7124523 0.7124569 0.7124571 0.7124580 0.7124589 0.7124601
# [8] 0.7124647 0.7124666 0.7124677 0.7124701 0.7124725 0.7124764 0.7124796
# [15] 0.7124797 0.7147240 0.7168759 0.7187756 0.7206672 0.7223919 0.7239199
# [22] 0.7254027 0.7266540 0.7278350 0.7288401 0.7297561 0.7305475 0.7312015
# [29] 0.7318841 0.7323994 0.7289631 0.7335943 0.7304489 0.7235622 0.7153189
# [36] 0.7070783 0.6993625 0.6923570 0.6860880 0.6805120 0.6755146 0.6710677
# [43] 0.6670146 0.6633822 0.6601632 0.6785349 0.6471623 0.6324021 0.6242296
# [50] 0.6185533 0.6144268 0.6112639 0.6086980 0.6064754 0.6046272 0.6030144
# [57] 0.6015476 0.6002817 0.5991710 0.5981460