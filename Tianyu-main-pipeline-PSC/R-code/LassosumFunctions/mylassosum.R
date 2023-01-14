mylassosum <- function(cor1, cor2,  fileName1, fileName2,
                       gamma=0.5,
                       lambda=exp(seq(log(0.001), log(0.05), length.out=5)), 
                       shrink=0.9, 
                       #ridge=F, thr= 10^(-4), 
                       chunk=T,
                       init=NULL, trace=0, maxiter=10000,
                       keep1=NULL, keep2=NULL, remove=NULL, extract1=NULL, extract2=NULL, exclude=NULL, distance=NULL,
                       chr=NULL, mem.limit=4*10^9,
                       LDblocks1=NULL, LDblocks2=NULL) {
  
  cor1 <- as.vector(cor1)
  cor2 <- as.vector(cor2)
  stopifnot(!any(is.na(cor1)))
  stopifnot(!any(is.na(cor2)))
  stopifnot(length(cor1) == length(cor2))

  print(paste0("distance is NULL:", is.null(distance)))
  parsed1<-parseselect(fileName1, extract=extract1, exclude=exclude, keep=keep1, remove=remove, distance=distance, chr=chr) #extract1 could be a logical vector? 
  parsed2<-parseselect(fileName2, extract=extract2, exclude=exclude, keep=keep2, remove=remove, distance=distance, chr=chr)
  
  ### Split/group into chunks###
  if(chunk==T){
    required.memory=as.numeric(parsed1$p) * as.numeric(parsed1$n) * 8 + as.numeric(parsed2$p) * as.numeric(parsed2$n) *8
    n.chunks<-ceiling(required.memory/mem.limit)
    print(paste0("Number of chunks: ",n.chunks))
    
    chunk.len1=ceiling(parsed1$p/n.chunks)
    chunk.len2=ceiling(parsed2$p/n.chunks)
    stopifnot(chunk.len1 == chunk.len2)
    chunk.len=chunk.len1
    print(paste0("chunk length: ",chunk.len))
    
    if(n.chunks > 1) {
      if(is.null(parsed1$extract)){
        extract1=lapply(1:n.chunks, function(i) logical.vector(min(((i-1)*chunk.len + 1), parsed1$P) : min((i*chunk.len), parsed1$P), parsed1$P))
      } else{
        extract1=lapply(1:n.chunks, function(i) logical.vector(parsed1$extract & (cumsum(parsed1$extract)<=i*chunk.len) & (cumsum(parsed1$extract) > (i-1)*chunk.len),parsed1$P))
      }
      if(is.null(parsed2$extract)){
        extract2=lapply(1:n.chunks, function(i) logical.vector(min(((i-1)*chunk.len + 1), parsed2$P) : min((i*chunk.len), parsed2$P), parsed2$P))
      } else{
        extract2=lapply(1:n.chunks, function(i) logical.vector(parsed2$extract & (cumsum(parsed2$extract)<=i*chunk.len) & (cumsum(parsed2$extract) > (i-1)*chunk.len),parsed2$P))
      }

      results.list <- lapply(1:n.chunks, function(i) {
        mylassosum(cor1=cor1, cor2=cor2, fileName1=fileName1, fileName2=fileName2, 
                   gamma=gamma, lambda=lambda, shrink=shrink, 
                   keep1=parsed1$keep, keep2=parsed2$keep, extract1=extract1[[i]],extract2=extract2[[i]], distance = NULL,
                   trace=trace, maxiter=maxiter, chunk=FALSE,
                   LDblocks1=LDblocks1, LDblocks2=LDblocks2,  
                   mem.limit=mem.limit)
      })
      print(paste0("length(results.list): ", length(results.list)) )
      print(paste0("class: ", sapply(results.list, "class") ))
      
      re.lassosum=merge.mylassosum(results.list)
#      return(do.call("merge.mylassosum", results.list))
      return(re.lassosum)
    }## n.chunks > 1
  }##chunk == T
  
  if(!is.null(parsed1$extract) & length(parsed1$extract) == length(cor1)) cor1=cor1[parsed1$extract]
  if(!is.null(parsed2$extract) & length(parsed2$extract) == length(cor2)) cor2=cor2[parsed2$extract]
  ### Split by LD region ###
  reindex <- function(x) {
    as.integer(x)-min(as.integer(x)) + 1
  }

  bimfile1=read.table2(paste0(fileName1,".bim"))
  bimfile2=read.table2(paste0(fileName2,".bim"))
  
  bimfile1$V1 = paste0("chr",bimfile1$V1)
  bimfile2$V1 = paste0("chr",bimfile2$V1)
  
#  if(is.null(extract)){
#    INDEX1 = factor(bimfile1$V1, levels=unique(bimfile1$V1))
#    INDEX2 = factor(bimfile2$V1, levels=unique(bimfile2$V1))
#  } else {
#    INDEX1 = factor(bimfile1$V1[parsed1$extract], levels=unique(bimfile1$V1[parsed1$extract]))
#    INDEX2 = factor(bimfile2$V1[parsed2$extract], levels=unique(bimfile2$V1[parsed2$extract]))
#  }


  if(is.null(extract1)){
    blocks1 <- splitgenome(CHR = bimfile1$V1, 
                           POS = bimfile1$V4,
                           ref.CHR = LDblocks1[,1], 
                           ref.breaks = LDblocks1[,3])
  } else {
    blocks1 <- splitgenome(CHR = bimfile1$V1[parsed1$extract], 
                           POS = bimfile1$V4[parsed1$extract],
                           ref.CHR = LDblocks1[,1], 
                           ref.breaks = LDblocks1[,3])
  }
  
  if(is.null(extract2)){
    blocks2 <- splitgenome(CHR = bimfile1$V1, 
                           POS = bimfile1$V4,
                           ref.CHR = LDblocks2[,1], 
                           ref.breaks = LDblocks2[,3])
  } else {
    blocks2 <- splitgenome(CHR = bimfile1$V1[parsed2$extract], 
                           POS = bimfile1$V4[parsed2$extract],
                           ref.CHR = LDblocks2[,1], 
                           ref.breaks = LDblocks2[,3])
  }
  
  #blocks1<-as.numeric(unlist(tapply(blocks1, INDEX1, reindex)))
  #blocks1<-as.integer(blocks1)-min(as.integer(blocks1)) + 1
  blocks1<-as.numeric(factor(blocks1, levels=unique(blocks1)))
  blocks2<-as.numeric(factor(blocks2, levels=unique(blocks2)))
  #blocks2<-as.integer(blocks2)-min(as.integer(blocks2)) + 1
  #blocks2<-as.numeric(unlist(tapply(blocks2, INDEX2, reindex)))
  
  el <- myelnet(lambda, shrink, gamma, fileName1, fileName2, cor1, cor2, 
                 #thr=thr,
                 parsed1, parsed2,
                 trace=trace, maxiter=maxiter, 
                 blocks1=blocks1, blocks2=blocks2)
  
  #--------------------------------------------
  #This section is not used yet
  #if(ridge) {
  #  if(is.null(blocks)) blocks <- rep(1, p)
  #  reps <- 1:max(blocks)
  #  svd <- lapply(1:reps, function(i) svd(X[,blocks==i], nu=0))
  #  invert <- lapply(1:reps, function(i) 1/(svd[[i]]$d^2 + shrink))
  #  VG <- lapply(1:reps, function(i) svd[[i]]$v %*% Diagonal(x=invert[[i]] - 1/shrink))
  #  VTr <- lapply(1:reps, function(i) t(svd[[i]]$v) %*% cor[blocks == i])
  #  VGVTr <- lapply(1:reps, function(i) VG[[i]] %*% VTr[[i]])
  #  Ridge <- lapply(1:reps, function(i) as.vector(VGVTr[[i]] + 1/shrink * cor[blocks==i]))
  #  Ridge <- unlist(Ridge)
  #} else Ridge <- NULL
  #----------------------------------------------
  
  nparams <- colSums(el$beta != 0) 
  
  toreturn <- list(lambda=lambda, 
                   beta=el$beta,
                   conv=el$conv,
                   pred1=el$pred1,
                   pred2=el$pred2,
                   loss=el$loss,
                   trainerror1 = el$trainerror1,
                   trainerror2 = el$trainerror2,
                   fbeta=el$fbeta,
                   sd1=el$sd1,
                   sd2=el$sd2,
                   shrink=shrink,
                   nparams=nparams, 
                   gamma=gamma
                   #,ridge=Ridge
                   )
  
  class(toreturn) <- "lassosum"
  
  print('loss in mylassosum.R')
  print(toreturn$trainerror1)
  print(toreturn$trainerror2)
  
  return(toreturn)
  
  #' @return A list with the following
  #' \item{lambda}{same as the lambda input}
  #' \item{beta}{A matrix of estimated coefficients}
  #' \item{conv}{A vector of convergence indicators. 1 means converged. 0 not converged.}
  #' \item{pred}{\eqn{=(1-s)X\beta}}
  #' \item{loss}{\eqn{=(1-s)\beta'X'X\beta/n - 2\beta'r}}
  #' \item{trainerror1}{training error for population 1}
  #' \item{trainerror2}{training error for population 2}
  #' \item{fbeta}{\eqn{=\beta'R\beta - 2\beta'r + 2\lambda||\beta||_1}}
  #' \item{sd}{The standard deviation of the reference panel SNPs}
  #' \item{shrink}{same as input}
  #' \item{nparams}{Number of non-zero coefficients}
  #' \item{gamma}{A parameter balancing two populations}
  #' \item{ridge}{ridge regression estimates}
  
}
