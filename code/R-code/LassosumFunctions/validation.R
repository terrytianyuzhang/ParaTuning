### Validation function ###
validate<-function(fit, fileName, extract=NULL, keep=NULL, distance=NULL, chunk=T, mem.limit=4*10^9,...){

  parsed<-parseselect(fileName, extract=extract, keep=keep, distance=distance)
  ### Split/group into chunks###
  if(chunk==T){
    required.memory=as.numeric(parsed$p) * as.numeric(parsed$n) * 8 
    n.chunks<-ceiling(required.memory/mem.limit)
    print(paste0("Number of chunks: ",n.chunks))
    
    chunk.len=ceiling(parsed$p/n.chunks)
    print(paste0("chunk length: ",chunk.len))
    
    if(n.chunks > 1) {
      if(is.null(parsed$extract)){
        extract=lapply(1:n.chunks, function(i) logical.vector(min(((i-1)*chunk.len + 1), parsed$P) : min((i*chunk.len), parsed$P), parsed$P))
        fit.list<-vector("list",length=n.chunks)
        for(i in 1:n.chunks) {
          fit.list[[i]] = fit
          fit.list[[i]]$beta=fit.list[[i]]$beta[extract[[i]],,drop=FALSE]
          }
      } else{
        extract=lapply(1:n.chunks, function(i) logical.vector(which((cumsum(parsed$extract)<=i*chunk.len) & (cumsum(parsed$extract) > (i-1)*chunk.len)),parsed$P))
        fit.list<-vector("list",length=n.chunks)
        for(i in 1:n.chunks) {
          fit.list[[i]] = fit
          fit.list[[i]]$beta = fit.list[[i]]$beta[which(parsed$extract) %in% which(extract[[i]]),,drop=FALSE]
          }
        }

      results.list <- lapply(1:n.chunks, function(i) {
        print(i)
        validate(fit=fit.list[[i]], fileName = fileName, extract=extract[[i]], chunk = FALSE, distance=NULL)
      })
      
      mysum <- function(...) {
        mat <- do.call("rbind", list(...))
        if(ncol(mat) > 0) return(as.vector(colSums(mat))) else
          return(numeric(0))
      }
      
      return(list(pgscore=do.call("mysum", lapply(results.list, function(x) as.vector(x$pgscore)))))
    }
  }
  #extract phenotypes
  phcovar <- parse.pheno.covar(pheno=NULL, covar=NULL, parsed=parsed, trace=2)
  pheno=phcovar$pheno
  pgscore=matrix(nrow=parsed$n, ncol=length(fit$lambda))
  
  if(is.null(parsed$extract)) {
    extract <- list(integer(0), integer(0))
    #col_skip_pos = integer(0); col_skip = integer(0)
  } else {
    #print(parsed$extract)
    extract <- selectregion(!parsed$extract)
    extract[[1]] <- extract[[1]] - 1
  }
  
  if(is.null(parsed$keep)) {
    keepbytes <- integer(0)
    keepoffset <- integer(0)
  } else {
    pos <- which(parsed$keep) - 1
    keepbytes <- floor(pos/4)
    keepoffset <- pos %% 4 * 2
  }
  #beta.id= which(extract.all) %in% which(parsed$extract) 
  #shrink=fit$shrink
  
  #r2=c()
  #AUC=c()
  library(pROC)
  for(i in 1:ncol(fit$beta)){
    #x1=x[,i,drop=F]
    #sd=fit$sd
    
    #pgs[,i]=(scale(test.data$G)/(sqrt(N-1)*sqrt(1-shrink))) %*% x1
    pgscore[,i]=pgs(fileName=fileName, N=parsed$N, P=parsed$P, shrink=fit$shrink, x=fit$beta[,i,drop=F], 
                    col_skip_pos=extract[[1]], col_skip=extract[[2]], keepbytes=keepbytes, keepoffset=keepoffset, 2)
    #r2[i]=cor(pgscore[,i], pheno, use="complete.obs") #remove missing values
    #AUC[i]=auc(roc(pheno, pgscore[,i]))
  }
  return(list(pgscore=pgscore
              #r2=r2,
              #AUC=AUC
              ))
}
