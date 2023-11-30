### Load functions ###
#' @title Function to read a text file 
#' @keywords internal
read.table2 <- function(file, header=F, data.table=F, check.names=TRUE, ...) {
  return(data.table::fread(file, header=header, data.table=data.table, 
                           check.names=check.names, ...))
}


#' @title Parse the keep/remove/extract/exclude/chr options
#' @details \code{keep} and \code{remove} can take one of three
#' formats:
#' 1. A logical vector indicating which indivduals to keep/remove,
#'
#' 2. A \code{data.frame} with two columns giving the FID and IID of the indivdiuals to keep/remove (matching those in the .fam file), or
#' 3. A character scalar giving the text file with the FID/IID.
#'
#' Note that these files should not contain headers.
#'
#' \code{extract} and \code{exclude} can be of format type __(1)__ or __(3)__ describe above, or a character vector of SNP ids (matching those in the .bim file).

#' @param bfile plink file stem
#' @param extract SNPs to extract
#' @param exclude SNPs to exclude
#' @param keep samples to keep
#' @param remove samples to remove
#' @param chr a vector of chromosomes
#' @param export Logical. Whether to include the \code{bim} and/or \code{fam} data in the returned object.
#' @param order.important Logical. Whether to strictly enforce that the order of \code{extract} and/or \code{keep} must match with \code{bim} or \code{fam}, respectively.
#'
#' @md
#' @export
parseselect <- function(bfile, extract=NULL, exclude=NULL,
                        keep=NULL, remove=NULL, chr=NULL,
                        export=FALSE, distance=NULL, order.important=FALSE) {
  
  #### Introduced for parsing multiple bfiles ####
  if(length(bfile) > 1) {
    if(!is.null(chr) || !is.null(extract) || !is.null(exclude)) {
      stop("bfile cannot be a vector if we are parsing extract/exclude/chr. yet.")
    }
    bfile <- bfile[1]
  }
  
  #### Checks ####
  stopifnot(is.character(bfile) && length(bfile) == 1)
  bedfile <- paste0(bfile, ".bed")
  bimfile <- paste0(bfile, ".bim")
  famfile <- paste0(bfile, ".fam")
  stopifnot(file.exists(bedfile))
  stopifnot(file.exists(bimfile))
  stopifnot(file.exists(famfile))
  
  if(grepl("^~", bfile)) {
    stop("Don't use '~' as a shortcut for the home directory.")
  }
  
  p <- P <- ncol.bfile(bfile)
  n <- N <- nrow.bfile(bfile)
  bim <- NULL
  fam <- NULL
  
  #### extract ####
  if(!is.null(extract)) {
    if(is.logical(extract)) {
      stopifnot(length(extract) == P)
    } else {
      if(is.null(attr(extract, "not.a.file", exact=T))) {
        if(is.character(extract) && length(extract) == 1 && file.exists(extract)) {
          ### I'm interpreting this as a filename
          SNPs <- read.table2(extract)
          stopifnot(ncol(SNPs)==1)
          extract <- SNPs[[1]]
        }
      } else {
        attr(extract, "not.a.file") <- NULL
      }
      if(is.vector(extract)) {
        if(is.null(distance)){
          Extract <- as.character(extract)
          bim <- read.table2(bimfile)
          extract <- bim$V2 %in% Extract
          if(order.important) {
            if(!all(bim$V2[extract] == Extract)) {
              stop("Order of extract SNPs does not match that in .bim file.")
            }
          }
        } else if(distance>0){
          Extract <- as.character(extract)
          bim <- read.table2(bimfile)
          bim.extract <- bim[bim$V2 %in% Extract,]
          bim.extract$start = bim.extract$V4 - distance
          bim.extract$start = ifelse(bim.extract$start<0, 0, bim.extract$start)
          bim.extract$end = bim.extract$V4 + distance
          extract<-rep(FALSE,nrow(bim))
          for(l in 1:nrow(bim.extract)){
            extract[bim$V4 > bim.extract$start[l] & bim$V4 < bim.extract$end[l] & bim.extract$V1[l]==bim$V1]<-TRUE
          }
        } else if(distance==0){
          Extract <- as.character(extract)
          bim <- read.table2(bimfile)
          extract <- bim$V2 %in% Extract
          if(order.important) {
            if(!all(bim$V2[extract] == Extract)) {
              stop("Order of extract SNPs does not match that in .bim file.")
            }
          }
        } else{
          stop("distance should be a non negative number ")
        }
      } else {
        stop("I don't know what to do with this type of input for extract")
      }
    }
    
    p <- sum(extract)
  }
  
  #### exclude ####
  if(!is.null(exclude)) {
    if(is.logical(exclude)) {
      stopifnot(length(exclude) == P)
    } else {
      if(is.null(attr(exclude, "not.a.file", exact=T))) {
        if(is.character(exclude) && length(exclude) == 1 && file.exists(exclude)) {
          ### I'm interpreting this as a filename
          SNPs <- read.table2(exclude)
          stopifnot(ncol(SNPs)==1)
          exclude <- SNPs[[1]]
        }
      } else {
        attr(exclude, "not.a.file") <- NULL
      }
      if(is.vector(exclude)) {
        exclude <- as.character(exclude)
        if(is.null(bim)) bim <- read.table2(bimfile)
        exclude <- bim$V2 %in% exclude
      } else {
        stop("I don't know what to do with this type of input for exclude")
      }
      
    }
    
    if(is.null(extract)) extract <- !exclude else
      extract <- extract & !exclude
    
    p <- sum(extract)
  }
  
  #### chr ####
  if(!is.null(chr)) {
    
    stopifnot(is.vector(chr))
    chr <- as.character(chr)
    
    if(is.null(bim)) bim <- read.table2(bimfile)
    bimchr <- bim$V1
    bimchr[bimchr==""]
    extract.chr <- bim$V1 %in% chr
    
    if(is.null(extract)) extract <- extract.chr else
      extract <- extract & extract.chr
    
    p <- sum(extract)
    
  }
  
  #### keep ####
  if(!is.null(keep)) {
    if(is.logical(keep)) {
      stopifnot(length(keep) == N)
    } else {
      if(is.character(keep) && length(keep) == 1 && file.exists(keep)) {
        ### I'm interpreting this as a filename
        keep <- read.table2(keep)
      }
      if(is.vector(keep)) {
        #keep <- as.data.frame(keep)
        #stopifnot(ncol(keep)==2)
        fam <- read.table2(famfile)
        famID <- fam[,1]
        #keepID <- keep
        keep <- famID %in% keep
        if(order.important) {
          if(!all(famID[keep] == keepID)) {
            stop("Order of keep doesn't match that in .fam file")
          }
        }
      } else {
        stop("I don't know what to do with this type of input for keep")
      }
      
    }
    n <- sum(keep)
  }
  
  #### remove ####
  if(!is.null(remove)) {
    if(is.logical(remove)) {
      stopifnot(length(remove) == N)
    } else {
      if(is.character(remove) && length(remove) == 1 && file.exists(remove)) {
        ### I'm interpreting this as a filename
        remove <- read.table2(remove)
      }
      if(inherits(remove, "data.frame")) {
        remove <- as.data.frame(remove)
        stopifnot(ncol(remove)==2)
        if(is.null(fam)) fam <- read.table2(famfile)
        famID <- paste(fam[,1], fam[,2], sep=".")
        removeID <- paste(remove[,1], remove[,2], sep=".")
        remove <- famID %in% removeID
      } else {
        stop("I don't know what to do with this type of input for remove")
      }
    }
    
    if(is.null(keep)) keep <- !remove else
      keep <- keep & !remove
    
    n <- sum(keep)
  }
  
  if(n==0) stop("No individuals left after keep/remove! Make sure the FID/IID are correct.")
  if(p==0) stop("No SNPs left after extract/exclude/chr! Make sure the SNP ids are correct.")
  
  if(!export) {
    return(list(keep=keep, extract=extract,
                N=N, P=P, n=n, p=p, bfile=bfile,
                bimfile=bimfile, famfile=famfile,
                bim=NULL, fam=NULL))
  } else {
    return(list(keep=keep, extract=extract,
                N=N, P=P, n=n, p=p, bfile=bfile,
                bimfile=bimfile, famfile=famfile,
                bim=bim, fam=fam))
  }
  #' @return A list containing:
  #' \item{keep}{Either NULL or a logical vector of which individuals to keep}
  #' \item{extract}{Either NULL or a logical vector of which SNPs to extract}
  #' \item{N}{Number of rows in the PLINK bfile}
  #' \item{P}{Number of columns in the PLINK bfile}
  #' \item{n}{Number of rows in the PLINK bfile after keep}
  #' \item{p}{Number of columns in the PLINK bfile after extract}
  #' \item{bfile}{File path to bfile stub}
  #' \item{bimfile}{File path to bim file}
  #' \item{famfile}{File path to fam file}
  #' \item{bim}{Either \code{NULL} or a data frame of bim data}
  #' \item{fam}{Either \code{NULL} or a data frame of fam data}
  
}

parseblocks <- function(vec) {
  
  #' @keywords internal
  if(is.factor(vec)) vec <- as.integer(vec)
  vec <- as.integer(factor(vec, levels=unique(vec)))
  blocks <- unique(vec)
  stopifnot(blocks == sort(blocks))
  stopifnot(min(blocks) == 1)
  stopifnot(max(blocks) == length(blocks))
  rle <- rle(vec)
  endvec <- cumsum(rle$lengths)
  startvec <- c(0, endvec[-length(endvec)])
  endvec <- endvec - 1 
  return(list("startvec"=startvec, "endvec"=endvec))
  
}




#' @title Obtains the number of column (SNPs) in a PLINK bfile
#' 
#' @param bfile Plink file stem
#' @return an integer with the number of columns
#' #@keywords internal
#' @export
ncol.bfile <- function(bfile) {
  bimfile <- paste0(bfile, ".bim")
  if(!file.exists(bimfile)) 
    stop(paste0("Cannot find ", bimfile)) 
  
  return(countLines(bimfile))
  # if(Sys.info()["sysname"] == "Windows") {
  # 	wc.output <- shell(paste("wc -l", bimfile), intern=T)
  # } else {
  # 	wc.output <- system(paste("wc -l", bimfile), intern=T)
  # }
  # return(as.numeric(strsplit(wc.output, split = "\\s+")[[1]][1]))
}

#' @title Obtains the number of individuals in a PLINK bfile
#' 
#' @param bfile Plink file stem
#' @return an integer with the number of rows
#' #@keywords internal
#' @export
nrow.bfile <- function(bfile) {
  famfile <- paste0(bfile, ".fam")
  if(!file.exists(famfile)) 
    stop(paste0("Cannot find ", famfile)) 
  return(countLines(famfile))
  # if(Sys.info()["sysname"] == "Windows") {
  # 	wc.output <- shell(paste("wc -l", famfile), intern=T)
  # } else {
  # 	wc.output <- system(paste("wc -l", famfile), intern=T)
  # }
  # return(as.numeric(strsplit(wc.output, split = "\\s+")[[1]][1]))
}

#' @title Internal function to parse extract
#'
#' @param keep A boolean vector of which position to keep
#' @keywords internal
#' 
selectregion <- function(keep) {
  
  pos.keep <- which(keep)
  end <- pos.keep[-1]
  start <- pos.keep[-length(pos.keep)]
  diff.keep <- end - start
  skip <- diff.keep > 1
  starts <- c(start[1], end[skip])
  rle <- rle(keep)
  # print(rle)
  lengths <- rle$lengths[rle$values]
  return(list(starts, lengths))
  
}



parse.pheno.covar <- function(pheno, covar, parsed, trace=0) {
  #' @keywords internal
  fam <- parsed[['fam']]
  keep <- parsed$keep
  # keep <- NULL
  pheno.df <- NULL
  
  update.keep <- function(old, new) {
    if(all(new)) {
      return(old)
    } else {
      if(is.null(old)) return(new) else {
        if(is.null(new)) return(old) else 
          return(old & new)
      }
    }
  }
  #### pheno ####
  if(!is.null(pheno) && is.character(pheno) && length(pheno) == 1) {
    if(file.exists(pheno)) pheno <- read.table2(pheno, header=T) else 
      stop(paste("Cannot find", pheno))
  }
  if(is.data.frame(pheno)) {
    if(ncol(pheno) != 3) {
      stop(paste("A pheno data.frame must have 3 columns exactly",
                 "with the first 2 with headers 'FID' and 'IID'"))
    }
    colnames <- colnames(pheno) 
    if(!all(colnames[1:2] == c("FID", "IID"))) {
      stop(paste("The first two columns of the pheno", 
                 "data.frame must have headers 'FID' and 'IID'"))
    }
    if(is.null(fam)) fam <- read.table2(parsed$famfile)
    rownames(fam) <- paste(fam$V1, fam$V2, sep="_")
    pheno.df <- pheno
    colnames(pheno.df)[3] <- "pheno"
    rownames(pheno) <- paste(pheno$FID, pheno$IID, sep="_")
    keep <- update.keep(keep, rownames(fam) %in% rownames(pheno))
    Pheno <- as.data.frame(pheno)[,3]
    names(Pheno) <- rownames(pheno)
  } else {
    if(!is.null(pheno)) {
      stopifnot(length(pheno) == parsed$n)
    } else {
      fam <- read.table2(parsed$famfile)
      if(is.null(parsed$keep)) pheno <- fam$V6 else 
        pheno <- fam$V6[parsed$keep]
    }
  }
  
  #### covar ####
  user.covar <- FALSE
  if(!is.null(covar) && is.character(covar) && length(covar) == 1) {
    if(file.exists(covar)) covar <- read.table2(covar, header=T) else 
      stop(paste("Cannot find", covar))
  }
  if(is.data.frame(covar) & all(colnames(covar)[1:2] == c("FID", "IID"))) {
    user.covar <- TRUE
    covar <- as.data.frame(covar)
    colnames <- colnames(covar) 
    if(is.null(fam)) fam <- read.table2(parsed$famfile)
    rownames(fam) <- paste(fam$V1, fam$V2, sep="_")
    rownames(covar) <- paste(covar$FID, covar$IID, sep="_")
    keep <- update.keep(keep, rownames(fam) %in% rownames(covar))
    Covar <- covar[,-(1:2), drop=FALSE]
  } else {
    if(!is.null(covar)) {
      if(is.vector(covar)) covar <- matrix(covar, ncol=1)
      if(is.matrix(covar)) covar <- as.data.frame(covar)
      Covar <- covar
    } 
  }
  
  #### updates ####
  parsed$keep <- update.keep(parsed$keep, keep)
  if(!is.null(parsed$keep)) parsed$n <- sum(parsed$keep)
  if(is.data.frame(pheno)) {
    if(!is.null(parsed$keep)) {
      names <- rownames(fam)[parsed$keep]
    } else {
      names <- rownames(fam)
    }
    pheno <- Pheno[names] 
    if(trace) {
      message(length(pheno), " out of ", length(Pheno), " samples kept in pheno.")
      # message(paste("Note that the order of best.pgs is the order given in the .fam file", 
      #               " rather than the pheno data.frame. Use v$best.pgs[v$order] to get", 
      #               " the pgs in the order of the phenotype."))
    }
    Order <- 1:length(pheno)
    names(Order) <- names
    pheno.df$order <- Order[names(Pheno)]
  } 
  
  if(user.covar) {
    if(!is.null(parsed$keep)) covar <- Covar[rownames(fam)[parsed$keep],,drop=F] else 
      covar <- Covar[rownames(fam),,drop=F]
    if(trace) message(nrow(covar), " out of ", nrow(Covar), " samples kept in covar.")
  } 
  
  if(length(pheno) == 0) {
    stop("No phenotype left. Perhaps the FID/IID do not match?")
  } else if(length(pheno) != parsed$n) {
    stop("The length of pheno does not match the number of samples.")
  }
  if(!is.null(covar) && nrow(covar) != parsed$n) {
    stop(paste("The dimension of the covar matrix does not match the number of samples used.", 
               "If your covariate is a data.frame with FID and IID, make sure they have headers."))
  }
  # if(sd(pheno, na.rm = TRUE) == 0) stop("There's no variation in phenotype")
  parsed$fam <- fam
  
  return(list(pheno=pheno, covar=covar, parsed=parsed, table=pheno.df))
  
}


#' @title Elastic net using summary statistics
#' @description Coordinate descent algorithm to solve: 
#' 0.5 x'X'Xx - x'b + lambda1 ||x||_1 + 0.5 lambda2 ||x||_2^2
#' Function to get elastic net solutions given X, a reference panel, and
#' b, regression coefficients
#' @keywords internal
myelnet <- function(lambda1, lambda2=0.0, gamma, fileName1, fileName2, b1, b2, 
                    #thr=1e-4, #lambda2 is s
                    parsed1, parsed2,
                    trace=0, maxiter=10000, 
                    blocks1=NULL, blocks2=NULL,
                    x=NULL) {
  #stopifnot(length(b1) == ncol(X1))
  #stopifnot(length(b2) == ncol(X2))
  #stopifnot(length(b1) == length(b2))
  #diag1 <- colSums(X1^2)
  #diag2 <- colSums(X2^2)
  
  N1 <- parsed1$N
  N2 <- parsed2$N
  stopifnot(parsed1$p==parsed2$p)
  P <- parsed1$P
  
  if(length(lambda2) > 1) {
    nlambda2 <- length(lambda2)
    for(i in 1:nlambda2) {
      result <- myelnet(lambda1, lambda2[i], gamma, fileName1, fileName2, b1, b2, 
                        #thr,
                        parsed1, parsed2,
                        trace, maxiter, x)
      result <- list(fit=result, lambda2=lambda2[i])
      if(i == 1) Result <- rep(result, nlambda2) else
        Result[i] <- result
      
    }
    return(Result)
  }
  
  #order <- order(lambda1, decreasing = T)
  #lambda1a <- lambda1[order]
  #conv <- lambda1a * NA
  len <- length(b1)
  #beta <- matrix(NA, len, length(lambda1))
  #pred1 <- matrix(NA, N1, length(lambda1))
  #pred2 <- matrix(NA, N2, length(lambda1))
  #loss <- rep(NA, length(lambda1))
  #fbeta <- loss
  
  if(is.null(x)) x <- b1 * 0.0 else { #initialize beta with 0s
    stopifnot(length(x) == len)
    x <- x + 0.0 # Making sure R creates a copy...
  }
  
  if(is.null(blocks1)) {
    Blocks1 <- list(startvec=0, endvec=len - 1)
    blocks1 <- rep(0, len)
  } else {
    Blocks1 <- parseblocks(blocks1)
    blocks1 <- blocks1 - 1 # for in C, pos starts from 0
    stopifnot(max(Blocks1$endvec)==len - 1) # leave it hear, the functions should be the same
  }
  
  if(is.null(blocks2)) {
    Blocks2 <- list(startvec=0, endvec=len - 1)
    blocks2 <- rep(0, len)
  } else {
    Blocks2 <- parseblocks(blocks2)
    blocks2 <- blocks2 - 1 # for in C, pos starts from 0
    stopifnot(max(Blocks2$endvec)==len - 1) # leave it hear, the functions should be the same
  }
  
  if(is.null(parsed1$extract)) {
    extract1 <- list(integer(0), integer(0))
    #col_skip_pos = integer(0); col_skip = integer(0)
  } else {
    #print(parsed1$extract)
    extract1 <- selectregion(!parsed1$extract)
    extract1[[1]] <- extract1[[1]] - 1
  }
  
  if(is.null(parsed1$keep)) {
    keepbytes1 <- integer(0)
    keepoffset1 <- integer(0)
  } else {
    pos <- which(parsed1$keep) - 1
    keepbytes1 <- floor(pos/4)
    keepoffset1 <- pos %% 4 * 2
  }
  
  if(is.null(parsed2$extract)) {
    extract2 <- list(integer(0), integer(0))
    #col_skip_pos = integer(0); col_skip = integer(0)
  } else {
    #print(parsed2$extract)
    extract2 <- selectregion(!parsed2$extract)
    extract2[[1]] <- extract2[[1]] - 1
  }
  
  if(is.null(parsed2$keep)) {
    keepbytes2 <- integer(0)
    keepoffset2 <- integer(0)
  } else {
    pos <- which(parsed2$keep) - 1
    keepbytes2 <- floor(pos/4)
    keepoffset2 <- pos %% 4 * 2
  }
  #X1 <- as.matrix(X1)
  #X2 <- as.matrix(X2)
  #yhat1 <- as.vector(X1 %*% x) #initialize PRS pred for population 1
  #yhat2 <- as.vector(X2 %*% x) #initialize PRS pred for population 2
  
  
  #  for(i in 1:length(lambda1a)) {
  #    if(trace > 0) cat("lambda1: ", lambda1a[i], "\n")
  #conv[i] <- myrepelnet(lambda1a[i], lambda2, gamma, diag1, diag2, X1, X2, b1, b2, 
  #thr,
  #                      x,yhat1, yhat2, trace-1,maxiter,
  #                      Blocks1$startvec, Blocks1$endvec, Blocks2$startvec, Blocks2$endvec, blocks1, blocks2)
  #if(conv[i] != 1) stop("Not converging...")
  results<- myrunElnet(lambda=lambda1, shrink=lambda2, gamma=gamma, 
                       fileName1=fileName1, fileName2=fileName2,
                       r1=b1, r2=b2, 
                       N1=N1, N2=N2, P=P, 
                       col_skip_pos1 = extract1[[1]], col_skip1 = extract1[[2]], 
                       keepbytes1 = keepbytes1, keepoffset1 = keepoffset1, 
                       col_skip_pos2 = extract2[[1]], col_skip2 = extract2[[2]], 
                       keepbytes2 = keepbytes2, keepoffset2 = keepoffset2, 
                       #thr, 
                       x, trace, maxiter, 
                       startvec1 = Blocks1$startvec, endvec1 = Blocks1$endvec,
                       startvec2 = Blocks2$startvec, endvec2 = Blocks2$endvec,
                       blocks1, blocks2)
  
  #    beta[,i] <- x
  #    pred1[,i] <- yhat1
  #    pred2[,i] <- yhat2
  #    loss[i] <- gamma * (sum(yhat1^2) - 2* sum(b1 * x)) + (1-gamma) * (sum(yhat2^2) - 2* sum(b2 * x)) 
  #    fbeta[i] <- loss[i] + 2* sum(abs(x))*lambda1a[i] + sum(x^2)*lambda2
  #  }
  
  
  #conv[order] <- conv
  #beta[,order] <- beta
  #pred1[,order] <- pred1
  #pred2[,order] <- pred2
  #loss[order] <- loss
  #fbeta[order] <- fbeta
  
  return(list(lambda1=lambda1, lambda2=lambda2, gamma=gamma, beta=results$beta, conv=results$conv, 
              pred1=results$pred1, pred2=results$pred2, loss=results$loss, fbeta=results$fbeta,
              sd1=results$sd1, sd2=results$sd2))
  
}


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
    }
  }
  
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
                   fbeta=el$fbeta,
                   sd1=el$sd1,
                   sd2=el$sd2,
                   shrink=shrink,
                   nparams=nparams, 
                   gamma=gamma
                   #,ridge=Ridge
  )
  
  class(toreturn) <- "lassosum"
  return(toreturn)
  
  #' @return A list with the following
  #' \item{lambda}{same as the lambda input}
  #' \item{beta}{A matrix of estimated coefficients}
  #' \item{conv}{A vector of convergence indicators. 1 means converged. 0 not converged.}
  #' \item{pred}{\eqn{=(1-s)X\beta}}
  #' \item{loss}{\eqn{=(1-s)\beta'X'X\beta/n - 2\beta'r}}
  #' \item{fbeta}{\eqn{=\beta'R\beta - 2\beta'r + 2\lambda||\beta||_1}}
  #' \item{sd}{The standard deviation of the reference panel SNPs}
  #' \item{shrink}{same as input}
  #' \item{nparams}{Number of non-zero coefficients}
  #' \item{gamma}{A parameter balancing two populations}
  #' \item{ridge}{ridge regression estimates}
  
}


splitgenome <- function(CHR, POS, ref.CHR, ref.breaks, details=T, right=TRUE) {
  #' Function to split a set of SNPs by their position using a reference
  #' Break points should be given for the reference by chromosome
  #' @keywords internal
  
  CHR <- as.character(CHR)
  ref.CHR <- as.character(ref.CHR)
  POS <- as.integer(POS)
  ref.breaks <- as.integer(ref.breaks)
  
  stopifnot(all(!is.na(POS)) && all(!is.na(ref.breaks)))
  stopifnot(all(POS >= 1) && all(ref.breaks >= 1))
  stopifnot(length(CHR) == length(POS))
  stopifnot(length(ref.CHR) == length(ref.breaks))
  
  
  chr <- (unique(CHR))
  chr.ref <- (unique(ref.CHR))
  included.chr <- chr %in% chr.ref
  # if(all(!included.chr)) stop("Cannot match the chromosomes. Make sure the notations are the same. e.g. 'chr1' vs 1 or chrX vs chr23.")
  if(!all(included.chr)) stop("Some chromosomes were not defined in the reference. Make sure the notations are the same. e.g. 'chr1' vs 1 or chrX vs chr23.")
  
  levels <- character(0)
  results <- character(length(POS))
  Details <- data.frame()
  for(C in chr.ref) {
    breaks <- sort(unique(ref.breaks[ref.CHR == C]))
    if(breaks[1] > 1) breaks <- c(1, breaks)
    if(breaks[length(breaks)] < Inf) breaks <- c(breaks, Inf)
    cut <- cut(POS[CHR == C],include.lowest = T,breaks = breaks, right=right)
    levels <- c(levels, paste0(C, "_", levels(cut)))
    cut <- paste0(C, "_", as.character(cut))
    results[CHR == C] <- cut
    if(details) {
      df <- data.frame(chr=C, start=breaks[-length(breaks)], end=breaks[-1])
      Details <- rbind(Details, df)
    }
  }
  results <- factor(results, levels = levels)
  
  if(details) {
    Details$counts <- as.integer(table(results))
    attr(results, "details") <- Details
  }
  
  return(results)
  
}


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


#' @title Merge lassosum results 
#' @description e.g. when calculated over different blocks/chromosomes
#' @method merge lassosum
#' @export
#' 
merge.mylassosum <- function(ll) {
  
  #  ll <- list(...)
  #  stopifnot(all(sapply(ll, "class") == "lassosum"))
  stopifnot(all(sapply(ll, function(x) all(x$lambda == ll[[1]]$lambda))))
  shrink <- sapply(ll, function(x) x$shrink)
  stopifnot(all(shrink == shrink[1]))
  gamma <- sapply(ll, function(x) x$gamma)
  stopifnot(all(gamma == gamma[1]))
  
  Cumsum <- function(...) {
    mat <- do.call("rbind", list(...))
    if(ncol(mat) > 0) return(as.vector(colSums(mat))) else
      return(numeric(0))
  }
  results <- list()
  results$lambda <- ll[[1]]$lambda
  results$beta <- do.call("rbind", lapply(ll, function(x) x$beta))
  results$conv <- do.call("pmin", lapply(ll, function(x) x$conv))
  pred1 <- do.call("Cumsum", lapply(ll, function(x) as.vector(x$pred1)))
  pred2 <- do.call("Cumsum", lapply(ll, function(x) as.vector(x$pred2)))
  results$pred1 <- matrix(pred1, ncol=length(results$lambda), nrow=nrow(ll[[1]]$pred1))
  results$pred2 <- matrix(pred2, ncol=length(results$lambda), nrow=nrow(ll[[1]]$pred2))
  results$loss <- do.call("Cumsum", lapply(ll, function(x) x$loss))
  results$fbeta <- do.call("Cumsum", lapply(ll, function(x) x$fbeta))
  results$sd1 <- do.call("c", lapply(ll, function(x) x$sd1))
  results$sd2 <- do.call("c", lapply(ll, function(x) x$sd2))
  results$shrink <- ll[[1]]$shrink
  results$nparams <- do.call("Cumsum", lapply(ll, function(x) x$nparams))
  results$gamma <- ll[[1]]$gamma
  class(results) <- "lassosum"
  return(results)
}