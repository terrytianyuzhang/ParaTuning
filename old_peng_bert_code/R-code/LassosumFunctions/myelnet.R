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
    print('loss in myelnet.R')
    print(results$loss)
    print(results$trainerror1)
    print(results$trainerror2)
    
  return(list(lambda1=lambda1, lambda2=lambda2, gamma=gamma, beta=results$beta, conv=results$conv, 
              pred1=results$pred1, pred2=results$pred2, loss=results$loss, fbeta=results$fbeta,
              trainerror1 = results$trainerror1, trainerror2 = results$trainerror2, 
              sd1=results$sd1, sd2=results$sd2))
  
}