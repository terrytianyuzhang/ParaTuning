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