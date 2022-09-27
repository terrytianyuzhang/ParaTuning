#####there are 129,897 SNPs in chr20
#####dim(numeric_ref_genotypes)
#####[1]   5000 129897

#####create a fake \hat \beta_0
set.seed(1)
chr20snpN <- 129897
beta0 <- rep(0, chr20snpN)
non_zero_snp <- sample(1:chr20snpN, 100)
beta0[non_zero_snp] <- 1

#####create a fake tauj
tau <- 1

#####sample size is 5000

#####LOAD R C V D MATRICES######
D <- get(load(file = 'SingularValueD.RData'))
V <- get(load(file = 'SingularVectorV.RData'))
C <- get(load(file = 'CovairanceMatrixC.RData'))
R <- get(load(file = 'EmpiricalCorrR.RData'))
B <- 100
n <- 2 #sample size
beta_pre <- c(1,0.9,1.1) #preliminary beta
####VARIANCE OF THE NOISE#####
tau <- sqrt(n^(-1)*(beta_pre[1] + crossprod(beta_pre, crossprod(C, beta_pre))))
####FOR EACH BOOSTRAP REPLICATION#####
cov_b <- 0
b <- 1

  ####GENERATE NOISE EPSILON######
  epsilon_tilde <- rnorm(n, mean = 0, sd = tau)
  ####CALCULATE X^T EPSILON###
  Xtepsilon <- V[,1:n] %*% (D[1:n] * epsilon_tilde)
  ####CALCULATE BOOTSTRAP SUMMARY STAT####
  R_b <- n*C %*% beta_pre + Xtepsilon
  ####CALCULATE BOOTSTRAP BETA#######
  #####xxxxxxxxxxxxx
  beta <- beta_pre
  cov_b_each <- crossprod(beta, Xtepsilon)
  cov_b <- cov_b + cov_b_each
####CALCULATE THE BOOTSTRAP COV####
cov_b <- cov_b/B
####SAVE FILES####
  
###return cov_b
