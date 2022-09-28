#####this piece of code is for studying how well we can approximte a design matrix
######using truncated SVD and approximation of the U matrix
library(MASS)
p <- 100
n <- 50
transform_m <- diag(1,p)
for(i in 1:(p-1)){
  transform_m[i,i+1] <- 1
}

# transform_m <- matrix(c(1,1,0,0,
#                         0,1,1,0,
#                         0,0,1,1,
#                         0,0,0,1), ncol = 4)
sigma <- crossprod(transform_m)
x <- mvrnorm(n, mu = rep(0, p), Sigma = sigma)

k <- min(n,p)  # number of eigenvalues to calculate, I should only need the first n
res <- eigs_sym(crossprod(x), k, which = "LM")  # "LM" is the default
D <- sqrt(pmax(res$values,0)) #square root of eigenvalues are the singualr values
V <- res$vectors

tau <- 1
epsilon_tilde <- rnorm(n, mean = 0, sd = tau)
app <- V %*% (D * epsilon_tilde[1:k])

true <- crossprod(x,epsilon_tilde)
plot(app, true)

uu <- svd(x)$u
vv <- svd(x)$v
# install.packages('RSpectra')
library(RSpectra)
library(lassosum)
library(data.table)
library(snpStats)
setwd('/home/tianyuz3/PRS/my_code/')

chr_index <- 'chr20'
ld.anc <- 'AFR' ###ancestor
if(ld.anc == 'AFR'){
  anc_initial <- 'YRI'
}

##############LOAD THE LD BLOCK BOUNDARY#####

##load the boundary information
LD <- read.table(system.file(paste0("data/Berisa.", 
                                      paste0(ld.anc,".hg38"), ".bed"), 
                               package="lassosum"), header=T)
LD <- data.table(LD)

##############LOAD THE DESIGN MATRICES#####

# Get the file paths for the reference data: -----
fam_path <- paste0("/raid6/Ron/prs/data/bert_sample/GWAS-Populations-SimulationInput/", anc_initial, "_reference_LDblocks/CHR/", anc_initial, "-",chr_index, ".fam")
bim_path <- paste0("/raid6/Ron/prs/data/bert_sample/GWAS-Populations-SimulationInput/", anc_initial, "_reference_LDblocks/CHR/", anc_initial, "-",chr_index, ".bim")
bed_path <- paste0("/raid6/Ron/prs/data/bert_sample/GWAS-Populations-SimulationInput/", anc_initial, "_reference_LDblocks/CHR/", anc_initial, "-",chr_index, ".bed")

# Read in the PLINK data -----
ref_snps_plink <- read.plink(bed_path, 
                             bim_path,
                             fam_path)


# Obtain the SnpMatrix object (genotypes) table from ref_snps_plink list
ref_genotypes <- ref_snps_plink$genotypes
print(ref_genotypes)
# A SnpMatrix with  503 rows and  22665064 columns
# Row names:  HG00096 ... NA20832 
# Col names:  rs537182016 ... rs781880 

####get a numeric version of it
numeric_ref_genotypes <- as(ref_genotypes, "numeric")

#Obtain the SNP information from ref_snps_plink list
ref_genobim <- data.table(ref_snps_plink$map)
colnames(ref_genobim) <- c("chr", "SNP", "gen.dist", "position", "A1", "A2")

######perform SVD block by block#########
SVD_big_list <- list()
LD_num_in_chr <- NROW(LD[chr == chr_index, ])
for(LD_index in 1:LD_num_in_chr){
  # LD_index <- 1
  print(LD_index)
  LD[chr == chr_index, ][3,]
  ###     chr start   stop
  #  1: chr20 79838 853837
  LD_lower <- as.numeric(LD[chr == chr_index, ][LD_index,2])
  LD_upper <- as.numeric(LD[chr == chr_index, ][LD_index,3])
  
  # print(head(ref_genobim))
  SNP_inthis_LD <- which(ref_genobim$position >= LD_lower & ref_genobim$position < LD_upper)
  
  
  #####slice the relevant SNPs
  X <- numeric_ref_genotypes[,SNP_inthis_LD]
  ###subtract column sums
  X <- X - matrix(rep(apply(X, MARGIN = 2, FUN = mean), 
                      NROW(X)), 
                  nrow = NROW(X), byrow = TRUE)
  C <- crossprod(X)/NROW(X)
  
  ######CALCULATE SVD###########
  k <- 200  # number of eigenvalues to calculate, I should only need the first n
  res <- eigs_sym(NROW(X)*C, k, which = "LM")  # "LM" is the default
  D <- sqrt(pmax(res$values,0)) #square root of eigenvalues are the singualr values
  V <- res$vectors
  
  temp_list <- list(chr_index = chr_index,
                    ancestor = ld.anc,
                    LD_lower = LD_lower, 
                    LD_upper = LD_upper,
                    C = C,
                    D = D,
                    V = V)
  SVD_big_list[[length(SVD_big_list)+1]] <- temp_list
}
save(SVD_big_list, file = paste0("/raid6/Tianyu/PRS/SVDdata/", anc_initial, "reference_LDblocks_", chr_index, "_SVD.RData"))
# 
# ######LOAD THE ORIGINAL LD MATRICES######
# setwd('/Users/tianyu/Documents/ParameterTuning/')
# # set.seed(123)
# # n = 100  # dimension of X
# # # Some random data
# # X = matrix(rnorm(n*n/2), nrow = n/2, ncol = n)
# X <- matrix(c(1,0,0.5,1,1,0.5), ncol = 3)
# beta0 <- c(1,1,1)
# Y <- X %*% beta0 + rnorm(2, 0, 0.05)
# # Make it symmetric
# C <- crossprod(X)/2
# R <- crossprod(X,Y)/2
# # # Show its largest 5 eigenvalues and associated eigenvectors
# # head(eigen(A)$values, 5)
# 
# ######CALCULATE SVD###########
# k <- 3  # number of eigenvalues to calculate, I should only need the first n
# res <- eigs_sym(crossprod(X), k, which = "LM")  # "LM" is the default
# D <- sqrt(pmax(res$values,0)) #square root of eigenvalues are the singualr values
# V <- res$vectors
# #### (XTX)V = VD
# ######SAVE THE FILES###########
# save(D, file = 'SingularValueD.RData')
# save(V, file = 'SingularVectorV.RData')
# save(C, file = 'CovairanceMatrixC.RData')
# save(R, file = 'EmpiricalCorrR.RData')
