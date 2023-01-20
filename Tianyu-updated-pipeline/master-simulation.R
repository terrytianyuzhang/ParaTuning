rm(list=ls()); gc()
options(stringsAsFactors = F)
# 600-609 h2 = 0.80 in both populations
# 700-709 h2=0.80 in CEU and 0.60 in YRI
# 800-809 H2 approx 0.8, betas the same in both populations

plink <- "/usr/local/bin/plink"
plink2 <- "/usr/local/bin/plink2"
# sims=800:809
sims <- 800:809 ##only consider one replicate

# for(i.sim in sims){
#   seed=sample(1e6,1)
#   print(i.sim)
#   # set up the parameters
#   source("simulation-parameters.varg.ld-blocks-ref-alt.R")
# }
# 
# # simulate the populations
# # I changed the number of nodes
# for(i.sim in sims){
#   print(i.sim)
#   source("simulate_population_multinode.varg.ld-blocks-ref-alt.R")
#   print(paste0("finished simulating population for ", i.sim))
# }
# 
# 
# # translate TUNE and TST into a bfile format and split by chromosome
# for(i.sim in sims){
#   print(i.sim)
#   source("split-pfile-by-chr-into-bfile.R")
# }
# 
# 
# 
# #run gwas on the training populations
# for(i.sim in sims){
#   print(i.sim)
#   source("run-gwas.R")
# }

# 
# #the PT section is skipped for now because PLINK is too old
# 
# # # do pruning based on gwas results in the training populations
# for(i.sim in sims){
#   print(i.sim)
#   source("PGS-PT-pruning.R")
# }
# print('finished PT pruning')
# 
# 
# # # testing part
# for(i.sim in sims){
#   print(i.sim)
#   source("PGS-PT-testing.R")
# }
# 
# print('finished PT testing')
# 
# 
# # weighted-PT
# for(i.sim in sims){
#   print(i.sim)
#   source("PGS-PT-weighted.R")
# }
# 
# print('finished PT-weighted')
# 
# 
# # # lassosum
# for(i.sim in sims){
#   print(i.sim)
#   source("PGS-lassosum-bestpara.R")
# }
# #
# # weighted lassosum
# for(i.sim in sims){
#   print(i.sim)
#   source("PGS-lassosum-weighted-bestpara.R")
# }
# 
# # fit joint lassosum
for(i.sim in sims){
  print(i.sim)
  source("PGS-JointLassosum.R")
}

# # get testing AUC of Joint Lassosum
for(i.sim in sims){
  print(i.sim)
  source("PGS-JointLassosum-testing.R")
}
# 
# library(devtools)
# install_github("terrytianyuzhang/TLPRS")
# 
# for(i.sim in sims){
#   print(i.sim)
#   source("PGS-TL.R")
# }

# if(1){
# # multi lassosum, this runs all sims at the same time
# source("BERT_ourx3_tuning-multi-sim.R")
# 
# # multi lassosum for testing
# for(i.sim in sims){
#   print(i.sim)
#   source("PGS-lassosum-multi.R")
# }
# 
# }
# 

