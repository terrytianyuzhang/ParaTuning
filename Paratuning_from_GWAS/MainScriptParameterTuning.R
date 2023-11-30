rm(list=ls()); gc()
options(stringsAsFactors = F)
# 600-609 h2 = 0.80 in both populations
# 700-709 h2=0.80 in CEU and 0.60 in YRI
# 800-809 H2 approx 0.8, betas the same in both populations
# sims=800:809
sims <- 600:609 ##only consider one replicate

# fit joint lassosum with real data
# for(i.sim in sims){
#   print(i.sim)
#   source("fit_JLS.R")
# }

# for(i.sim in sims){
#   print(i.sim)
#   source("split-pfile-by-chr-into-bfile.R")
# }

# get the testing error of models trained with different hyper-parameters
# for(i.sim in sims){
#   print(i.sim)
#   source("test_JLS.R")
# }
# 
# # generate synthetic data for parameter tuning
for(i.sim in sims){
  print(i.sim)
  source("GenerateSyntheticData.R")
}
#
# # fit joint lassosum with synthetic data
for(i.sim in sims){
  print(i.sim)
  source("PGS-JL-SyntheticDataTraining.R")
}
#
# get testing AUC of Joint Lassosum, trained with synthetic data
for(i.sim in sims){
  print(i.sim)
  source("PGS-JL-SyntheticDataTesting.R")
}
