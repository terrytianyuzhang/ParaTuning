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

# simulate the populations
# I changed the number of nodes
for(i.sim in sims){
  print(i.sim)
  source("simulate_population_multinode.varg.ld-blocks-ref-alt.R")
  print(paste0("finished simulating population for ", i.sim))
}


# translate TUNE and TST into a bfile format and split by chromosome
for(i.sim in sims){
  print(i.sim)
  source("split-pfile-by-chr-into-bfile.R")
}



#run gwas on the training populations
for(i.sim in sims){
  print(i.sim)
  source("run-gwas.R")
}

# # fit joint lassosum with real data
# for(i.sim in sims){
#   print(i.sim)
#   source("PGS-JointLassosum.R")
# }
# 
# # # get testing AUC of Joint Lassosum
# for(i.sim in sims){
#   print(i.sim)
#   source("PGS-JointLassosum-testing.R")
# }
# 
# # generate synthetic data for parameter tuning
# for(i.sim in sims){
#   print(i.sim)
#   source("GenerateSyntheticData.R")
# }
# 
# # fit joint lassosum with synthetic data
# for(i.sim in sims){
#   print(i.sim)
#   source("PGS-JL-SyntheticDataTraining.R")
# }
# 
# # get testing AUC of Joint Lassosum, trained with synthetic data
# for(i.sim in sims){
#   print(i.sim)
#   source("PGS-JL-SyntheticDataTesting.R")
# }
