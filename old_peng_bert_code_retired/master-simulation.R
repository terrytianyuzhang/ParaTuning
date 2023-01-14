rm(list=ls()); gc()
options(stringsAsFactors = F)
# 600-609 h2 = 0.80 in both populations
# 700-709 h2=0.80 in CEU and 0.60 in YRI
# 800-809 H2 approx 0.8, betas the same in both populations
sims=800:809
for(i.sim in sims){
  seed=sample(1e6,1)
  print(i.sim)
  # set up the parameters
  source("simulation-parameters.varg.ld-blocks-ref-alt.R")
}
  
#simulate the populations
for(i.sim in sims){
  print(i.sim)
  source("simulate_population_multinode.varg.ld-blocks-ref-alt.R")
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

# do pruning based on gwas results in the training populations
for(i.sim in sims){
  print(i.sim)
  source("PGS-PT-pruning.R")
}  

# testing part
for(i.sim in sims){
  print(i.sim)
  source("PGS-PT-testing.R")
}  

# weighted-PT
for(i.sim in sims){
  print(i.sim)
  source("PGS-PT-weighted.R")
}  

# lassosum
for(i.sim in sims){
  print(i.sim)
  source("PGS-lassosum.R")
} 
  
# weighted lassosum
for(i.sim in sims){
  print(i.sim)
  source("PGS-lassosum-weighted.R")
} 

# multi lassosum, this runs all sims at the same time
source("BERT_ourx3_tuning-multi-sim.R")

# multi lassosum for testing
for(i.sim in sims){
  print(i.sim)
  source("PGS-lassosum-multi.R")
} 




