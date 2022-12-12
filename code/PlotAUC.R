setwd('/Users/tianyu/Documents/ParaTuning/data/Simulation-AUC')


results <- data.frame()


for(i.sim in 800:801){
#######PT###########
tune.sets <- c('CEU.TUNE', 'YRI.TUNE')
for(tune.set in tune.sets){
  PTraw <- get(load(paste0("No", i.sim, tune.set, "weightedPT-AUC.RData")))
  temp.results <- data.frame(method = 'PT-CEU',
                             i.sim = i.sim,
                             test.population = gsub('.TUNE', '', tune.set),
                             AUC = PTraw[PTraw$wght == 1,]$AUC)
  results <- rbind(results, temp.results)
  
  temp.results <- data.frame(method = 'PT-weighted',
                             i.sim = i.sim,
                             test.population = gsub('.TUNE', '', tune.set),
                             AUC = max(PTraw$AUC))
  results <- rbind(results, temp.results)
}

#######Lassosum###########
test.sets <- c('CEU.TST', 'YRI.TST')
test.set <- test.sets[1]
for(test.set in test.sets){
  lassosumraw <- get(load(paste0("No", i.sim, test.set, "weightedLassosum-AUC.RData")))
  
  temp.results <- data.frame(method = 'Lsm-CEU',
                             i.sim = i.sim,
                             test.population = gsub('.TST', '', test.set),
                             AUC = lassosumraw[lassosumraw$wght == 1,]$AUC)
  results <- rbind(results, temp.results)
  
  temp.results <- data.frame(method = 'Lsm-weighted',
                             i.sim = i.sim,
                             test.population = gsub('.TST', '', test.set),
                             AUC = max(lassosumraw$AUC))
  results <- rbind(results, temp.results)
  
}

###Joint Lassosum###
GAMMA <- c(0.2,0.5,0.8)
test.sets <- c('CEU.TST', 'YRI.TST')
for(test.set in test.sets){
  all.AUC <- c()
  for(gamma in GAMMA){
    JLraw <- get(load(paste0("No", i.sim, sprintf("-gamma-%.2f",gamma),"-", test.set, "JointLassosum-AUC.RData")))
    all.AUC <- c(all.AUC, JLraw)
  }
  temp.results <- data.frame(method = 'JL',
                             i.sim = i.sim,
                             test.population = gsub('.TST', '', test.set),
                             AUC = max(all.AUC))
  results <- rbind(results, temp.results)
}

}#i.sim

##########violin plots#####
library(ggplot2)
results <- rbind(results,results,results,results)
results$AUC <- results$AUC + rnorm(NROW(results), 0, 0.01)
ggplot(results[results$test.population == 'YRI',]) +
  geom_violin(aes(x=method, y=AUC, fill = as.factor(method)),
              trim=FALSE)+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9",
                             "#F6B419", "#56B419"))+
  theme_bw()
  
  # geom_point(aes(x=method, y=AUC))
  
ggplot(results[results$test.population == 'YRI',]) +
  geom_violin(aes(x=method, y=AUC))
