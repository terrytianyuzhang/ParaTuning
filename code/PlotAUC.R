library(data.table)


setwd('/Users/tianyu/Documents/ParaTuning/data/')
rawSimulationAUC <- get(load('simulations-Results-6xx-7xx-8xx01-13-2023.RData'))
length(rawSimulationAUC)
rawSimulationAUC[[1]]

clean_simulation_AUC <- data.table()
for(setting_index in 1:length(rawSimulationAUC)){
  
  pnt_result <- rawSimulationAUC[[setting_index]]$pnt
  pnt_result_formated <- pnt_result[,c('TRN','TST','auc')]
  pnt_result_formated$setting <- names(rawSimulationAUC)[setting_index]
  pnt_result_formated$method <- 'PnT'
  
  clean_simulation_AUC <- rbind(clean_simulation_AUC, 
                                pnt_result_formated)
  
  wpnt_result <- rawSimulationAUC[[setting_index]]$wpnt
  wpnt_result_formated <- wpnt_result[,c('TST','auc')]
  wpnt_result_formated$TRN <- 'Both'
  wpnt_result_formated$setting <- names(rawSimulationAUC)[setting_index]
  wpnt_result_formated$method <- 'PnT-wt'
  
  clean_simulation_AUC <- rbind(clean_simulation_AUC, 
                                wpnt_result_formated)
  
  ls_result <- rawSimulationAUC[[setting_index]]$ls
  ls_result_formated <- ls_result[,c('TRN','TST','auc')]
  ls_result_formated$setting <- names(rawSimulationAUC)[setting_index]
  ls_result_formated$method <- 'LS'
  
  clean_simulation_AUC <- rbind(clean_simulation_AUC, 
                                ls_result_formated)
  
  wls_result <- rawSimulationAUC[[setting_index]]$wls
  wls_result_formated <- wls_result[,c('TST','auc')]
  wls_result_formated$TRN <- 'Both'
  wls_result_formated$setting <- names(rawSimulationAUC)[setting_index]
  wls_result_formated$method <- 'LS-wt'
  
  clean_simulation_AUC <- rbind(clean_simulation_AUC, 
                                wls_result_formated)
  
  jls_result <- rawSimulationAUC[[setting_index]]$jls
  jls_result_formated <- jls_result[,c('TST','auc')]
  jls_result_formated$TRN <- 'Both'
  jls_result_formated$setting <- names(rawSimulationAUC)[setting_index]
  jls_result_formated$method <- 'JLS'
  
  clean_simulation_AUC <- rbind(clean_simulation_AUC, 
                                jls_result_formated)
  
  tl_result <- rawSimulationAUC[[setting_index]]$tl
  tl_result_formated <- tl_result[,c('TST','auc')]
  tl_result_formated$TRN <- 'Both'
  tl_result_formated$setting <- names(rawSimulationAUC)[setting_index]
  tl_result_formated$method <- 'TL'
  
  clean_simulation_AUC <- rbind(clean_simulation_AUC, 
                                tl_result_formated)
  # clean_simulation_AUC <- rbind(clean_simulation_AUC,
  #                               data.frame(setting = names(rawSimulationAUC)[setting_index],
  #                                          TrainData = )
}

clean_simulation_AUC[, method := factor(method,
                                        levels = c("JLS", "TL", "LS-wt", "PnT-wt", "LS", "PnT"))]
clean_simulation_AUC[, method_TRN := factor(paste0(method,'(',TRN,')'),
                     levels = c("JLS(Both)", "TL(Both)", "LS-wt(Both)", "PnT-wt(Both)", 
                                "LS(CEU)", "PnT(CEU)", "LS(YRI)", "PnT(YRI)"))]

##########violin plots#####
library(ggplot2)
library(ggpubr)
mainYRI <- ggplot(clean_simulation_AUC[TST == 'YRI',]) +
  geom_violin(aes(x=method_TRN, y=auc, fill = method_TRN),
              trim=T, alpha = 0.4)+
  geom_point(aes(x=method_TRN, y=auc,
                 color = method_TRN),  alpha = 0.8,
             position = position_jitter(width = 0.05))+
  scale_fill_manual(values=c("#9b5fe0", "#16a4d8", "#60dbe8",
                             "#8bd346", "#efdf48","#f9a52c",
                             "#d64e12","#B6BeB2"))+
  scale_color_manual(values=c("#9b5fe0", "#16a4d8", "#60dbe8",
  "#8bd346", "#efdf48","#f9a52c","#d64e12","#B6BeB2"))+
  ylim(c(0.5,0.9))+
  # xlab("Method")+
  xlab(NULL)+
  ylab("YRI tst. AUC")+
  theme_bw()
  # theme(axis.text.x=element_blank())

mainCEU <- ggplot(clean_simulation_AUC[TST == 'CEU',]) +
  geom_violin(aes(x=method_TRN, y=auc, fill = method_TRN),
              trim=T, alpha = 0.4)+
  geom_point(aes(x=method_TRN, y=auc,
                 color = method_TRN),  alpha = 0.8,
             position = position_jitter(width = 0.05))+
  scale_fill_manual(values=c("#9b5fe0",
                             # "#16a4d8",
                             "#60dbe8", "#8bd346", "#efdf48","#f9a52c","#d64e12","#B6BeB2"))+
  scale_color_manual(values=c("#9b5fe0",
                              # "#16a4d8", 
                              "#60dbe8", "#8bd346", "#efdf48","#f9a52c","#d64e12","#B6BeB2"))+
  ylim(c(0.5,0.9))+
  xlab("Method")+
  ylab("CEU tst. AUC")+
  theme_bw()



# ggsave(filename = "Sim8xxTestingAUC.pdf",
#        width = 18,
#        height = 5,
#        units = "cm")
mainCombined<- ggarrange(
  mainYRI,mainCEU,  
  nrow = 2,
  # labels = c("A","B"),
  common.legend = TRUE,legend = 'none'
) 
mainCombined
ggsave(filename = "Sim8xx-testing-AUC.pdf",
       device = "pdf",
       width = 180, height = 90, units = "mm")

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
