library(data.table)
library(ggplot2)

normalize_PRS <- function(x){
  x <- (x - mean(x))/ (sd(x))
  return(x)
}

setwd('/Users/tianyu/Documents/ParaTuning/data/Simulation-PRS/')
####
rawPRS_YRI <- get(load('YRI-PRS-6xx-7xx-8xx01-11-2023.RData'))
rawPRS_CEU <- get(load('CEU-PRS-6xx-7xx-8xx01-11-2023.RData'))

#####
simulation_index_sets <- list(600:609, 700:709, 800:809)
#####
for(simulation_setting in 1:length(simulation_index_sets)){
  simulation_index_set <- simulation_index_sets[[simulation_setting]]
  #####
  FDR_all_simulations <- data.table() ##all the replicates for one simulation setting
  #####
  for(simulation_index in simulation_index_set){
    
    simulaton_name <- paste0('Sim-', simulation_index)
    
    #####
    rawPRS_YRI_oneSimulation <- data.table(rawPRS_YRI[[simulaton_name]])
    #####
    rawPRS_YRI_oneSimulation[, population := 'YRI']
    rawPRS_YRI_oneSimulation[, jls := normalize_PRS(jls)]
    rawPRS_YRI_oneSimulation[, wls := normalize_PRS(wls)]
    rawPRS_YRI_oneSimulation[, ls.cis := normalize_PRS(ls.cis)]
    rawPRS_YRI_oneSimulation[, ls.trans := normalize_PRS(ls.trans)]
    rawPRS_YRI_oneSimulation[, ls.CEU := ls.trans]
    rawPRS_YRI_oneSimulation[, ls.YRI := ls.cis]
    rawPRS_YRI_oneSimulation[, PHENO1 := factor(PHENO1)]
    #####
    
    rawPRS_CEU_oneSimulation <- data.table(rawPRS_CEU[[simulaton_name]])
    #####
    rawPRS_CEU_oneSimulation[, population := 'CEU']
    rawPRS_CEU_oneSimulation[, jls := normalize_PRS(jls)]
    rawPRS_CEU_oneSimulation[, wls := normalize_PRS(wls)]
    rawPRS_CEU_oneSimulation[, ls.cis := normalize_PRS(ls.cis)]
    rawPRS_CEU_oneSimulation[, ls.trans := normalize_PRS(ls.trans)]
    rawPRS_CEU_oneSimulation[, ls.CEU := ls.cis]
    rawPRS_CEU_oneSimulation[, ls.YRI := ls.trans]
    rawPRS_CEU_oneSimulation[, PHENO1 := factor(PHENO1)]
    #####
    PRS_both_population <- rbind(rawPRS_YRI_oneSimulation,
                                 rawPRS_CEU_oneSimulation)
    #####
    PRS_both_population[, is.case := as.numeric((PHENO1 == 2))]
    
    
    #####
    #####
    ##NEXT SECTION PLOTS THE FDR AGAINST PRS THRESHOLD
    #####
    #####
    
    FDR_all_method <- data.table()
    #####
    FDR_by_threshold <- data.table()
    
    for(threshold in seq(-4, 4, length = 50)){
      temp_table <- PRS_both_population[jls >= threshold, 
                             .(total_positive = .N,
                               true_positive = sum(is.case),
                               threshold = threshold),
                             by = population]
      FDR_by_threshold <- rbind(FDR_by_threshold,
                                temp_table)
    }
    
    FDR_by_threshold[, FDR:= 1 - true_positive/total_positive]
    FDR_by_threshold[, method:= 'JLssum']
    
    #####
    
    FDR_all_method <- rbind(FDR_all_method, FDR_by_threshold)
    
    ###########
    FDR_by_threshold <- data.table()
    
    for(threshold in seq(-4, 4, length = 50)){
      
      temp_table <- PRS_both_population[ls.CEU >= threshold, 
                                        .(total_positive = .N,
                                          true_positive = sum(is.case),
                                          threshold = threshold),
                                        by = population]
      FDR_by_threshold <- rbind(FDR_by_threshold,
                                temp_table)
    }
    
    FDR_by_threshold[, FDR:= 1-true_positive/total_positive]
    FDR_by_threshold[, method:= 'LS-CEU-train']
    
    #####
    
    FDR_all_method <- rbind(FDR_all_method, FDR_by_threshold)
    
    ###########
    FDR_by_threshold <- data.table()
  
    for(threshold in seq(-4, 4, length = 50)){
  
      temp_table <- PRS_both_population[ls.YRI >= threshold,
                                        .(total_positive = .N,
                                          true_positive = sum(is.case),
                                          threshold = threshold),
                                        by = population]
      FDR_by_threshold <- rbind(FDR_by_threshold,
                                temp_table)
    }
  
    FDR_by_threshold[, FDR:= 1-true_positive/total_positive]
    FDR_by_threshold[, method:= 'LS-YRI-train']
  
    #####
  
    FDR_all_method <- rbind(FDR_all_method, FDR_by_threshold)
  
    ####
    FDR_all_method$simulation_index <- simulation_index
    #####
    FDR_all_simulations <- rbind(FDR_all_simulations, FDR_all_method)
  }##simulation_index
  ###########
  
  FDR_average_simulation <- FDR_all_simulations[, .(mean_FDR = mean(FDR),
                                                    mean_total_positive = mean(total_positive),
                                                    mean_true_positive = mean(true_positive)),
                                                by = .(population, threshold, method)]
  FDR_average_simulation[, mean_false_positive := mean_total_positive - mean_true_positive]
  #####
  ggplot()+
    geom_line(aes(x = threshold, y = mean_FDR, 
                  group = as.factor(paste0(population,'-',method)), 
                  color = as.factor(paste0(population,'-',method))),
              data = FDR_average_simulation)+
    geom_point(aes(x = threshold, y = mean_FDR, 
                  group = as.factor(paste0(population,'-',method)), 
                  shape = as.factor(paste0(population,'-',method)),
                  color = as.factor(paste0(population,'-',method))),
               # show.legend = FALSE,
               data = FDR_average_simulation)+
    # geom_line(aes(x = threshold, y = total_positive,
    #               group = as.factor(paste0(population,'-',method)),
    #               color = as.factor(paste0(population,'-',method))),
    #           data = FDR_all_method)+
    theme_bw()+
    # theme(legend.position = 'bottom')+
    xlim(c(-2.5, 2.5))+
    xlab('PGS threshold')+
    ylab('FDR')+
    scale_color_manual(name=" ",
                       values=c("#d72631", "#a2d5c6", 
                                "#077b8a", "#5c3c92",
                                "#070601", "#C2C5C6"))+
    scale_shape_manual(name=" ",
                       values = 1:6)
  
  ###
  filename_prefix <- c('6xx', '7xx', '8xx')
  FDR_plot_file <- paste0('Sim',filename_prefix[simulation_setting],
                          '-testing-FDR.pdf')
  ggsave(filename = FDR_plot_file,
         device = "pdf",
         width = 160, height = 90, units = "mm")
  ###
}#simulation_setting

####
ggplot()+
  geom_line(aes(x = threshold, y = mean_false_positive, 
                group = as.factor(paste0(population,'-',method)), 
                color = as.factor(paste0(population,'-',method))),
            data = FDR_average_simulation[population == 'CEU',])+
  geom_point(aes(x = threshold, y = mean_false_positive, 
                 group = as.factor(paste0(population,'-',method)), 
                 shape = as.factor(paste0(population,'-',method)),
                 color = as.factor(paste0(population,'-',method))),
             show.legend = FALSE,
             data = FDR_average_simulation[population == 'CEU',])+
  theme_bw()+
  xlim(c(-2.5, 2.5))+
  xlab('PGS threshold')+
  ylab('# of false positive')+
  scale_color_manual(name=" ",
                     # labels= 
                     values=c("#d72631", "#a2d5c6", 
                              "#077b8a", "#5c3c92",
                              "#070601", "#C2C5C6"))


# 
# 
# #####
# ggplot()+
#   geom_histogram(aes(x=jls, color=PHENO1, group = PHENO1),
#                  alpha=0.1, position = 'identity',
#                  data = rawPRS_YRI_oneSimulation)+
#   theme_bw()+
#   xlim(c(-4,4))+
#   ggtitle('Joint Lassosum')
# 
# ggplot()+
#   geom_histogram(aes(x=wls, color=PHENO1, group = PHENO1),
#                  alpha=0.1, position = 'identity',
#                  data = rawPRS_YRI_oneSimulation)+
#   theme_bw()+
#   xlim(c(-4,4))+
#   ggtitle('Weighted Lassosum')
# 
# ggplot()+
#   geom_histogram(aes(x=ls.trans, color=PHENO1, group = PHENO1),
#                  alpha=0.1, position = 'identity',
#                  data = rawPRS_YRI_oneSimulation)+
#   theme_bw()+
#   xlim(c(-4,4))+
#   ggtitle('CEU Lassosum')
# 
# ggplot()+
#   geom_histogram(aes(x=ls.cis, color=PHENO1, group = PHENO1),
#                  alpha=0.1, position = 'identity',
#                  data = rawPRS_YRI_oneSimulation)+
#   theme_bw()+
#   xlim(c(-4,4))+
#   ggtitle('YRI Lassosum')
# #####
# 
# 
# 
# 
# 
# 
# 
# ############get the false discovery rate plot
# test.inf[, is.case := as.numeric((PHENO1 == 2))]
# test.inf[, is.case := as.numeric((PHENO1 == 2))]
# 
# FDR.table.bothmethods <- data.table()
# 
# FDR.table <- data.table()
# for(thr in seq(-4, 4, length = 50)){
#   
#   temp.table <- test.inf[JLPGS >= thr, 
#                          .(all.pos = .N,
#                            true.pos = sum(is.case),
#                            thr = thr),
#                          by = pop]
#   FDR.table <- rbind(FDR.table,temp.table)
# }
# FDR.table[, FDR:= 1-true.pos/all.pos]
# FDR.table[, method:= 'JLssum']
# 
# FDR.table.bothmethods <- rbind(FDR.table.bothmethods, FDR.table)
# 
# FDR.table <- data.table()
# for(thr in seq(-4, 4, length = 50)){
#   
#   temp.table <- test.inf[LPGS >= thr, 
#                          .(all.pos = .N,
#                            true.pos = sum(is.case),
#                            thr = thr),
#                          by = pop]
#   FDR.table <- rbind(FDR.table,temp.table)
# }
# FDR.table[, FDR:= 1-true.pos/all.pos]
# FDR.table[, method:= 'Lssum']
# FDR.table.bothmethods <- rbind(FDR.table.bothmethods, FDR.table)
# 
# ggplot()+
#   geom_line(aes(x = thr, y = FDR, 
#                 group = as.factor(paste0(pop,'-',method)), 
#                 color = as.factor(paste0(pop,'-',method))),
#             data = FDR.table.bothmethods)+
#   theme_bw()+
#   xlab('PGS threshold')+
#   scale_color_manual(name=" ",
#                      # labels= 
#                      values=c("#6D90AD","#ADC0ED",
#                               "#B8689a","#F8A8Ca"))
# 
# 
# 
# 
# head(rawSimulationPRS[[30]])
# 
# work.dir <- '/raid6/Tianyu/PRS/SimulationPipeline/Work/Sim-801'
# test.inf.2pop <- data.frame()
# 
# for(anc in c('YRI', 'CEU')){
# ####TRUE OUTCOME###
# test.inf <- fread(paste0(work.dir, '/', anc,'.TST/', anc,'.TST.psam'),
#                   header = T)
# 
# ###Joint Lassosum
# ###read in AUC
# JL.AUC <- get(load(paste0(work.dir, 
#                           '/JointLassoSum/JointLassosum--gamma-0.80-', anc,'.TST-AUC.Rdata')))
# best.beta.index <- which.max(JL.AUC) #find best lambda
# ###read in PGS
# JL.PGS <- get(load(paste0(work.dir, 
#                           '/JointLassoSum/JointLassosum--gamma-0.80-', anc,'.TST-PGS.Rdata')))
# JL.best.PGS <- JL.PGS[, best.beta.index]
# 
# ##merge to the test data.frame
# test.inf$JLPGS <- JL.best.PGS
# 
# 
# ###Lassosum 
# ##read in PGS
# L.PGS <- fread(paste0(work.dir,
#                          '/', anc,'.TST/Lassosum/', anc,'.TST/CEU.TRN-', anc,'.TST-lassosum.score'),
#                header = T)
# test.inf$LPGS <- L.PGS$best.pgs
# 
# ###record testing population
# test.inf$pop <- anc
# 
# test.inf.2pop <- rbind(test.inf.2pop, test.inf)
# }
# 
# 
# 
# ####save the result for plotting
# save(test.inf.2pop, 
#      file = paste0(work.dir, '/HistPlot-PGS.RData'))
# 
# #######local plotting
# library(ggplot2)
# library(data.table)
# test.inf <- data.table(get(load('/Users/tianyu/Documents/ParaTuning/data/HistPlot-PGS.RData')))
# test.inf[,PHENO1 := as.factor(PHENO1)]
# 
# normalize.PGS <- function(x){
#   x <- (x - mean(x))/ (sd(x))
#   return(x)
# }
# 
# ###normalize the PGS scores
# # test.inf[, JLPGS:= (JLPGS - mean(JLPGS))/ (sd(JLPGS))]
# # test.inf[, LPGS:= (LPGS - mean(LPGS))/ (sd(LPGS))]
# 
# test.inf[, JLPGS := normalize.PGS(JLPGS), by = pop]
# test.inf[, LPGS := normalize.PGS(LPGS), by = pop]
# 
# ggplot()+
#   geom_histogram(aes(x=JLPGS, color=PHENO1),
#                  alpha=0.1, position = 'identity',
#                  data = test.inf[pop == 'YRI',])+
#   theme_bw()+
#   xlim(c(-4,4))
# 
# ggplot()+
#   geom_histogram(aes(x=LPGS, color=PHENO1),
#                  alpha=0.1, position = 'identity',
#                  data = test.inf[pop == 'YRI',])+
#   theme_bw()+
#   xlim(c(-4,4))
# 
# ggplot()+
#   geom_histogram(aes(x=LPGS), fill = 'red', color = 'red',
#                  alpha=0.5, position = 'identity',
#                  data = test.inf[PHENO1 == 1,])+
#   geom_histogram(aes(x=JLPGS), fill = 'blue', color = 'blue',
#                  alpha=0.5, position = 'identity',
#                  data = test.inf[PHENO1 == 1,])+
#   theme_bw()+
#   xlim(c(-4,4))
# 
# ############get the false discovery rate plot
# test.inf[, is.case := as.numeric((PHENO1 == 2))]
# test.inf[, is.case := as.numeric((PHENO1 == 2))]
# 
# FDR.table.bothmethods <- data.table()
# 
# FDR.table <- data.table()
# for(thr in seq(-4, 4, length = 50)){
#   
# temp.table <- test.inf[JLPGS >= thr, 
#                        .(all.pos = .N,
#                          true.pos = sum(is.case),
#                          thr = thr),
#                       by = pop]
# FDR.table <- rbind(FDR.table,temp.table)
# }
# FDR.table[, FDR:= 1-true.pos/all.pos]
# FDR.table[, method:= 'JLssum']
# 
# FDR.table.bothmethods <- rbind(FDR.table.bothmethods, FDR.table)
# 
# FDR.table <- data.table()
# for(thr in seq(-4, 4, length = 50)){
#   
#   temp.table <- test.inf[LPGS >= thr, 
#                          .(all.pos = .N,
#                            true.pos = sum(is.case),
#                            thr = thr),
#                          by = pop]
#   FDR.table <- rbind(FDR.table,temp.table)
# }
# FDR.table[, FDR:= 1-true.pos/all.pos]
# FDR.table[, method:= 'Lssum']
# FDR.table.bothmethods <- rbind(FDR.table.bothmethods, FDR.table)
# 
# ggplot()+
#   geom_line(aes(x = thr, y = FDR, 
#                 group = as.factor(paste0(pop,'-',method)), 
#                 color = as.factor(paste0(pop,'-',method))),
#             data = FDR.table.bothmethods)+
#   theme_bw()+
#   xlab('PGS threshold')+
#   scale_color_manual(name=" ",
#                    # labels= 
#                    values=c("#6D90AD","#ADC0ED",
#                             "#B8689a","#F8A8Ca"))
# 
# FDR.table.bothmethods[thr >= -0.1 & thr <= 0,]
# FDR.table.bothmethods[thr >= 0.9 & thr <= 1.1,]
