library(data.table)
library(pROC)
library(ggplot2)
library(ggpubr)


setwd('/Users/tianyu/Documents/ParaTuning/data/Different_overlap_percentage/')

all_result <- data.table()

rawSimulationAUC <- readRDS('YRI-PGS-low-h2_2023-04-03.RDS')[['YRI']]
# rawSimulationAUC <- readRDS(paste0('YRI-PGS-overlap-', overlap_percentage, '_2023-03-30.RDS'))

for(repeat_index in 0:(length(rawSimulationAUC)-1)){
  repitition_name <- paste0('rep.', repeat_index)
  current_label <- rawSimulationAUC[[repitition_name]]$PHENO
  
  method_names <- names(rawSimulationAUC[[repitition_name]])[c(-1, -2)]
  
  for(method_name_index in 1:length(method_names)){
    method_name <- method_names[method_name_index]
    current_score <- rawSimulationAUC[[repitition_name]][[method_name]]
    current_auc <- roc(current_label, current_score)$auc
    
    temp_result <- data.frame(auc = current_auc,
                              method = method_name,
                              overlap_percentage = overlap_percentage,
                              repeart = repeat_index)
    all_result <- rbind(all_result, temp_result)
  }
}

all_result[, method := as.factor(method)]
levels(all_result$method) <- c('JLS', 'LS(C2Y)', 'LS(Y2Y)', 'TL', 'WLS')
all_result[, method := factor(method,
                              levels = c('LS(Y2Y)','LS(C2Y)', 'WLS', 'TL', 'JLS'))]

color_vector <- c("#9b5fe0", "#16a4d8", "#8bd346", "#f9a52c","#d64e12",
                  "#B6BeB2", "#efdf48","#60dbe8")

low_heritability_auc_plot <- ggplot(all_result) +
  geom_violin(aes(x=method, y=auc, fill = method),
              trim=T, alpha = 0.4)+
  geom_point(aes(x=method, y=auc,
                 color = method),  alpha = 0.8,
             position = position_jitter(width = 0.05))+
  scale_fill_manual(values= color_vector)+
  scale_color_manual(values= color_vector)+
  ylim(c(0.5,0.9))+
  # xlab("Method")+
  xlab(NULL)+
  ylab(NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  # ggtitle(paste0("low heritability setting"))+
  theme(legend.position = 'none')

##add the performance of lsum CEU
for_average_lsum_ceu <- readRDS('YRI-PGS-low-h2_2023-04-03.RDS')
lsum_ceu_aucs <- NULL
for(repeat_index in 1:length(for_average_lsum_ceu[['CEU']])){
  current_CEU_data <- for_average_lsum_ceu[['CEU']][[repeat_index]]
  lsum_ceu_aucs <- c(lsum_ceu_aucs, roc(current_CEU_data$PHENO, current_CEU_data$LS.CEU)$auc)
}
lsum_ceu_aucs_average <- mean(lsum_ceu_aucs)
low_heritability_auc_plot <- low_heritability_auc_plot + geom_hline(yintercept = lsum_ceu_aucs_average,
                                                                      linetype = "dashed", alpha = 0.5)
# low_heritability_auc_plot
# ggsave(filename = "low_heritability_AUC.pdf",
#        device = "pdf",
#        width = 200, height = 90, units = "mm")

###AUC plot for the 7xx setting
setwd('/Users/tianyu/Documents/ParaTuning/data/')
rawSimulationAUC <- get(load('simulations-Results-6xx-7xx-8xx01-13-2023.RData'))

clean_simulation_AUC <- data.table()
for(setting_index in 11:20){
  ls_result <- rawSimulationAUC[[setting_index]]$ls
  ls_result_formated <- ls_result[,c('TRN','TST','auc')]
  ls_result_formated$setting <- names(rawSimulationAUC)[setting_index]
  ls_result_formated$method <- 'LS'
  
  clean_simulation_AUC <- rbind(clean_simulation_AUC, 
                                ls_result_formated)
  
  wls_result <- rawSimulationAUC[[setting_index]]$wls
  wls_result_formated <- wls_result[,c('TST','auc')]
  wls_result_formated$TRN <- 'Y&C'
  wls_result_formated$setting <- names(rawSimulationAUC)[setting_index]
  wls_result_formated$method <- 'WLS'
  
  clean_simulation_AUC <- rbind(clean_simulation_AUC, 
                                wls_result_formated)
  
  jls_result <- rawSimulationAUC[[setting_index]]$jls
  jls_result_formated <- jls_result[,c('TST','auc')]
  jls_result_formated$TRN <- 'Y&C'
  jls_result_formated$setting <- names(rawSimulationAUC)[setting_index]
  jls_result_formated$method <- 'JLS'
  
  clean_simulation_AUC <- rbind(clean_simulation_AUC, 
                                jls_result_formated)
  
  tl_result <- rawSimulationAUC[[setting_index]]$tl
  tl_result_formated <- tl_result[,c('TST','auc')]
  tl_result_formated$TRN <- 'Y&C'
  tl_result_formated$setting <- names(rawSimulationAUC)[setting_index]
  tl_result_formated$method <- 'TL'
  
  clean_simulation_AUC <- rbind(clean_simulation_AUC, 
                                tl_result_formated)
}

all_result <- clean_simulation_AUC[TST == 'YRI',]
all_result[method == 'LS', method := paste0(method, TRN)]
all_result[, method := as.factor(method)]
levels(all_result$method) <- c('JLS', 'LS(C2Y)', 'LS(Y2Y)', 'TL', 'WLS')
all_result[, method := factor(method,
                              levels = c('LS(Y2Y)','LS(C2Y)', 'WLS', 'TL', 'JLS'))]

color_vector <- c("#9b5fe0", "#16a4d8", "#8bd346", "#f9a52c","#d64e12",
                  "#B6BeB2", "#efdf48","#60dbe8")

series_700_auc_plot <- ggplot(all_result) +
  geom_violin(aes(x=method, y=auc, fill = method),
              trim=T, alpha = 0.4)+
  geom_point(aes(x=method, y=auc,
                 color = method),  alpha = 0.8,
             position = position_jitter(width = 0.05))+
  scale_fill_manual(values= color_vector)+
  scale_color_manual(values= color_vector)+
  ylim(c(0.5,0.9))+
  # xlab("Method")+
  xlab(NULL)+
  ylab(NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  # ggtitle(paste0("low heritability setting"))+
  theme(legend.position = 'none')

lsum_ceu_aucs_average <- mean(clean_simulation_AUC[TST == 'CEU' & method == 'LS' & TRN == 'CEU',]$auc)
series_700_auc_plot <- series_700_auc_plot + geom_hline(yintercept = lsum_ceu_aucs_average,
                                                                    linetype = "dashed", alpha = 0.5)

mainCombined<- ggarrange(
  low_heritability_auc_plot, series_700_auc_plot,
  # mainCEU,  
  nrow = 1, ncol = 2,
  # labels = c("A","B"),
  common.legend = TRUE, legend = 'none'
) 


mainCombined <- annotate_figure(mainCombined, 
                                left = grid::textGrob("AUC", rot = 90, hjust = -0.5, vjust = 0.5, gp = gpar(cex = 1)),
                                bottom = grid::textGrob("Method", hjust = 0, gp = gpar(cex = 1)))
ggsave(filename = "low_heritability_AUC.pdf",
       device = "pdf",
       width = 200, height = 90, units = "mm")

#######FDR plot
normalize_PRS <- function(x){
  x <- (x - mean(x))/ (sd(x))
  return(x)
}

rawSimulationAUC <- readRDS('YRI-PGS-low-h2_2023-04-03.RDS')

FDR_all_simulations <- data.table()
for(testing_population in c('YRI','CEU')){
  # testing_population <- 'YRI'
  for(repeat_index in 0:(length(rawSimulationAUC[[testing_population]])-1)){
    
    repitition_name <- paste0('rep.', repeat_index)
    
    raw_PRS_onesimulation <- data.table(rawSimulationAUC[[testing_population]][[repitition_name]])
    raw_PRS_onesimulation <- melt(raw_PRS_onesimulation, 
                                  id.vars = c("ID", "PHENO"), 
                                  variable.name = "method", 
                                  value.name = "pgs_score")
    raw_PRS_onesimulation[, PHENO:= factor(PHENO)]
    raw_PRS_onesimulation[, is.case := as.numeric((PHENO == 2))]
    raw_PRS_onesimulation[, pgs_score := normalize_PRS(pgs_score), by = "method"]
    
    method_names <- names(rawSimulationAUC[[testing_population]][[repitition_name]])[c(-1, -2)]
    
    
    ####
    FDR_all_method <- data.table()
    for(method_name_index in 1:length(method_names)){
      method_name <- method_names[method_name_index]
      #####
      FDR_by_threshold <- data.table()
      
      for(threshold in seq(-4, 4, length = 50)){
        temp_table <- raw_PRS_onesimulation[method == method_name & pgs_score >= threshold, 
                                            .(total_positive = .N,
                                              true_positive = sum(is.case),
                                              threshold = threshold)]
        FDR_by_threshold <- rbind(FDR_by_threshold,
                                  temp_table)
      }
      
      FDR_by_threshold[, FDR:= 1 - true_positive/total_positive]
      FDR_by_threshold[, method:= method_name]
      
      #####
      
      FDR_all_method <- rbind(FDR_all_method, FDR_by_threshold)
    }#one simulation
    FDR_all_method$simulation_index <- repeat_index
    FDR_all_method$testing_population <- testing_population
    FDR_all_simulations <- rbind(FDR_all_simulations, FDR_all_method)
  }#all simulations in one setting
}#both populations


FDR_average_simulation <- FDR_all_simulations[, .(mean_FDR = mean(FDR),
                                                  mean_total_positive = mean(total_positive),
                                                  mean_true_positive = mean(true_positive)),
                                              by = .(threshold, method, testing_population)]
FDR_average_simulation[, mean_false_positive := mean_total_positive - mean_true_positive]
FDR_average_simulation[, plot_legend := as.factor(paste0(method, testing_population))]

change_name <- function(old_name){
  if(old_name == 'JLSYRI'){
    new_name <- 'JLS(Y&C2Y)'
  }else if(old_name == 'LS.CEUYRI'){
    new_name <- 'LS(C2Y)'
  }else if(old_name == 'TLYRI'){
    new_name <- 'TL(Y&C2Y)'
  }else if(old_name == 'LS.CEUCEU'){
    new_name <- 'LS(C2C)'
  }else if(old_name == 'LS.YRIYRI'){
    new_name <- 'LS(Y2Y)'
  }else  if(old_name == 'WLSYRI'){
    new_name <- 'WLS(Y&C2Y)'
  }
  return(new_name)
}
for(row_index in 1:nrow(FDR_average_simulation)){
  FDR_average_simulation[row_index, plot_legend := as.factor(change_name(plot_legend))]
}

#####
break_names <- method_names
low_heritability_FDR <- ggplot( data = FDR_average_simulation[plot_legend != 'WLS(Y&C2Y)'])+
  geom_line(aes(x = threshold, y = mean_FDR, 
                group = plot_legend, 
                color = plot_legend))+
  geom_point(aes(x = threshold, y = mean_FDR, 
                 group = plot_legend, 
                 shape = plot_legend,
                 color = plot_legend))+
  theme_bw()+
  # theme(legend.position = 'bottom')+
  xlim(c(-2.5, 2.5))+
  ylim(c(0, 0.5))+
  xlab(NULL)+
  ylab(NULL)+
  scale_color_manual(name=" ",
                     values=c("#a2d5c6", "#077b8a", "#5c3c92",
                              "#d72631","#070601", "#C2C5C6"))+
  scale_shape_manual(name=" ",
                     values = 1:6)+
  # ggtitle(paste0("Overlap percentage = ", current_overlap_percentage))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# ggsave(filename = "low_heritability_FDR.pdf",
#        device = "pdf",
#        width = 120, height = 90, units = "mm")
rawSimulationAUC <- readRDS('YRI-PGS-different-h2_2023-04-07.RDS')

FDR_all_simulations <- data.table()
for(testing_population in c('YRI','CEU')){
  # testing_population <- 'YRI'
  for(repeat_index in 0:(length(rawSimulationAUC[[testing_population]])-1)){
    
    repitition_name <- paste0('rep.', repeat_index)
    
    raw_PRS_onesimulation <- data.table(rawSimulationAUC[[testing_population]][[repitition_name]])
    raw_PRS_onesimulation <- melt(raw_PRS_onesimulation, 
                                  id.vars = c("ID", "PHENO"), 
                                  variable.name = "method", 
                                  value.name = "pgs_score")
    raw_PRS_onesimulation[, PHENO:= factor(PHENO)]
    raw_PRS_onesimulation[, is.case := as.numeric((PHENO == 2))]
    raw_PRS_onesimulation[, pgs_score := normalize_PRS(pgs_score), by = "method"]
    
    method_names <- names(rawSimulationAUC[[testing_population]][[repitition_name]])[c(-1, -2)]
    
    
    ####
    FDR_all_method <- data.table()
    for(method_name_index in 1:length(method_names)){
      method_name <- method_names[method_name_index]
      #####
      FDR_by_threshold <- data.table()
      
      for(threshold in seq(-4, 4, length = 50)){
        temp_table <- raw_PRS_onesimulation[method == method_name & pgs_score >= threshold, 
                                            .(total_positive = .N,
                                              true_positive = sum(is.case),
                                              threshold = threshold)]
        FDR_by_threshold <- rbind(FDR_by_threshold,
                                  temp_table)
      }
      
      FDR_by_threshold[, FDR:= 1 - true_positive/total_positive]
      FDR_by_threshold[, method:= method_name]
      
      #####
      
      FDR_all_method <- rbind(FDR_all_method, FDR_by_threshold)
    }#one simulation
    FDR_all_method$simulation_index <- repeat_index
    FDR_all_method$testing_population <- testing_population
    FDR_all_simulations <- rbind(FDR_all_simulations, FDR_all_method)
  }#all simulations in one setting
}#both populations


FDR_average_simulation <- FDR_all_simulations[, .(mean_FDR = mean(FDR),
                                                  mean_total_positive = mean(total_positive),
                                                  mean_true_positive = mean(true_positive)),
                                              by = .(threshold, method, testing_population)]
FDR_average_simulation[, mean_false_positive := mean_total_positive - mean_true_positive]
FDR_average_simulation[, plot_legend := as.factor(paste0(method, testing_population))]

change_name <- function(old_name){
  if(old_name == 'JLSYRI'){
    new_name <- 'JLS(Y&C2Y)'
  }else if(old_name == 'LS.CEUYRI'){
    new_name <- 'LS(C2Y)'
  }else if(old_name == 'TLYRI'){
    new_name <- 'TL(Y&C2Y)'
  }else if(old_name == 'LS.CEUCEU'){
    new_name <- 'LS(C2C)'
  }else if(old_name == 'LS.YRIYRI'){
    new_name <- 'LS(Y2Y)'
  }else  if(old_name == 'WLSYRI'){
    new_name <- 'WLS(Y&C2Y)'
  }
  return(new_name)
}
for(row_index in 1:nrow(FDR_average_simulation)){
  FDR_average_simulation[row_index, plot_legend := as.factor(change_name(plot_legend))]
}

#####
break_names <- method_names
series_700_FDR <- ggplot( data = FDR_average_simulation[plot_legend != 'WLS(Y&C2Y)'])+
  geom_line(aes(x = threshold, y = mean_FDR, 
                group = plot_legend, 
                color = plot_legend))+
  geom_point(aes(x = threshold, y = mean_FDR, 
                 group = plot_legend, 
                 shape = plot_legend,
                 color = plot_legend))+
  theme_bw()+
  # theme(legend.position = 'bottom')+
  xlim(c(-2.5, 2.5))+
  ylim(c(0, 0.5))+
  xlab(NULL)+
  ylab(NULL)+
  scale_color_manual(name=" ",
                     values=c("#a2d5c6", "#077b8a", "#5c3c92",
                              "#d72631","#070601", "#C2C5C6"))+
  scale_shape_manual(name=" ",
                     values = 1:6)+
  # ggtitle(paste0("Overlap percentage = ", current_overlap_percentage))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


FDRCombined<- ggarrange(
  low_heritability_FDR, series_700_FDR,
  # mainCEU,  
  nrow = 1, ncol = 2,
  # labels = c("A","B"),
  common.legend = TRUE, legend = 'top'
) 

FDRCombined <- annotate_figure(FDRCombined, 
                                left = grid::textGrob("FDR", rot = 90, hjust = -0.5, vjust = 0.5, gp = gpar(cex = 1)),
                                bottom = grid::textGrob("PGS Threshold", hjust = 0.5, gp = gpar(cex = 1)))
ggsave(filename = "low_heritability_FDR.pdf",
       device = "pdf",
       width = 200, height = 90, units = "mm")

 