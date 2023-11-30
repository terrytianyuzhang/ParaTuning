library(data.table)
library(pROC)
library(ggplot2)
library(ggpubr)


setwd('/Users/tianyu/Documents/ParaTuning/data/Different_overlap_percentage/')


plot_list <- list()
for(population_name in c('YRI', 'CEU')){
  all_result <- data.table()
  rawSimulationAUC <- readRDS('YRI-CEU-PGS-equal-n_2023-04-05.RDS')[[population_name]]
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
  
  if(population_name == 'YRI'){
    all_result[, method := as.factor(method)]
    levels(all_result$method) <- c('JLS', 'LS(C2Y)', 'LS(Y2Y)', 'TL', 'WLS')
    all_result[, method := factor(method,
                                  levels = c('LS(Y2Y)','LS(C2Y)', 'WLS', 'TL', 'JLS'))]
    
    
  }else{
    all_result[, method := as.factor(method)]
    levels(all_result$method) <- c('JLS', 'LS(C2C)', 'LS(Y2C)', 'WLS')
    all_result[, method := factor(method,
                                  levels = c('LS(Y2C)','LS(C2C)', 'WLS', 'JLS'))]
    
    
  }
  
  color_vector <- c("#9b5fe0", "#16a4d8", "#8bd346", "#d64e12",
                    "#B6BeB2", "#efdf48","#60dbe8")
  plot_list[[population_name]] <- ggplot(all_result[method != 'TL',]) +
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
    ylab(paste0("AUC (", population_name, ")"))+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    # ggtitle(paste0("equal sample size"))+
    theme(legend.position = 'none')
}

same_n_combined <- ggarrange(
  plot_list[[1]], plot_list[[2]] ,
  # mainCEU,  
  nrow = 1, ncol = 2,
  # labels = c("A","B"),
  common.legend = TRUE, legend = 'none'
) 
same_n_combined

ggsave(filename = "same_sample_size_AUC.pdf",
       device = "pdf",
       width = 200, height = 90, units = "mm")

##




###FDR PLOT

normalize_PRS <- function(x){
  x <- (x - mean(x))/ (sd(x))
  return(x)
}

rawSimulationAUC <- readRDS('YRI-CEU-PGS-equal-n_2023-04-05.RDS')

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
  }else  if(old_name == 'LS.YRICEU'){
    new_name <- 'LS(Y2C)'
  }else  if(old_name == 'JLSCEU'){
    new_name <- 'JLS(Y&C2C)'
  }else{
    new_name <- old_name
    }
  return(new_name)
}
for(row_index in 1:nrow(FDR_average_simulation)){
  FDR_average_simulation[row_index, plot_legend := as.factor(change_name(plot_legend))]
}
FDR_average_for_plot <- FDR_average_simulation[plot_legend %in% c('LS(Y2Y)', 'LS(Y2C)','LS(C2C)','LS(C2Y)', 'JLS(Y&C2C)', 'JLS(Y&C2Y)')]

break_names <-  c('LS(C2Y)','LS(Y2C)', 'LS(Y2Y)','LS(C2C)','JLS(Y&C2Y)', 'JLS(Y&C2C)')

ggplot( data = FDR_average_for_plot)+
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
  xlab('PGS threshold')+
  ylab('FDR')+
  scale_color_manual(name=" ",
                     values=c("#d72631","#a2d5c6", "#077b8a", "#5c3c92",
                              "#070601", "#C2C5C6"),
                     breaks = break_names)+
  scale_shape_manual(name=" ",
                     values = 1:6,
                     breaks = break_names)+
  # ggtitle(paste0("Overlap percentage = ", current_overlap_percentage))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(filename = "same_sample_size_FDR.pdf",
       device = "pdf",
       width = 120, height = 90, units = "mm")

