library(data.table)
library(ggplot2)
library(ggpubr)

normalize_PRS <- function(x){
  x <- (x - mean(x))/ (sd(x))
  return(x)
}
setwd('/Users/tianyu/Documents/ParaTuning/data/Different_overlap_percentage/')

# plot_list <- list()
fairness_data <- data.frame()
setting_index <- 1
for(current_overlap_percentage in c(100, 80, 60, 40)){
  print(current_overlap_percentage)
  rawSimulationAUC <- readRDS(paste0('YRI-PGS-overlap-', current_overlap_percentage, '_2023-03-31.RDS'))
  
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
          raw_PRS_onesimulation[, predicted_case := as.numeric(pgs_score >= threshold)]
          temp_table <- raw_PRS_onesimulation[method == method_name,
                                              .(true_positive = sum(predicted_case & is.case),
                                                false_positive = sum(predicted_case & !is.case),
                                                true_negative = sum(!predicted_case & !is.case),
                                                false_negative = sum(!predicted_case & is.case),
                                                threshold = threshold)]
# 
#           temp_table <- raw_PRS_onesimulation[method == method_name & pgs_score >= threshold,
#                                               .(total_positive = .N,
#                                                 true_positive = sum(is.case),
#                                                 threshold = threshold)]
          FDR_by_threshold <- rbind(FDR_by_threshold, temp_table)
        }
        
        ##total_positive is the number of positive subjects the algorithm predicts 
        # FDR_by_threshold[, FDR:= 1 - true_positive/total_positive]
        FDR_by_threshold[, FDR := false_positive/(true_positive + false_positive)]
        FDR_by_threshold[, FPR:= false_positive/(false_positive + true_negative)]
        FDR_by_threshold[, FNR:= false_negative/(true_positive + false_negative)]
        FDR_by_threshold[, method:= method_name]
        
        #####
        
        FDR_all_method <- rbind(FDR_all_method, FDR_by_threshold)
      }#one simulation
      FDR_all_method$simulation_index <- repeat_index
      FDR_all_method$testing_population <- testing_population
      FDR_all_simulations <- rbind(FDR_all_simulations, FDR_all_method)
    }#all simulations in one setting
  }#both populations


  # FDR_average_simulation <- FDR_all_simulations[, .(mean_FDR = mean(FDR),
  #                                                   mean_total_positive = mean(total_positive),
  #                                                   mean_true_positive = mean(true_positive)),
  #                                               by = .(threshold, method, testing_population)]
  # FDR_average_simulation[, mean_false_positive := mean_total_positive - mean_true_positive]
  
  FDR_average_simulation <- FDR_all_simulations[, .(mean_FDR = mean(FDR),
                                                    mean_FPR = mean(FPR),
                                                    mean_FNR = mean(FNR),
                                                    mean_true_positive = mean(true_positive),
                                                    mean_false_positive = mean(false_positive),
                                                    mean_true_negative = mean(true_negative),
                                                    mean_false_negative = mean(false_negative)),
                                                by = .(threshold, method, testing_population)]
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
  }
  for(row_index in 1:nrow(FDR_average_simulation)){
    FDR_average_simulation[row_index, plot_legend := as.factor(change_name(plot_legend))]
  }
  
  FDR_average_simulation$overlap_percentage <- current_overlap_percentage
  fairness_data <- rbind(fairness_data, FDR_average_simulation)
  # #####
  # panel_labels <- c('A','B','C','D')
  # break_names <- method_names
  # plot_list[[setting_index]] <- ggplot( data = FDR_average_simulation[plot_legend != 'WLS(Y&C2Y)'])+
  #   geom_line(aes(x = threshold, y = mean_FDR, 
  #                 group = plot_legend, 
  #                 color = plot_legend))+
  #   geom_point(aes(x = threshold, y = mean_FDR, 
  #                  group = plot_legend, 
  #                  shape = plot_legend,
  #                  color = plot_legend))+
  #   theme_bw()+
  #   # theme(legend.position = 'bottom')+
  #   xlim(c(-2.5, 2.5))+
  #   ylim(c(0, 0.5))+
  #   xlab(NULL)+
  #   ylab(NULL)+
  #   scale_color_manual(name=" ",
  #                      values=c("#a2d5c6", "#077b8a", "#5c3c92",
  #                               "#d72631","#070601", "#C2C5C6"))+
  #   scale_shape_manual(name=" ",
  #                      values = 1:6)+
  #   # ggtitle(paste0("Overlap percentage = ", current_overlap_percentage))+
  #   theme(panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank())+
  #   annotate("text", 
  #            x = -2.5, y = 0,
  #            label = paste0("Overlap ", current_overlap_percentage, "%"),
  #            hjust = 0, vjust = 0)+
  #   annotate("text", 
  #            x = 2.4, y = 0.49,
  #            label = panel_labels[setting_index],
  #            hjust = 1.8, vjust = 0)
  # setting_index <- setting_index + 1
}#all settings

fairness_data_file <- '/Users/tianyu/Documents/ParaTuning/data/Different_overlap_percentage/fairness_data.RData'
save(fairness_data, file = fairness_data_file)

fairness_data_file_csv <- fairness_data_file <- '/Users/tianyu/Documents/ParaTuning/data/Different_overlap_percentage/fairness_data.csv'
write.csv(fairness_data[method != 'WLS', ], file = fairness_data_file_csv)

fairness_data <- get(load(fairness_data_file))
to_present_data_file <- '/Users/tianyu/Documents/ParaTuning/data/Different_overlap_percentage/fairness_data_at_2_overlap_60.csv'
to_present_data <- fairness_data[overlap_percentage == 60 & method != 'WLS'& threshold > 2 & threshold < 2.05,]
write.csv(to_present_data, file = to_present_data_file)
# mainCombined<- ggarrange(
#   plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], 
#   # mainCEU,  
#   nrow = 2, ncol = 2,
#   # labels = c("A","B"),
#   common.legend = TRUE
# ) 
# mainCombined
# 
# mainCombined <- annotate_figure(mainCombined, 
#                                 left = grid::textGrob("FDR", rot = 90, hjust = 0, vjust = 0.5, gp = gpar(cex = 1)),
#                                 bottom = grid::textGrob("PGS Threshold", hjust = 0.3, gp = gpar(cex = 1)))
# mainCombined
# ggsave(filename = "Different_overlap_FDR.pdf",
#        device = "pdf",
#        width = 180, height = 180, units = "mm")
# 
# 
# ###report some numbers
# FDR_average_simulation[threshold > 2 & threshold < 2.05,]
