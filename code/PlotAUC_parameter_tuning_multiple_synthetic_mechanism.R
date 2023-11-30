library(data.table)
library(ggplot2)
library(ggpubr)
#####

setwd('/Users/tianyu/Documents/ParaTuning/data/')
cross_validation_auc <- get(load('simulations-Results-6xx-7xx-8xx01-13-2023.RData'))

cross_validation_auc_data_frame <- data.table()
for(repeat_name in names(cross_validation_auc)){
  
  temporary_data_frame <- cross_validation_auc[[repeat_name]][['jls']]
  formatted_temporary_data_frame <- data.frame(ancestry = temporary_data_frame$TST,
                                               cross_validation_auc = temporary_data_frame$auc,
                                               simulation_index = gsub("Sim-", "", repeat_name))
  cross_validation_auc_data_frame <- rbind(cross_validation_auc_data_frame,
                                           formatted_temporary_data_frame)
  
}

#####



synthetic_parameter_table <- data.frame()
synthetic_parameter_table <- rbind(synthetic_parameter_table, 
                                   data.frame(gamma = "0.20", lambda = "0.0150"))
synthetic_parameter_table <- rbind(synthetic_parameter_table, 
                                   data.frame(gamma = "0.50", lambda = "0.0070"))
synthetic_parameter_table <- rbind(synthetic_parameter_table, 
                                   data.frame(gamma = "0.80", lambda = "0.0070"))
synthetic_parameter_table <- rbind(synthetic_parameter_table, 
                                   data.frame(gamma = "0.80", lambda = "0.0025"))

#####

setwd('/Users/tianyu/Documents/ParaTuning/data/parameter_tuning_auc_600/')

for(synthetic_parameter_index in 1:nrow(synthetic_parameter_table)){
    
  synthetic_data_tuning_auc_file <- paste0("synthetic_data_tuning_auc_gamma_", 
                                           synthetic_parameter_table[synthetic_parameter_index,]$gamma,
                                           "_lambda_",
                                           synthetic_parameter_table[synthetic_parameter_index,]$lambda,
                                           ".RData")
  raw_parameter_tuning_auc <- get(load(synthetic_data_tuning_auc_file))
  # raw_parameter_tuning_auc <- get(load('synthetic_data_tuning_auc_gamma_0.20_lambda_0.0150.RData'))
  
  #####
  
  auc_data_frame <- data.table()
  for(repeat_index in 1:length(raw_parameter_tuning_auc)){
    
    temporary_data_frame <- raw_parameter_tuning_auc[[repeat_index]]
    temporary_data_frame$simulation_index <- names(raw_parameter_tuning_auc)[repeat_index]
    
    auc_data_frame <- rbind(auc_data_frame,
                            temporary_data_frame)
    
  }
  
  if(synthetic_parameter_index == 1){
    auc_data_frame_column_name <- colnames(auc_data_frame)
    auc_data_frame_column_name[which(auc_data_frame_column_name == "tuned_auc")] <- paste0("tuned_auc",
                                                                                           "_gamma_", 
                                                                                           synthetic_parameter_table[synthetic_parameter_index,]$gamma,
                                                                                           "_lambda_",
                                                                                           synthetic_parameter_table[synthetic_parameter_index,]$lambda)
    colnames(auc_data_frame) <- auc_data_frame_column_name
  }else{
    auc_data_frame <- auc_data_frame[, c("ancestry", "tuned_auc", "simulation_index")]
    auc_data_frame_column_name <- colnames(auc_data_frame)
    auc_data_frame_column_name[which(auc_data_frame_column_name == "tuned_auc")] <- paste0("tuned_auc",
                                                                                           "_gamma_", 
                                                                                           synthetic_parameter_table[synthetic_parameter_index,]$gamma,
                                                                                           "_lambda_",
                                                                                           synthetic_parameter_table[synthetic_parameter_index,]$lambda)
    colnames(auc_data_frame) <- auc_data_frame_column_name
  }
  
  
  
  cross_validation_auc_data_frame <- merge(auc_data_frame, cross_validation_auc_data_frame,
                          by = c("simulation_index", "ancestry"),
                          all.x = T)
  
  cross_validation_auc_data_frame[, cross_validation_auc := pmin(cross_validation_auc, oracle_auc)]

}
###

variable_to_melt <- c("cross_validation_auc")
for(synthetic_parameter_index in 1:nrow(synthetic_parameter_table)){
  variable_to_melt <- c(variable_to_melt, 
                        paste0("tuned_auc",
                               "_gamma_", 
                               synthetic_parameter_table[synthetic_parameter_index,]$gamma,
                               "_lambda_",
                               synthetic_parameter_table[synthetic_parameter_index,]$lambda))
}


cross_validation_auc_data_frame_long <- melt(cross_validation_auc_data_frame,
                                             id.vars = c("simulation_index", "ancestry", "oracle_auc"),
                                             # measure.vars = c("oracle_auc", "tuned_auc"),
                                             measure.vars = variable_to_melt,
                                             variable.name = "auc_source", value.name = "auc")

#######
color_vector <- c("#d72631", "#e8dd72","#a8c66c","#1b6535","#f9a52c", "#efdf48","#d64e12")
label_vector <- c("CV","Synthetic (I)",
                  "Synthetic (II)","Synthetic (III)",
                  "Synthetic (IV)")
shape_vector <- c(10, 15:18)

unique_auc_source_name <- unique(cross_validation_auc_data_frame_long$auc_source)
group1_auc_source_name <- unique_auc_source_name[c(1,2,3)]
group2_auc_source_name <- unique_auc_source_name[c(1,4,5)]

both_tuning_scatterplot_group1 <- ggplot(cross_validation_auc_data_frame_long[auc_source %in% group1_auc_source_name, ])+
  geom_point(aes(x = oracle_auc, y = auc,
                 shape = auc_source, color = auc_source),
             alpha = 0.6, size = 3)+
  annotate("text", x = 0.8, y = 0.72,
           label = "YRI",
           hjust = 0, vjust = 0)+
  annotate("text", x = 0.87, y = 0.82,
           label = "CEU",
           hjust = 0, vjust = 0)+
  geom_abline(slope = 1, intercept = 0, size = 0.3)+
  scale_shape_manual(name = NULL,
                     labels = label_vector,
                     values= shape_vector,
                     drop = TRUE)+
  scale_color_manual(name = NULL,
                     labels = label_vector,
                     values= color_vector,
                     drop = TRUE)+
  xlim(c(0.6, 0.9))+
  ylim(c(0.6, 0.9))+
  xlab(NULL)+
  ylab(NULL)+
  # xlab('Oracle Testing AUC')+
  # ylab('Testing AUC of Tuned Models')+
  theme_bw()+
  theme(legend.position = c(0.2, 0.84))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

both_tuning_scatterplot_group2 <- ggplot(cross_validation_auc_data_frame_long[auc_source %in% group2_auc_source_name,])+
  geom_point(aes(x = oracle_auc, y = auc,
                 shape = auc_source, color = auc_source),
             alpha = 0.6, size = 3)+
  # geom_rect(aes(xmin = 0.72, xmax = 0.8, ymin = 0.68, ymax = 0.79), 
  #           alpha = 0, color = "#B0B0B0",)+
  # geom_rect(aes(xmin = 0.83, xmax = 0.87, ymin = 0.82, ymax = 0.87),
  # alpha = 0, color = "#B0B0B0",)+
  annotate("text", x = 0.8, y = 0.72,
           label = "YRI",
           hjust = 0, vjust = 0)+
  annotate("text", x = 0.87, y = 0.82,
           label = "CEU",
           hjust = 0, vjust = 0)+
  geom_abline(slope = 1, intercept = 0, size = 0.3)+
  scale_shape_manual(name = NULL,
                     labels = label_vector[c(1,4,5)],
                     values= shape_vector[c(1,4,5)],
                     drop = TRUE)+
  scale_color_manual(name = NULL,
                     labels = label_vector[c(1,4,5)],
                     values= color_vector[c(1,4,5)],
                     drop = TRUE)+
  xlim(c(0.6, 0.9))+
  ylim(c(0.6, 0.9))+
  xlab(NULL)+
  ylab(NULL)+
  theme_bw()+
  theme(legend.position = c(0.2, 0.84))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# YRI_tuning_scatterplot <- ggplot(cross_validation_auc_data_frame_long[ancestry == 'YRI'])+
#   geom_point(aes(x = oracle_auc, y = auc,
#                  shape = auc_source, color = auc_source),
#              alpha = 0.8, size = 3)+
#   geom_abline(slope = 1, intercept = 0, size = 0.3)+
#   scale_shape_manual(name = 'Method',
#                      labels = label_vector,
#                      values= shape_vector)+
#   scale_color_manual(name = 'Method',
#                      labels = label_vector,
#                      values=color_vector)+
#   xlim(c(0.69, 0.78))+
#   ylim(c(0.69, 0.78))+
#   xlab('Oracle Testing AUC')+
#   ylab('Testing AUC of Tuned Models')+
#   theme_bw()+
#   theme(legend.position = 'bottom')
# 
# CEU_tuning_scatterplot <- ggplot(cross_validation_auc_data_frame_long[ancestry == 'CEU'])+
#   geom_point(aes(x = oracle_auc, y = auc,
#                  shape = auc_source, color = auc_source),
#              alpha = 0.8, size = 3)+
#   geom_abline(slope = 1, intercept = 0, size = 0.3)+
#   scale_shape_manual(name = 'Method',
#                      labels = label_vector,
#                      values= shape_vector)+
#   scale_color_manual(name = 'Method',
#                      labels = label_vector,
#                      values=color_vector)+
#   xlim(c(0.81, 0.875))+
#   ylim(c(0.81, 0.875))+
#   xlab('Oracle Testing AUC')+
#   ylab('Testing AUC of Tuned Models')+
#   theme_bw()+
#   theme(legend.position = 'bottom')

two_population_scatterplot <- ggarrange(
  both_tuning_scatterplot_group1,both_tuning_scatterplot_group2,  
  nrow = 1,
  # labels = c("A","B"),
  common.legend = FALSE
) 

two_population_scatterplot <- annotate_figure(two_population_scatterplot, 
                                left = grid::textGrob("AUC of Selected Models", rot = 90, hjust = 0.5, vjust = 0.5, gp = grid::gpar(cex = 1)),
                                bottom = grid::textGrob("Oracle AUC", hjust = 0.3, gp = grid::gpar(cex = 1)))

two_population_scatterplot

ggsave(filename = "Sim6xx-parameter-tuning-scatter.pdf",
       device = "pdf",
       width = 200, height = 95, units = "mm")

