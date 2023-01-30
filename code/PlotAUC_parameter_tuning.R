library(data.table)
library(ggplot2)


setwd('/Users/tianyu/Documents/ParaTuning/data/parameter_tuning_auc/')
raw_parameter_tuning_auc <- get(load('synthetic_data_tuning_auc_gamma_0.20_lambda_0.00150.RData'))

#####

auc_data_frame <- data.table()
for(repeat_index in 1:length(raw_parameter_tuning_auc)){
  
  temporary_data_frame <- raw_parameter_tuning_auc[[repeat_index]]
  temporary_data_frame$simulation_index <- names(raw_parameter_tuning_auc)[repeat_index]
  
  auc_data_frame <- rbind(auc_data_frame,
                          temporary_data_frame)
  
}

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

auc_data_frame <- merge(auc_data_frame, cross_validation_auc_data_frame,
                        by = c("simulation_index", "ancestry"),
                        all.x = T)

###

auc_data_frame_long <- melt(auc_data_frame,
                           id.vars = c("simulation_index", "ancestry"),
                           measure.vars = c("oracle_auc", "tuned_auc"),
                           variable.name = "auc_source", value.name = "auc")

#######


ggplot(auc_data_frame)+
  geom_point(aes(x = oracle_auc, y = tuned_auc,
                 shape = ancestry))+
  geom_abline(slope = 1, intercept = 0)+
  scale_color_manual(values=c("#9b5fe0",
                              # "#16a4d8", 
                              "#60dbe8", "#8bd346", "#efdf48","#f9a52c","#d64e12","#B6BeB2"))+
  xlim(c(0.5, 0.9))+
  ylim(c(0.5, 0.9))+
  theme_bw()

ggsave(filename = "Sim8xx-parameter-tuning-scatter.pdf",
       device = "pdf",
       width = 120, height = 100, units = "mm")

