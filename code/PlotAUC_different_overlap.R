library(data.table)
library(pROC)
library(ggplot2)
library(ggpubr)


setwd('/Users/tianyu/Documents/ParaTuning/data/Different_overlap_percentage/')

all_result <- data.table()
for(overlap_percentage in c(100, 80, 60, 40)){
  # for(overlap_percentage in c(40)){ 
  rawSimulationAUC <- readRDS(paste0('YRI-PGS-overlap-', overlap_percentage, '_2023-03-30.RDS'))
  
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
}
all_result[, method := as.factor(method)]
levels(all_result$method) <- c('JLS', 'LS(C2Y)', 'LS(Y2Y)', 'TL', 'WLS')
all_result[, method := factor(method,
                              levels = c('LS(Y2Y)','LS(C2Y)', 'WLS', 'TL', 'JLS'))]

plot_list <- list()
setting_index <- 1
color_vector <- c("#9b5fe0", "#16a4d8", "#8bd346", "#f9a52c","#d64e12",
                  "#B6BeB2", "#efdf48","#60dbe8")
panel_labels <- c('A','B','C','D')
for(current_overlap_percentage in c(100, 80, 60, 40)){
plot_list[[setting_index]] <- ggplot(all_result[overlap_percentage == current_overlap_percentage,]) +
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
  annotate("text", 
           x = "TL", y = 0.5,
           label = paste0("Overlap ", current_overlap_percentage, "%"),
           hjust = 0, vjust = 0)+
  annotate("text", 
         x = "LS(Y2Y)", y = 0.89,
         label = panel_labels[setting_index],
         hjust = 1.8, vjust = 0)
  # ggtitle(paste0("Overlap ", current_overlap_percentage, "%"))

  ##add the performance of lsum CEU
  for_average_lsum_ceu <- readRDS(paste0('YRI-PGS-overlap-', overlap_percentage, '_2023-03-31.RDS'))
  lsum_ceu_aucs <- NULL
  for(repeat_index in 1:length(for_average_lsum_ceu[['CEU']])){
    current_CEU_data <- for_average_lsum_ceu[['CEU']][[repeat_index]]
    lsum_ceu_aucs <- c(lsum_ceu_aucs, roc(current_CEU_data$PHENO, current_CEU_data$LS.CEU)$auc)
  }
  lsum_ceu_aucs_average <- mean(lsum_ceu_aucs)
  plot_list[[setting_index]] <- plot_list[[setting_index]] + geom_hline(yintercept = lsum_ceu_aucs_average,
                                                                        linetype = "dashed", alpha = 0.5)


if(setting_index %% 2 == 0){
  plot_list[[setting_index]] <- plot_list[[setting_index]] +
    ylab("")
}
setting_index <- setting_index + 1
}

mainCombined<- ggarrange(
  plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], 
  # mainCEU,  
  nrow = 2, ncol = 2,
  # labels = c("A","B"),
  common.legend = TRUE, legend = 'none'
) 


mainCombined <- annotate_figure(mainCombined, 
                                left = grid::textGrob("AUC", rot = 90, hjust = -0.5, vjust = 0.5, gp = gpar(cex = 1)),
                                bottom = grid::textGrob("Method", hjust = 0, gp = gpar(cex = 1)))
mainCombined
ggsave(filename = "Different_overlap_AUC.pdf",
       device = "pdf",
       width = 200, height = 180, units = "mm")

# plot_current <- plot_list[[1]] 
# overlap_percentage <- 40


####equal beta case
all_result <- data.table()
# for(overlap_percentage in c(40)){ 
rawSimulationAUC <- readRDS(paste0('YRI-PGS-equal-beta_2023-03-30.RDS'))

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
                              repeart = repeat_index)
    all_result <- rbind(all_result, temp_result)
  }
}

all_result[, method := as.factor(method)]
levels(all_result$method) <- c('JLsum', 'Lsum(CEU)', 'Lsum(YRI)', 'TL', 'Lsum-wt')
all_result[, method := factor(method,
                              levels = c('Lsum(YRI)','Lsum(CEU)', 'Lsum-wt', 'TL', 'JLsum'))]


color_vector <- c("#9b5fe0", "#16a4d8", "#8bd346", "#f9a52c","#d64e12",
                  "#B6BeB2", "#efdf48","#60dbe8")
equal_beta_plot <- ggplot(all_result) +
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
  ylab("Testing AUC")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ggtitle(paste0("Equal beta setting"))+
  theme(legend.position = 'none')

##add the performance of lsum CEU
for_average_lsum_ceu <- readRDS('YRI-PGS-equal-beta_2023-03-31.RDS')
lsum_ceu_aucs <- NULL
for(repeat_index in 1:length(for_average_lsum_ceu[['CEU']])){
  current_CEU_data <- for_average_lsum_ceu[['CEU']][[repeat_index]]
  lsum_ceu_aucs <- c(lsum_ceu_aucs, roc(current_CEU_data$PHENO, current_CEU_data$LS.CEU)$auc)
}
lsum_ceu_aucs_average <- mean(lsum_ceu_aucs)
equal_beta_plot <- equal_beta_plot + geom_hline(yintercept = lsum_ceu_aucs_average,
                                                linetype = "dashed", alpha = 0.5)

equal_beta_plot
ggsave(filename = "Equal_beta_AUC.pdf",
       device = "pdf",
       width = 200, height = 90, units = "mm")

