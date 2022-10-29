results <- data.frame()
anc <- "CEU"
auc <- get(load('/Users/tianyu/Documents/ParaTuning/data/PGSnPHENO_auc_CEU.RData'))
temp.results <- data.frame(lambda.index = 1:10,
                           population = rep(anc, length(auc)),
                           auc = auc,
                           setting = rep("cv", length(auc)))
results <- rbind(results, temp.results)

anc <- "YRI"
auc <- get(load('/Users/tianyu/Documents/ParaTuning/data/PGSnPHENO_auc_YRI.RData'))
temp.results <- data.frame(lambda.index = 1:10,
                           population = rep(anc, length(auc)),
                           auc = auc,
                           setting = rep("cv", length(auc)))
results <- rbind(results, temp.results)

anc <- "CEU"
auc <- get(load('/Users/tianyu/Documents/ParaTuning/data/PGSnPHENO_org_auc_CEU.RData'))
temp.results <- data.frame(lambda.index = 1:10,
                           population = rep(anc, length(auc)),
                           auc = auc,
                           setting = rep("truth", length(auc)))
results <- rbind(results, temp.results)

anc <- "YRI"
auc <- get(load('/Users/tianyu/Documents/ParaTuning/data/PGSnPHENO_org_auc_YRI.RData'))
temp.results <- data.frame(lambda.index = 1:10,
                           population = rep(anc, length(auc)),
                           auc = auc,
                           setting = rep("truth", length(auc)))
results <- rbind(results, temp.results)

library(ggplot2)
ggplot()+
  geom_line(aes(x = lambda.index, y = auc, group = as.factor(paste0(population ,setting)),
                color = as.factor(paste0(population ,setting))), 
            data = results)+
  geom_point(aes(x = lambda.index, y = auc, group = as.factor(paste0(population ,setting)),
                color = as.factor(paste0(population ,setting))), 
            data = results)
