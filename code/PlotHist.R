library(data.table)

work.dir <- '/raid6/Tianyu/PRS/SimulationPipeline/Work/Sim-801'
test.inf.2pop <- data.frame()

for(anc in c('YRI', 'CEU')){
####TRUE OUTCOME###
test.inf <- fread(paste0(work.dir, '/', anc,'.TST/', anc,'.TST.psam'),
                  header = T)

###Joint Lassosum
###read in AUC
JL.AUC <- get(load(paste0(work.dir, 
                          '/JointLassoSum/JointLassosum--gamma-0.80-', anc,'.TST-AUC.Rdata')))
best.beta.index <- which.max(JL.AUC) #find best lambda
###read in PGS
JL.PGS <- get(load(paste0(work.dir, 
                          '/JointLassoSum/JointLassosum--gamma-0.80-', anc,'.TST-PGS.Rdata')))
JL.best.PGS <- JL.PGS[, best.beta.index]

##merge to the test data.frame
test.inf$JLPGS <- JL.best.PGS


###Lassosum 
##read in PGS
L.PGS <- fread(paste0(work.dir,
                         '/', anc,'.TST/Lassosum/', anc,'.TST/CEU.TRN-', anc,'.TST-lassosum.score'),
               header = T)
test.inf$LPGS <- L.PGS$best.pgs

###record testing population
test.inf$pop <- anc

test.inf.2pop <- rbind(test.inf.2pop, test.inf)
}



####save the result for plotting
save(test.inf.2pop, 
     file = paste0(work.dir, '/HistPlot-PGS.RData'))

#######local plotting
library(ggplot2)
library(data.table)
test.inf <- data.table(get(load('/Users/tianyu/Documents/ParaTuning/data/HistPlot-PGS.RData')))
test.inf[,PHENO1 := as.factor(PHENO1)]

normalize.PGS <- function(x){
  x <- (x - mean(x))/ (sd(x))
  return(x)
}

###normalize the PGS scores
# test.inf[, JLPGS:= (JLPGS - mean(JLPGS))/ (sd(JLPGS))]
# test.inf[, LPGS:= (LPGS - mean(LPGS))/ (sd(LPGS))]

test.inf[, JLPGS := normalize.PGS(JLPGS), by = pop]
test.inf[, LPGS := normalize.PGS(LPGS), by = pop]

ggplot()+
  geom_histogram(aes(x=JLPGS, color=PHENO1),
                 alpha=0.1, position = 'identity',
                 data = test.inf[pop == 'YRI',])+
  theme_bw()+
  xlim(c(-4,4))

ggplot()+
  geom_histogram(aes(x=LPGS, color=PHENO1),
                 alpha=0.1, position = 'identity',
                 data = test.inf[pop == 'YRI',])+
  theme_bw()+
  xlim(c(-4,4))

ggplot()+
  geom_histogram(aes(x=LPGS), fill = 'red', color = 'red',
                 alpha=0.5, position = 'identity',
                 data = test.inf[PHENO1 == 1,])+
  geom_histogram(aes(x=JLPGS), fill = 'blue', color = 'blue',
                 alpha=0.5, position = 'identity',
                 data = test.inf[PHENO1 == 1,])+
  theme_bw()+
  xlim(c(-4,4))

############get the false discovery rate plot
test.inf[, is.case := as.numeric((PHENO1 == 2))]
test.inf[, is.case := as.numeric((PHENO1 == 2))]

FDR.table.bothmethods <- data.table()

FDR.table <- data.table()
for(thr in seq(-4, 4, length = 50)){
  
temp.table <- test.inf[JLPGS >= thr, 
                       .(all.pos = .N,
                         true.pos = sum(is.case),
                         thr = thr),
                      by = pop]
FDR.table <- rbind(FDR.table,temp.table)
}
FDR.table[, FDR:= 1-true.pos/all.pos]
FDR.table[, method:= 'JLssum']

FDR.table.bothmethods <- rbind(FDR.table.bothmethods, FDR.table)

FDR.table <- data.table()
for(thr in seq(-4, 4, length = 50)){
  
  temp.table <- test.inf[LPGS >= thr, 
                         .(all.pos = .N,
                           true.pos = sum(is.case),
                           thr = thr),
                         by = pop]
  FDR.table <- rbind(FDR.table,temp.table)
}
FDR.table[, FDR:= 1-true.pos/all.pos]
FDR.table[, method:= 'Lssum']
FDR.table.bothmethods <- rbind(FDR.table.bothmethods, FDR.table)

ggplot()+
  geom_line(aes(x = thr, y = FDR, 
                group = as.factor(paste0(pop,'-',method)), 
                color = as.factor(paste0(pop,'-',method))),
            data = FDR.table.bothmethods)+
  theme_bw()+
  xlab('PGS threshold')+
  scale_color_manual(name=" ",
                   # labels= 
                   values=c("#6D90AD","#ADC0ED",
                            "#B8689a","#F8A8Ca"))

FDR.table.bothmethods[thr >= -0.1 & thr <= 0,]
FDR.table.bothmethods[thr >= 0.9 & thr <= 1.1,]
