#### make some naive GWAS data
library(data.table)
ceu.trn <- fread('/raid6/Ron/prs/data/bert_sample/CEU.TRN.PHENO1.glm.logistic.hybrid')
names(ceu.trn)[1] <- 'CHR'
ceu.trn <- ceu.trn[CHR %in% c(20,21),]

ceu.tosave <- ceu.trn[, .(CHR, OR, P)]
ceu.tosave$n <- 20000
ceu.tosave[OR > 1, OR :=2]
ceu.tosave[OR < 1, OR := 0.5]

fwrite(ceu.tosave, file = '/raid6/Tianyu/PRS/BootData/CEU.TRN.SUBSET')

yri.trn <- fread('/raid6/Ron/prs/data/bert_sample/YRI.TRN.PHENO1.glm.logistic.hybrid')
names(yri.trn)[1] <- 'CHR'
yri.trn <- yri.trn[CHR %in% c(20,21),]

yri.tosave <- yri.trn[, .(CHR, OR, P)]
yri.tosave$n <- 4000
yri.tosave[OR > 1, OR :=2]
yri.tosave[OR < 1, OR := 0.5]

fwrite(yri.tosave, file = '/raid6/Tianyu/PRS/BootData/YRI.TRN.SUBSET')
