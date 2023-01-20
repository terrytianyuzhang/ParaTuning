##play with gwas
library(data.table)

ParameterTuningDirectory <- "/raid6/Tianyu/PRS/SimulationPipeline/Work/Sim-800/ParameterTuningData/"
#####
chr <- 20
unlabeled_fam_file <- paste0(ParameterTuningDirectory,
                             "CHR/YRI-chr", chr, "synthetic-train.fam")
unlabeled_fam <- fread(unlabeled_fam_file)

#####
anc <- 'YRI'
TrainSampleIndexRFile <- paste0(ParameterTuningDirectory, "/", anc, "-synthetic-train-index.Rdata")
train.index <- get(load(TrainSampleIndexRFile))

####
SyntheticYFile <- paste0(ParameterTuningDirectory, '/',
                         anc,'-SyntheticY', '.RData')
SyntheticY <- get(load(SyntheticYFile))
SyntheticY <- SyntheticY[train.index]

####
# unlabeled_fam$V6 <- rbinom(NROW(unlabeled_fam),size = 1, prob = 0.5) + 1
unlabeled_fam$V6 <- SyntheticY + 1
fwrite(unlabeled_fam, file = unlabeled_fam_file, sep = ' ', col.names = FALSE)

plink2 <- "/usr/local/bin/plink2"
dir.create(paste0(ParameterTuningDirectory, 'Assoc/'),
           showWarnings = F,recursive = T)
plink2.command=paste(plink2,"--nonfounders","--allow-no-sex",
                     # "--pfile",paste(set.dir,set,sep=""),
                    # "--pfile", "/raid6/Tianyu/PRS/SimulationPipeline/Work/Sim-800/YRI.TRN/YRI.TRN",
                    "--bfile", paste0(ParameterTuningDirectory, "CHR/YRI-chr", chr,"synthetic-train"),
                     "--glm","allow-no-covars","omit-ref",
                     # "--out",paste(set.dir,"Assoc/",set,sep=""),
                    # "--out", "/raid6/Tianyu/PRS/SimulationPipeline/Work/Sim-800/CEU.TRN/Assoc/CEU.TRN",
                    "--out", paste0(ParameterTuningDirectory, "Assoc/YRI-chr", chr,"synthetic-train"),
                     sep=" ")
system(plink2.command)

# > paste(set.dir,set,sep="")
# [1] "/raid6/Tianyu/PRS/SimulationPipeline/Work/Sim-800/CEU.TRN/CEU.TRN"
# > paste(set.dir,"Assoc/",set,sep="")
# [1] "/raid6/Tianyu/PRS/SimulationPipeline/Work/Sim-800/CEU.TRN/Assoc/CEU.TRN"

######
all_GWAS <- data.table()
chr <- 20
chr20_GWAS <- fread(paste0(ParameterTuningDirectory, 
                           "Assoc/YRI-chr", chr,
                           "synthetic-train.PHENO1.glm.logistic.hybrid"),
                    header = T)
all_GWAS <- rbind(all_GWAS, chr20_GWAS)

chr <- 21
chr21_GWAS <- fread(paste0(ParameterTuningDirectory, 
                           "Assoc/YRI-chr", chr,
                           "synthetic-train.PHENO1.glm.logistic.hybrid"),
                    header = T)
all_GWAS <- rbind(all_GWAS, chr21_GWAS)

fwrite(all_GWAS, paste0(ParameterTuningDirectory, 
                        "Assoc/YRI-",
                        "synthetic-train.PHENO1.glm.logistic.hybrid"),
       sep = ' ')

