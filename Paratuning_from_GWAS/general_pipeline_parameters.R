###This piece of information contains the location of important files, softwares
###They do not change between replicates, so we collect here for easier reference

##location of the plink software
plink <- "/usr/local/bin/plink"
plink2 <- "/usr/local/bin/plink2"

##the directories of the main simulation and the parameter tuning pipeline
main_simulation_pipeline_directory <- "/raid6/Tianyu/PRS/SimulationPipeline/"
parameter_tuning_pipeline_directory <- "/raid6/Tianyu/PRS/Paratuning-pipeline/"
data_directory <- paste0(main_simulation_pipeline_directory,"/Data/") 
reference_genotype_directory <- paste0(data_directory,"/Reference-LDblocks/")     

####
main.dir <- "/raid6/Tianyu/PRS/SimulationPipeline/"
work.dir <- paste0("/raid6/Tianyu/PRS/parameter_tuning_from_GWAS/Sim-", i.sim, "/") #"/raid6/Tianyu/PRS/SimulationPipeline/Work/Sim-800/"

####
ld.populations=c("EUR.hg38","AFR.hg38"); names(ld.populations)=c("CEU","YRI")

####
# information on the number of samples
sample_sizes <- list()
sample_sizes[["CEU"]]=data.frame(n.case=10000,n.control=10000)
sample_sizes[["YRI"]]=data.frame(n.case=2000,n.control=2000)

####method-specific parameters for joint lassosum
GAMMA = c(0.2, 0.5, 0.8)
lambda=exp(seq(log(0.0025), log(0.025), length.out=10)) 
shrink=.9

####parameters used when generating the synthetic data
chrs <- 1:22
TrainTestNFold <- 5
# the first setting
# gammaGenerateData <- 0.8
# lambdaIndexGenerateData <- 5

# the second setting
# gammaGenerateData <- 0.5
# lambdaIndexGenerateData <- 5

# the third setting
# gammaGenerateData <- 0.2
# lambdaIndexGenerateData <- 8

# the forth setting
gammaGenerateData <- 0.8
lambdaIndexGenerateData <- 1

