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
ld.populations=c("EUR.hg38","AFR.hg38"); names(ld.populations)=c("CEU","YRI")

####
# information on the number of samples
sample_sizes <- list()
sample_sizes[["CEU"]]=data.frame(n.case=10000,n.control=10000)
sample_sizes[["YRI"]]=data.frame(n.case=2000,n.control=2000)