library(data.table)

#####

source("general_pipeline_parameters.R")

print('the directory of the main pipeline is')
print(main_simulation_pipeline_directory)

print('the directory of the parameter tuning pipeline is')
print(parameter_tuning_pipeline_directory)

print('the gammas under consideration are')
print(GAMMA)

print('the lambdas under consideration are')
print(lambda)

####

ancestries <- c("CEU", "YRI")

####


parameter_tuning_result <- list()
original_testing_result <- list()
best_and_tuned_auc <- list()


for(simulation_index in as.character(800:802)){
  
  jointlassosum_parameter_tuning_directory <- paste0(main_simulation_pipeline_directory, "Work/Sim-",
                                                      simulation_index, "/ParameterTuningData/JointLassosum/")
  
  ####load the synthetic data tuning AUC
  
  
  for(ancestry in ancestries){
    synthetic_auc_one_ancestry <- list()
    for(gamma in GAMMA){
      
      synthetic_auc_file <- paste0(jointlassosum_parameter_tuning_directory,
                                   "JointLassosum-", sprintf("-gamma-%.2f",gamma),"-",
                                   ancestry, "-synthetic-validation-AUC.Rdata")
      synthetic_auc <- get(load(synthetic_auc_file))
      
      ####
      
      synthetic_auc_one_ancestry <- rbind(synthetic_auc_one_ancestry, synthetic_auc)
      
    }
    
    ####
    rownames(synthetic_auc_one_ancestry) <- as.character(GAMMA)
    
    ####
    parameter_tuning_result[[simulation_index]][[ancestry]] <- synthetic_auc_one_ancestry
    
  }
  
  ####FINISH: load the synthetic data tuning AUC
  
  
  ####load the original (true) data testing AUC
  
  
  jointlassosum_directory <- paste0(main_simulation_pipeline_directory, "Work/Sim-",
                                    simulation_index, "/JointLassoSum/")
  
  for(ancestry in ancestries){
    synthetic_auc_one_ancestry <- list()
    for(gamma in GAMMA){
      
      synthetic_auc_file <- paste0(jointlassosum_directory,
                                   "JointLassosum-", sprintf("-gamma-%.2f",gamma),"-",
                                   ancestry, ".TST-AUC.Rdata")
      synthetic_auc <- get(load(synthetic_auc_file))
      
      ####
      
      synthetic_auc_one_ancestry <- rbind(synthetic_auc_one_ancestry, synthetic_auc)
      
    }
    
    ####
    rownames(synthetic_auc_one_ancestry) <- as.character(GAMMA)
    
    ####
    original_testing_result[[simulation_index]][[ancestry]] <- synthetic_auc_one_ancestry
    
  }
  ####FINISH: load the original (true) data testing AUC
  
  
  ####best AUC versus the selected AUC
  
  best_and_tuned_auc_one_setting <- data.frame()
  for(ancestry in ancestries){
    best_tuned_parameter_index <- which.max(parameter_tuning_result[[simulation_index]][[ancestry]])
    best_parameter_index <- which.max(original_testing_result[[simulation_index]][[ancestry]])
    best_and_tuned_auc_one_setting_one_ancestry <- data.frame(ancestry = ancestry,
                                                              oracle_auc = as.numeric(original_testing_result[[simulation_index]][[ancestry]][best_parameter_index]),
                                                              tuned_auc = as.numeric(original_testing_result[[simulation_index]][[ancestry]][best_tuned_parameter_index]))
    best_and_tuned_auc_one_setting <- rbind(best_and_tuned_auc_one_setting,
                                            best_and_tuned_auc_one_setting_one_ancestry)
  }
  best_and_tuned_auc[[simulation_index]] <- best_and_tuned_auc_one_setting
  
  ####FINISH: best AUC versus the selected AUC

}
