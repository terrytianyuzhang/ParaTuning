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

####

# gammaGenerateData <- 0.5
# lambdaIndexGenerateData <- 5

ancestries <- c("CEU", "YRI")

####


parameter_tuning_result <- list()
original_testing_result <- list()
best_and_tuned_auc <- list()


for(simulation_index in as.character(800:809)){
  
  jointlassosum_parameter_tuning_directory <- paste0(main_simulation_pipeline_directory, "Work/Sim-",
                                                      simulation_index, "/ParameterTuningData",
                                                     "_gamma_", sprintf("%.2f",gammaGenerateData), 
                                                     "_lambda_", sprintf("%.4f",lambda[lambdaIndexGenerateData]),
                                                     "/JointLassosum/")
  
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
    
    determine_best_parameter <- function(testing_auc_table,
                                         rule = "99%"){
      
      number_candidate_gamma <- nrow(testing_auc_table)

      favorable_auc_each_gamma <- rep(NA, number_candidate_gamma)
      
      for(gamma_index in 1:number_candidate_gamma){
        
        auc_for_one_gamme <- unlist(testing_auc_table[gamma_index,])
        
        best_auc <- max(auc_for_one_gamme)
        favorable_auc <- 0.99 * best_auc
        
        favorable_lambda_index <- which(auc_for_one_gamme >= favorable_auc)
        favorable_lambda_index_most_regularized <- max(favorable_lambda_index)
        
        favorable_auc_each_gamma[gamma_index] <- auc_for_one_gamme[favorable_lambda_index_most_regularized]
      }
      
      favorable_auc <- max(favorable_auc_each_gamma)
      
      favorable_combination_index <- which(testing_auc_table == favorable_auc)
      
      return(favorable_combination_index)
    }
    
    best_tuned_parameter_index <- determine_best_parameter(parameter_tuning_result[[simulation_index]][[ancestry]])
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

#########SAVE THE RESULTS

parameter_tuning_result_directory <- paste0(parameter_tuning_pipeline_directory,
                                            'Data/')
dir.create(parameter_tuning_result_directory,
           showWarnings = F,recursive = T)

#####

parameter_tuning_result_file <- paste0(parameter_tuning_result_directory, "synthetic_data_tuning_auc",
                                       "_gamma_", sprintf("%.2f",gammaGenerateData), 
                                       "_lambda_", sprintf("%.4f",lambda[lambdaIndexGenerateData]),
                                       ".RData")

save(best_and_tuned_auc, file = parameter_tuning_result_file)


##

parameter_tuning_result_complete_file <- paste0(parameter_tuning_result_directory, "synthetic_data_tuning_auc_complete",
                                                "_gamma_", sprintf("%.2f",gammaGenerateData), 
                                                "_lambda_", sprintf("%.4f",lambda[lambdaIndexGenerateData]),
                                                ".RData")

save(parameter_tuning_result, file = parameter_tuning_result_complete_file)
#########FINISH:SAVE THE RESULTS




