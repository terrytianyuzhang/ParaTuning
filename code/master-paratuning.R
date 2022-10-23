print('original model fitting')
source("tianyu-combined-lassosum-training.R")

print('generate boot')
source("generate-boot-data-calib.R")

print('fit boostrap models')
source("combined-lassosum-boot-training.R")

