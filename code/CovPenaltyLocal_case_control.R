####locally verify if the covariance penalization works 
####i will normalize everything just as in lassosum
####it seems like my implementation in the previous versions is not strictly following the proposal in the paper
if(1){
#####generate the original training data #####
library(MASS)
library(glmnet)
normalize.mine <- function(x){
  return((x - mean(x))/(sqrt(sum((x - mean(x))^2))))
}
# set.seed(2019)
n <- 200
s <- 5
p <- 700
prev. <- 0.3
beta <- c(100, 0.5, 0.1, -0.5, -0.1)
# x <- mvrnorm(n, mu = rep(0, s), Sigma = diag(1,s))

# risk.score <- x %*% beta
# risk.score <- exp(risk.score)/(1 + exp(risk.score))
# risk.score <- (risk.score - min(risk.score))/ (max(risk.score) - min(risk.score))

gen.casecontrol <- function(targetn, s, p, prev., beta){
  ##generate case-control data set
  # targetn <- 100
  # 
  # s <- 5
  # p <- 700
  # prev. <- 0.05
  
  n <- ceiling(1.2*(targetn*0.5)/prev.)
  beta <- c(1, 0.5, 0.1, -0.5, -0.1)
  x <- mvrnorm(n, mu = rep(0, s), Sigma = diag(1,s))
  
  risk.score <- x %*% beta
  risk.score <- (risk.score - min(risk.score))/ (max(risk.score) - min(risk.score))
  risk.score <- risk.score/(1/(2*prev.))
  y <- rbinom(n, size = 1, prob = risk.score)
  
  ###select a subset of case and control
  control.index <- (1:n)[y == 0]
  control.choose <- sample(control.index, size = targetn/2)
  
  case.index <- (1:n)[y == 1]
  case.choose <- sample(case.index, size = targetn/2)
  
  covariate.matrix <- rbind(x[control.choose,], x[case.choose,])
  covariate.matrix <- cbind(covariate.matrix, 
                            mvrnorm(targetn, mu = rep(0, p - s), Sigma = diag(1,p - s)))
  ###
  y <- c(rep(0,targetn/2), rep(1, targetn/2))
  return(list(x = covariate.matrix,
         y = y))
}

train.data <- gen.casecontrol(targetn = n, s = s,
                p = p, prev. = prev.,
                beta = beta)
x <- train.data$x
y <- train.data$y

####normalize x, y###

for(i in 1:p){
  x[,i] <- normalize.mine(x[,i])
}
y <- normalize.mine(y)
###generate y
# y <- rbinom(n, size = 1, prob = risk.score)
# summary(y)

###attach useless x

# x <- cbind(x, mvrnorm(n, mu = rep(0, p - s), Sigma = diag(1,p - s)))
#############################################

#####generate the testing data #####
testn <- 5e3
test.data <- gen.casecontrol(targetn = testn, s = s,
                              p = p, prev. = prev.,
                              beta = beta)
testx <- test.data$x
testy <- test.data$y

for(i in 1:p){
  testx[,i] <- normalize.mine(testx[,i])
}
# testx <- mvrnorm(testn, mu = rep(0, s), Sigma = diag(1,s))
# 
# risk.score.test <- testx %*% beta
# risk.score.test <- exp(risk.score.test)/(1 + exp(risk.score.test))
# risk.score.test <- (risk.score.test - min(risk.score.test))/ (max(risk.score.test) - min(risk.score.test))
# 
# ###generate y
# testy <- rbinom(testn, size = 1, prob = risk.score.test)
# summary(testy)
# 
# ###attach useless x
# testx <- cbind(testx, mvrnorm(testn, mu = rep(0, p - s), Sigma = diag(1,p - s)))
# #############################################



#######train the model with glmnet###########


glmnet.fit <- glmnet(x = x, y = y)
summary(glmnet.fit)
dim(glmnet.fit$beta)
glmnet.fit$beta[1:10,1:5]

##beta0 corresponds to a smallish lambda
beta0 <- glmnet.fit$beta[, ceiling(NCOL(glmnet.fit$beta)/2)]
lambda.org <- glmnet.fit$lambda

###for each beta_lambda, calculate training error
trn.est <- rep(0, length(lambda.org))
for(i in 1:length(lambda.org)){
  temp.mse <- sum((y - x %*% glmnet.fit$beta[,i])^2)
  trn.est[i] <- temp.mse
}
#############################################

#######generate boostrap data################

###size of reference panel
refn <- n
newx <- mvrnorm(refn, mu = rep(0, p), Sigma = diag(1,p))

for(i in 1:p){
  newx[,i] <- normalize.mine(newx[,i])
}


bootn <- refn
##there are repeated covariates
if(refn < bootn){
  newx <- newx[sample(1:refn, size = bootn, replace = T),]
}else{
  newx <- newx[sample(1:refn, size = bootn, replace = F),]
  
  }

risk.score.new <- newx %*% beta0
risk.score.new <- (risk.score.new - min(risk.score.new))/ (max(risk.score.new) - min(risk.score.new))
# risk.score.new <- 0.05 + (risk.score.new - 0.05)/20 ####case-control

####
B <- 5
boot.Y <- matrix(0, nrow = bootn, ncol = B)
for(b in 1:B){
  
  ###generate y
  newy <- rbinom(bootn, size = 1, prob = risk.score.new)

  newy <- normalize.mine(newy)
  boot.Y[, b] <- newy
  
}

cov.est <- rep(0, length(lambda.org))


for(i in 1:length(lambda.org)){
  score.est <- matrix(0, nrow = bootn, ncol = B)
  for(b in 1:B){
    
    #####use the original lambda
    # glmnet.boot <- glmnet(x = newx, y = boot.Y[,b], 
    #                       lambda = lambda.org[i])  
    
    ####use the square-root rule
    glmnet.boot <- glmnet(x = newx, y = boot.Y[,b], 
                          lambda = lambda.org[i] * sqrt(n/bootn))  
    
    beta.boot <- glmnet.boot$beta
    
    risk.score.boot <- as.numeric(newx %*% beta.boot)
    
    score.est[, b] <- risk.score.boot
    # temp.cov <- temp.cov + cov(risk.score.boot, boot.Y[,b])
  }
  # temp.cov <- temp.cov/B
  
  temp.cov <- 0
  for(j in 1:bootn){
      ###for each individual, calculate covariance
      temp.cov <- temp.cov + cov(score.est[j,], boot.Y[j,])
  }
  
  
  cov.est[i] <- temp.cov
}
#############################################


#####use AUC to determine which is the best lambda#####
glmnet.cvfit <- cv.glmnet(x = x, y = y,
                          lambda = lambda.org)
print(which.min(glmnet.cvfit$cvm))
cv.index <- which.min(glmnet.cvfit$cvm)
  
glmnet.auc.cvfit <- cv.glmnet(x = x, y = y,
                          type.measure = "auc",
                          family  = "binomial")
print(cv.auc.index <- which.max(glmnet.auc.cvfit$cvm))
glmnet.cvfit$cvm


#############################################


###########combine training error and covaraince est##########
cov.est.norm <- cov.est/cov.est[length(lambda.org)] * trn.est[length(lambda.org)]
print(which.min(2*cov.est + trn.est))
cov.index <- which.min(2*cov.est + trn.est)
#############################################

#############compare the difference###################
library(pROC)

cv.auc.test.score <- testx %*% glmnet.auc.cvfit$glmnet.fit$beta[, cv.auc.index]
cv.test.score <- testx %*% glmnet.fit$beta[,cv.index]
cov.test.score <- testx %*% glmnet.fit$beta[,cov.index]

auc(testy, as.numeric(cv.auc.test.score))
auc(testy, as.numeric(cv.test.score))
auc(testy, as.numeric(cov.test.score))

plot(cov.est)
plot(trn.est)
plot(2*cov.est + trn.est)
#############################################
}
#############################################
#############################################
#############################################
#############################################



#################simulation ##############
repeats <- 30
results <- data.frame()
for(rep.i in 1:repeats){
  #####generate the original training data #####
  library(MASS)
  # set.seed(2019)
  n <- 500
  s <- 5
  x <- mvrnorm(n, mu = rep(0, s), Sigma = diag(1,s))
  beta <- c(1, 0.5, 0.1, -0.5, -0.1)
  
  risk.score <- x %*% beta
  risk.score <- exp(risk.score)/(1 + exp(risk.score))
  risk.score <- (risk.score - min(risk.score))/ (max(risk.score) - min(risk.score))
  
  ###generate y
  y <- rbinom(n, size = 1, prob = risk.score)
  summary(y)
  
  ###attach useless x
  p <- 10
  x <- cbind(x, mvrnorm(n, mu = rep(0, p - s), Sigma = diag(1,p - s)))
  #############################################
  
  #####generate the testing data #####
  testn <- 5e3
  testx <- mvrnorm(testn, mu = rep(0, s), Sigma = diag(1,s))
  
  risk.score.test <- testx %*% beta
  risk.score.test <- exp(risk.score.test)/(1 + exp(risk.score.test))
  risk.score.test <- (risk.score.test - min(risk.score.test))/ (max(risk.score.test) - min(risk.score.test))
  
  ###generate y
  testy <- rbinom(testn, size = 1, prob = risk.score.test)
  summary(testy)
  
  ###attach useless x
  testx <- cbind(testx, mvrnorm(testn, mu = rep(0, p - s), Sigma = diag(1,p - s)))
  #############################################
  
  
  
  #######train the model with glmnet###########
  library(glmnet)
  
  normalize.mine <- function(x){
    return((x - mean(x))/(sqrt(sum((x - mean(x))^2))))
  }
  
  ####normalize x, y###
  y <- normalize.mine(y)
  
  for(i in 1:p){
    x[,i] <- normalize.mine(x[,i])
  }
  
  glmnet.fit <- glmnet(x = x, y = y)
  summary(glmnet.fit)
  dim(glmnet.fit$beta)
  head(glmnet.fit$beta)
  
  ##beta0 corresponds to a smallish lambda
  beta0 <- glmnet.fit$beta[, ceiling(NCOL(glmnet.fit$beta)/2)]
  lambda.org <- glmnet.fit$lambda
  
  ###for each beta_lambda, calculate training error
  trn.est <- rep(0, length(lambda.org))
  for(i in 1:length(lambda.org)){
    temp.mse <- sum((y - x %*% glmnet.fit$beta[,i])^2)
    trn.est[i] <- temp.mse
  }
  #############################################
  
  #######generate boostrap data################
  
  ###size of reference panel
  refn <- n
  newx <- mvrnorm(refn, mu = rep(0, p), Sigma = diag(1,p))
  
  for(i in 1:p){
    newx[,i] <- normalize.mine(newx[,i])
  }
  
  
  bootn <- refn
  ##there are repeated covariates
  if(refn < bootn){
    newx <- newx[sample(1:refn, size = bootn, replace = T),]
  }else{
    newx <- newx[sample(1:refn, size = bootn, replace = F),]
    
  }
  
  risk.score.new <- newx %*% beta0
  risk.score.new <- (risk.score.new - min(risk.score.new))/ (max(risk.score.new) - min(risk.score.new))
  
  ####
  B <- 5
  boot.Y <- matrix(0, nrow = bootn, ncol = B)
  for(b in 1:B){
    
    ###generate y
    newy <- rbinom(bootn, size = 1, prob = risk.score.new)
    
    newy <- normalize.mine(newy)
    boot.Y[, b] <- newy
    
  }
  
  cov.est <- rep(0, length(lambda.org))
  
  
  for(i in 1:length(lambda.org)){
    score.est <- matrix(0, nrow = bootn, ncol = B)
    for(b in 1:B){
      
      #####use the original lambda
      # glmnet.boot <- glmnet(x = newx, y = boot.Y[,b], 
      #                       lambda = lambda.org[i])  
      
      ####use the square-root rule
      glmnet.boot <- glmnet(x = newx, y = boot.Y[,b], 
                            lambda = lambda.org[i] * sqrt(n/bootn))  
      
      beta.boot <- glmnet.boot$beta
      
      risk.score.boot <- as.numeric(newx %*% beta.boot)
      
      score.est[, b] <- risk.score.boot
      # temp.cov <- temp.cov + cov(risk.score.boot, boot.Y[,b])
    }
    # temp.cov <- temp.cov/B
    
    temp.cov <- 0
    for(j in 1:bootn){
      ###for each individual, calculate covariance
      temp.cov <- temp.cov + cov(score.est[j,], boot.Y[j,])
    }
    
    
    cov.est[i] <- temp.cov
  }
  #############################################
  
  
  #####use AUC to determine which is the best lambda#####
  glmnet.cvfit <- cv.glmnet(x = x, y = y,
                            lambda = lambda.org)
  print(which.min(glmnet.cvfit$cvm))
  cv.index <- which.min(glmnet.cvfit$cvm)
  
  # glmnet.cvfit <- cv.glmnet(x = x, y = y,
  #                           lambda = lambda.org,
  #                           type.measure = "auc",
  #                           family  = "binomial")
  # print(which.max(glmnet.cvfit$cvm))
  # glmnet.cvfit$cvm
  
  
  #############################################
  
  
  ###########combine training error and covaraince est##########
  cov.est.norm <- cov.est/cov.est[length(lambda.org)] * trn.est[length(lambda.org)]
  print(which.min(2*cov.est + trn.est))
  cov.index <- which.min(2*cov.est + trn.est)
  #############################################
  
  #############compare the difference###################
  library(pROC)
  
  cv.test.score <- testx %*% glmnet.fit$beta[,cv.index]
  cov.test.score <- testx %*% glmnet.fit$beta[,cov.index]
  
  cv.auc <- auc(testy, as.numeric(cv.test.score))
  cov.auc <- auc(testy, as.numeric(cov.test.score))
  
  plot(cov.est)
  plot(trn.est)
  plot(2*cov.est + trn.est)
  #############################################
  
temp.results <- data.frame(rep.i = rep.i,
                           cv.auc = cv.auc,
                           cov.auc = cov.auc,
                           n = n,
                           p = p,
                           s = s,
                           refn = refn)
results <- rbind(results, temp.results)
}

save(results, 
     file = '/Users/tianyu/Documents/ParaTuning/data/cov_penalty_local_study_times2_normalize_paper_3.RData')

results <- data.table(results)
results[, diff := cov.auc - cv.auc]
ggplot(results) +
  geom_histogram(aes(x=cv.auc), 
                 fill="blue", alpha=0.8, binwidth = 0.01)+
  geom_histogram(aes(x=cov.auc), 
                 fill="red", alpha=0.8, binwidth = 0.01)+
  ggtitle(paste0("n=",results$n[1]," p=", results$p[1]))+
  xlab('testing auc')

ggplot(results) +
  geom_histogram(aes(x=diff), 
                 fill="blue", alpha=0.8, binwidth = 0.01)+
  ggtitle(paste0("n=",results$n[1]," p=", results$p[1]))+
  xlab('testing auc diff')


########################################
########################################
########################################

#######load existing data and plot#####
for(i in 1:6){
results <- get(load(paste0("/Users/tianyu/Documents/ParaTuning/data/cov_penalty_local_study_times2_normalize_",i,".RData")))
print(results$n)
print(results$p)

results <- data.table(results)
results[, diff := cov.auc - cv.auc]
print(ggplot(results) +
  geom_histogram(aes(x=cv.auc), 
                 fill="blue", alpha=0.8, binwidth = 0.01)+
  geom_histogram(aes(x=cov.auc), 
                 fill="red", alpha=0.8, binwidth = 0.01)+
  ggtitle(paste0("n=",results$n[1]," p=", results$p[1]))+
    xlab('testing auc')
  )

print(ggplot(results) +
  geom_histogram(aes(x=diff), 
                 fill="blue", alpha=0.8, binwidth = 0.01)+
    ggtitle(paste0("n=",results$n[1]," p=", results$p[1]))+
    xlab('testing auc: Cov pen. minus cv'))
}


