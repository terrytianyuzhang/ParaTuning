####locally verify if the covariance penalization works 

#####generate the original training data #####
library(MASS)
# set.seed(2019)
n <- 100
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
p <- 300
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
glmnet.fit <- glmnet(x = x, y = y)
summary(glmnet.fit)
dim(glmnet.fit$beta)
head(glmnet.fit$beta)

##beta0 corresponds to a smallish lambda
beta0 <- glmnet.fit$beta[, NCOL(glmnet.fit$beta) - 3]
lambda.org <- glmnet.fit$lambda

###for each beta_lambda, calculate training error
trn.est <- rep(0, length(lambda.org))
for(i in 1:length(lambda.org)){
  temp.mse <- mean((y - x %*% glmnet.fit$beta[,i])^2)
  trn.est[i] <- temp.mse
}
#############################################

#######generate boostrap data################

###size of reference panel
refn <- n/2
newx <- mvrnorm(refn, mu = rep(0, p), Sigma = diag(1,p))


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
B <- 20
boot.Y <- matrix(0, nrow = bootn, ncol = B)
for(b in 1:B){
  
  ###generate y
  newy <- rbinom(bootn, size = 1, prob = risk.score.new)
  boot.Y[, b] <- newy
  
}

cov.est <- rep(0, length(lambda.org))

for(i in 1:length(lambda.org)){
  temp.cov <- 0
  for(b in 1:B){
    
    #####use the original lambda
    # glmnet.boot <- glmnet(x = newx, y = boot.Y[,b], 
    #                       lambda = lambda.org[i])  
    
    ####use the square-root rule
    glmnet.boot <- glmnet(x = newx, y = boot.Y[,b], 
                          lambda = lambda.org[i] * sqrt(n/bootn))  
    
    beta.boot <- glmnet.boot$beta
    
    risk.score.boot <- as.numeric(newx %*% beta.boot)
    
    temp.cov <- temp.cov + cov(risk.score.boot, boot.Y[,b])
  }
  temp.cov <- temp.cov/B
  
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

auc(testy, cv.test.score)
auc(testy, cov.test.score)

plot(cov.est)
plot(trn.est)
#############################################

#############################################
#############################################
#############################################
#############################################



#################simulation ##############
repeats <- 30
results <- data.frame()
for(rep.i in 1:repeats){
####locally verify if the covariance penalization works 

#####generate the original training data #####
library(MASS)
# set.seed(2019)
n <- 100
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
p <- 700
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
glmnet.fit <- glmnet(x = x, y = y, nlambda = 30)
summary(glmnet.fit)
dim(glmnet.fit$beta)
head(glmnet.fit$beta)

##beta0 corresponds to a smallish lambda
beta0 <- glmnet.fit$beta[, NCOL(glmnet.fit$beta)/2]
lambda.org <- glmnet.fit$lambda

###for each beta_lambda, calculate training error
trn.est <- rep(0, length(lambda.org))
for(i in 1:length(lambda.org)){
  temp.mse <- mean((y - x %*% glmnet.fit$beta[,i])^2)
  trn.est[i] <- temp.mse
}
#############################################

#######generate boostrap data################

###size of reference panel
refn <- n
newx <- mvrnorm(refn, mu = rep(0, p), Sigma = diag(1,p))


bootn <- n
##there are repeated covariates
if(refn < bootn){
  newx <- newx[sample(1:refn, size = bootn, replace = T),]
}else{
  newx <- newx[sample(1:refn, size = bootn, replace = F),]
  
}

risk.score.new <- newx %*% beta0
risk.score.new <- (risk.score.new - min(risk.score.new))/ (max(risk.score.new) - min(risk.score.new))

####
B <- 20
boot.Y <- matrix(0, nrow = bootn, ncol = B)
for(b in 1:B){
  
  ###generate y
  newy <- rbinom(bootn, size = 1, prob = risk.score.new)
  boot.Y[, b] <- newy
  
}

cov.est <- rep(0, length(lambda.org))

for(i in 1:length(lambda.org)){
  temp.cov <- 0
  for(b in 1:B){
    
    #####use the original lambda
    # glmnet.boot <- glmnet(x = newx, y = boot.Y[,b], 
    #                       lambda = lambda.org[i])  
    
    ####use the square-root rule
    glmnet.boot <- glmnet(x = newx, y = boot.Y[,b], 
                          lambda = lambda.org[i] * sqrt(n/bootn))  
    
    beta.boot <- glmnet.boot$beta
    
    risk.score.boot <- as.numeric(newx %*% beta.boot)
    
    temp.cov <- temp.cov + cov(risk.score.boot, boot.Y[,b])
  }
  temp.cov <- temp.cov/B
  
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

cv.auc <- auc(testy, cv.test.score)
cov.auc <- auc(testy, cov.test.score)

plot(cov.est)
plot(trn.est)
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
     file = '/Users/tianyu/Documents/ParaTuning/data/cov_penalty_local_study_times2_5.RData')

results <- data.table(results)
ggplot(results) +
  geom_histogram(aes(x=cv.auc), 
                 fill="blue", alpha=0.8, binwidth = 0.01)+
  geom_histogram(aes(x=cov.auc), 
                 fill="red", alpha=0.8, binwidth = 0.01)
