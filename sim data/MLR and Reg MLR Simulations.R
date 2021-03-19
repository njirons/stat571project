## This file is to run simulations and microbiome application for
## the MLR, MLR + Lasso, MLR +  Ridge, MLR+ Electric Net methods
## This file requires: xtrain#.RData, ytrain#.RData,
##                     beta#_new.RData
##                    For Scenario # = 1, 2, 3
##                    xtrain.RData, ytrain.RData
##                    xtest.RData, ytest.RData
##                    For microbiome data
## This file outputs: scen#list.RData for scenarios 1,2,3 
##                    list of data.frames with the 4 performance metrics
##                    for the MLR, MLR+Lasso,Ridge,EN methods

library(glmnet)
library("MGLM")
library(purrr)
library(nnet)
library(mosaic)
library(dplyr)

##

dir_link <- "/Users/Leah/Desktop/UW Stuff/Biost 571/Biost 571 Group Project/Code and Data/Sim Datasets/"



# file_string <- c(sprintf("x%s_new.RData",1:3),
#                  sprintf("y%s_new.RData",1:3),
#                  sprintf("beta%s_new.RData",1:3))
file_string <- c(sprintf("xtrain%s.RData",1:3),
                 sprintf("ytrain%s.RData",1:3),
                 sprintf("xtest%s.RData",1:3),
                 sprintf("ytest%s.RData",1:3),
                 sprintf("beta%s_new.RData",1:3))


## Load simulation datasets
for(i in c(file_string)){
  load(paste0(dir_link,i))
}

set.seed(2021)



## MLR (using nnet)
sim_nnet <- function(x.train, x.test, y.train, y.test, betamat,nscene){
  ncat <- ncol(y.train)
  beta <- betamat
  
  # make last response category the reference
  # fit w/o intercept cuz already included in x.train
  # Keep maxit large so model converges
  plain.mod<-plain.mod<-multinom(y.train[,c(ncat,1:(ncat-1))]~x.train-1, MaxNWts=1200, maxit=1000)
  
  # match names of test dataset to model
  colnames(x.test) <- colnames(coef(plain.mod))
  # predict.nnet requires test dataset same size, will remove extra rows later
  plain.prob0<-predict(plain.mod, rbind(x.test, x.test, x.test, x.test,x.test, x.test)[1:nrow(x.train),], "probs") #
  plain.prob<-plain.prob0[1:nrow(x.test),c(2:ncat,1)]
  #plain.prob<-cbind(plain.prob0[1:nrow(x.test),-1], plain.prob0[1:nrow(x.test),1])
  yhats.plain<- rowSums(y.test)*plain.prob
  
  data.frame(dist = nscene,
             sample = nrow(x.train),
             response = ncat,
             method = "MLRnnet",
             MSE= mean((t(coef(plain.mod))-beta)^2),
             MAE= mean(abs(t(coef(plain.mod))-beta)),
             MSPE=mean((y.test-yhats.plain)^2),
             MAPE = mean(abs(y.test-yhats.plain)))
}

testnnet<-sim_nnet(xtrain1[[d]][[j]],xtest1[[d]][[j]], 
               ytrain1[[d]][[j]],ytest1[[d]][[j]], 
               beta1[[d]][[j]],1) 

## Penalized MLR (lasso, ridge, EN)  with glmnet
sim_metric <- function(x.train, x.test, y.train, y.test, betamat,nscene){
  ncat <- ncol(y.train)
  beta <- betamat
  alpha_vec <- c(0, 0.5, 1)
  metric_mat <- matrix(ncol=4, nrow=3)
  kfold<-5
  
  ## MLR + Regularization
  for(i in 1:3){
    cv.mod<-cv.glmnet(x.train[,-1], # x and y must be matrices (exclude intercept)
                      ## if too many response variables, model doesn't converge
                      y.train, 
                      family = "multinomial",
                      ## lasso: alpha=1
                      ## elastic net: 0< alpha < 1
                      ## ridge: alpha= 0
                      alpha = alpha_vec[i],
                      type.measure= "mse", 
                      intercept=TRUE,
                      nfold= kfold,
                      grouped=FALSE,
                      maxit= 1e+05,
                      thresh=1e-5,
                      nlambda=50) 
    coef.mod<- coef(cv.mod, s="lambda.min")
    coef.mat0 <- sapply(coef.mod, as.matrix)
    coef.mat <- (coef.mat0-coef.mat0[,ncat])[,-ncat]
    pred.mod<- predict(cv.mod, 
                       ## newx is the test data (even tho I used the 5 rows from training data)
                       newx=x.test[,-1],
                       ## link outputs yhats
                       type="response",
                       ## s determines which lambda to use
                       ## could use lambda.min, lambda.1se, or random lambda choice
                       s="lambda.min") 
    yhats<- rowSums(y.test)*pred.mod[,,1]
    
    # mse, mae, mspe, mape
    metric_mat[i,] <- c(mean((coef.mat-beta)^2), 
                          mean(abs(coef.mat-beta)), 
                          mean((y.test-yhats)^2), 
                          mean(abs(y.test-yhats)))
  }

  newdf<-cbind(data.frame(dist = nscene,
                          sample = nrow(x.train),
                          response = ncat,
                          method = c("Ridge","EN","Lasso")), metric_mat)
  names(newdf)[5:8] <- c("MSE", "MAE","MSPE","MAPE") #("MSE" , "MAE" , "MSPE" , "MAPE")
  newdf

}


####### Running Things ##########3

### Scenario 1, 9 different combinations, 50 datasets each
## Running Scenario 1 MLR w/ nnet  (increased weights), takes ~ 30 minutes
tstart <- Sys.time()
scen1nnet_list<- lapply(1:9,   function(d) map_dfr(1:50, function(j) 
  sim_nnet(xtrain1[[d]][[j]],xtest1[[d]][[j]], 
           ytrain1[[d]][[j]],ytest1[[d]][[j]], 
           beta1[[d]][[j]],1)) )
tend <- Sys.time()
tend-tstart  

## Running Scenario 1 MLR + Lasso/Ridge/EN, takes ~ 42 minutes
tstart<-Sys.time()
scen1glmnet_list<- lapply(1:9,   function(d) map_dfr(1:50, function(j) 
  sim_metric(xtrain1[[d]][[j]],xtest1[[d]][[j]], 
             ytrain1[[d]][[j]],ytest1[[d]][[j]], 
             beta1[[d]][[j]],1)) )
tend <- Sys.time()
tend-tstart

## Combining
scen1list <- lapply(1:9, function(i) rbind(scen1glmnet_list[[i]],scen1nnet_list[[i]])%>% arrange(method))
#save(scen1list, file="/Users/Leah/Desktop/UW Stuff/Biost 571/Biost 571 Group Project/Code and Data/Sim Datasets/scen1list.RData")

### Scenario 2, 3 different combinations, 50 datasets each


## MLR ~ 10 min
rho_vec <- c(0.1, 0.5, 0.9)
scen2nnet_list<- lapply(1:3,   function(d) cbind(data.frame(rho=rho_vec[d]),
                                                 map_dfr(1:50, function(j) sim_nnet(xtrain2[[d]][[j]],xtest2[[d]][[j]], 
                                                                                    ytrain2[[d]][[j]],ytest2[[d]][[j]], 
                                                                                    beta2[[d]][[j]],2)) ))

## MLR + Lasso/RR/EN, ~ 15 min
tstart<-Sys.time()
scen2glmnet_list<- lapply(1:3,   function(d) cbind(data.frame(rho=rho_vec[d]),
                                                   map_dfr(1:50, function(j) sim_metric(xtrain2[[d]][[j]],
                                                                                        xtest2[[d]][[j]], 
                                                                                        ytrain2[[d]][[j]],
                                                                                        ytest2[[d]][[j]], 
                                                                                        beta2[[d]][[j]],2)) ))

scen2list <- lapply(1:3, function(i) rbind(scen2glmnet_list[[i]],scen2nnet_list[[i]])%>% arrange(method))
tend <- Sys.time()
tend-tstart

#save(scen2list, file="/Users/Leah/Desktop/UW Stuff/Biost 571/Biost 571 Group Project/Code and Data/Sim Datasets/scen2list.RData")





### Scenario 3, 9 different combinations, 50 datasets each

## Running Scenario 3 MLR w/ nnet  (increased weights), ~ 35 min
tstart <- Sys.time()
scen3nnet_list<- lapply(1:9,   function(d) map_dfr(1:50, function(j) 
  sim_nnet(xtrain3[[d]][[j]],xtest3[[d]][[j]], 
           ytrain3[[d]][[j]],ytest3[[d]][[j]], 
           beta3[[d]][[j]],3)) )
tend <- Sys.time()
tend-tstart  

## MLR + Lasso/RR/EN, ~42 minutes
tstart<-Sys.time()

scen3glmnet_list<- lapply(1:9,   function(d) map_dfr(1:50, function(j) 
  sim_metric(xtrain3[[d]][[j]],xtest3[[d]][[j]], 
             ytrain3[[d]][[j]],ytest3[[d]][[j]], 
             beta3[[d]][[j]],3)) )
tend <- Sys.time()
tend-tstart

scen3list <- lapply(1:9, function(i) rbind(scen3glmnet_list[[i]],scen3nnet_list[[i]])%>% arrange(method))
#save(scen3list, file="/Users/Leah/Desktop/UW Stuff/Biost 571/Biost 571 Group Project/Code and Data/Sim Datasets/scen3list.RData")



########### Microbiome Application #######

dir_link2 <- "/Users/Leah/Desktop/UW Stuff/Biost 571/Biost 571 Group Project/Code and Data/microbiome data/"


file_string2 <- c(paste0("x",c("train.RData","test.RData")),
                  paste0("y",c("train.RData","test.RData")))


## Load simulation datasets
for(i in file_string2){
  load(paste0(dir_link2,i))
}

sim_nnet_app <- function(x.train, x.test, y.train, y.test, nscene){
  ncat <- ncol(y.train)
  
  # make last response category the reference
  # fit w/o intercept cuz already included in x.train
  plain.mod<-multinom(y.train[,c(ncat,1:(ncat-1))]~x.train-1, MaxNWts=1200)
  
  # match names of test dataset to model
  colnames(x.test) <- colnames(coef(plain.mod))
  # predict.nnet requires test dataset same size, will remove extra rows later
  plain.prob0<-predict(plain.mod, rbind(x.test, x.test, x.test, x.test,x.test, x.test)[1:nrow(x.train),], "probs") #
  plain.prob<-plain.prob0[1:nrow(x.test),c(2:ncat,1)]
  yhats.plain<- rowSums(y.test)*plain.prob
  
  data.frame(dist = nscene,
             sample = nrow(x.train),
             response = ncat,
             method = "MLRnnet",
             MSPE=mean((y.test-yhats.plain)^2),
             MAPE = mean(abs(y.test-yhats.plain)))
}

sim_nnet_app(xtrain[[j]], xtest[[j]], ytrain[[j]], ytest[[j]], "app")
ncat <- ncol(ytrain[[1]])
multinom(ytrain[[j]][,5:8]~xtrain[[j]]-1, MaxNWts=1300)  #, #MaxNWts=1200)

sim_metric_app <- function(x.train, x.test, y.train, y.test, nscene){
  
  # nobs <- nrow(xmat)
  ncat <- ncol(y.train)
  
  alpha_vec <- c(0, 0.5, 1)
  metric_mat <- matrix(ncol=2, nrow=3)
  kfold<-3
  
  
  
  ## MLR + Regularization
  for(i in 1:3){
    cv.mod<-cv.glmnet(x.train[,-1], # x and y must be matrices (exclude intercept)
                      ## if too many response variables, model doesn't converge
                      y.train, 
                      family = "multinomial",
                      ## lasso: alpha=1
                      ## elastic net: 0< alpha < 1
                      ## ridge: alpha= 0
                      alpha = alpha_vec[i],
                      type.measure= "mse", 
                      intercept=TRUE,
                      nfold= kfold,
                      grouped=FALSE,
                      maxit= 1e+05,
                      thresh=1e-5,
                      nlambda=30) 
    coef.mod<- coef(cv.mod, s="lambda.min")
    coef.mat0 <- sapply(coef.mod, as.matrix)
    coef.mat <- (coef.mat0-coef.mat0[,ncat])[,-ncat]
    pred.mod<- predict(cv.mod, 
                       ## newx is the test data (even tho I used the 5 rows from training data)
                       newx=x.test[,-1],
                       ## link outputs yhats
                       type="response",
                       ## s determines which lambda to use
                       ## could use lambda.min, lambda.1se, or random lambda choice
                       s="lambda.min") 
    yhats<- rowSums(y.test)*pred.mod[,,1]
    
    # mse, mae, mspe, mape
    metric_mat[i,] <- c( mean((y.test-yhats)^2), 
                           mean(abs(y.test-yhats)))
}
  newdf<-cbind(data.frame(dist = nscene,
                          sample = nrow(x.train),
                          response = ncat,
                          method = c("MLRnnet","Ridge","EN","Lasso","MLRglmnet")), metric_mat)
  names(newdf)[5:6] <- c("MSPE","MAPE") #("MSE" , "MAE" , "MSPE" , "MAPE")
  newdf
}

j<-2
sim_metric_app(xtrain[[j]], xtest[[j]], ytrain[[j]], ytest[[j]], "app")
undebug(sim_metric_app)
cv.mod<-cv.glmnet(rbind(xtrain[[j]],xtest[[j]]) [,-1], # x and y must be matrices (exclude intercept)
                  ## if too many response variables, model doesn't converge
                  rbind(ytrain[[j]], ytest[[j]]), 
                  family = "multinomial",
                  ## lasso: alpha=1
                  ## elastic net: 0< alpha < 1
                  ## ridge: alpha= 0
                  alpha = 0,
                  type.measure= "mse", 
                  intercept=TRUE,
                  nfold= 3,
                  grouped=FALSE,
                  maxit= 1e+06,
                  thresh=1e-7,
                  nlambda=30) 





