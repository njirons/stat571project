## Fitting and predicting MLR and MLR+Lasso/Ridge/Elastic Net

#setwd("/Users/Leah/Desktop/UW Stuff/Biost 571/Biost 571 Group Project/Code and Data")
#setwd("~/Documents/Classes/Winter 2021/stat 571/final")
library(readr)
library(tidyr)
library(dplyr)
library(nnet)
library(mlogit)
library(reshape2)
library(foreign)
library(stats4)
library(glmnet)
library(nnet)

######### Loading Microbiome Data #######

covariates <- read_csv("GvHD_Covariates_571.csv")
# covariates <- na.omit(covariates)
# View(covariates)
# glimpse(covariates)

taxon <- read.csv("GvHD_Microbiome_Data_571.csv")

# taxon %>% group_by(patientID) %>% filter(sample_day==0)
# taxon <- taxon[,-(1:5)]


taxon <- taxon %>% group_by(patientID) %>% 
  summarise(across(names(taxon)[-(1:5)], sum))
taxon <- taxon[,-1]

m <- 138   # Using all people with complete.case covariates 
J <- 500
x <- na.omit(covariates %>% select(sex, agvhgrd,agvhskn,agvhlvr,agvhgut,status,
                                   donsex,txage,donage,race,donrel))
subjects <- sample(1:dim(x)[1],size=m, replace=FALSE)

x$status <- x$status == "Remission"
x$donage <- x$donage - mean(x$donage)
x$txage <- x$txage - mean(x$txage)
x$sex <- x$sex == "Female"
x$donsex <- x$donsex == "Female"
x$race <- x$race == "Caucasian"
x$donrel <- x$donrel == "Not Related"
# x$agvhgrd[is.na(x$agvhgrd)] <- 0
# x[is.na(x)] <- 0
x$intercept <- 1
K <- ncol(x)

x <- x[subjects,]
taxa <- taxon[subjects,1:J]

y <- as.matrix(taxa)
x <- as.matrix(x)

############ MLR (slight work in progress) ####################
ntaxa0 <- 70

fit.plain = multinom(y[,1:ntaxa0] ~ x)
# Comments:
# The function below cannot handle linear separabilty either.
# It works well with y[,1:70], but fails at y[,1:80] or above.

dfx <- as.data.frame(x)
names(dfx) <- colnames(coef(fit.plain))[-1]
pred.plain <- predict(fit.plain, newdata = dfx, type = "probs")

# Comments:
# The package nnet requires that we input a data frame into the
# function predict(). Hence we create a data frame called dfx
# Note: for now we set our test set as the training set dfx, but eventually we should
# create a new dataset as the test set
# For some reason, this function requires the same number of rows for its 
# newdata (the test data) as the training dataset (so it works when that's true)
# Or there's some glitch but we've a lotta time on it and can't figure it out
## Function might need column names or need to be fit as a dataset? 

yhats.plain<- rowSums(y[,1:ntaxa0])*pred.plain

mse.yhat.plain <- mean((y[,1:ntaxa0]-yhats.plain)^2)

mse.prob.plain<- mean((y[,1:ntaxa0]/rowSums(y[,1:ntaxa0])-pred.plain)^2)


## Comments:
## glmnet can predict probabilities for each response variable (the y hats were a little suspect)
## To compare true y counts to the predicted y probabilities we can either divide true y counts by rowsum
## or multiply predicted y probabilities by true y row sums
## Used both methods to get mses based on yhat and probability

############# MLR w/ Regularization ###############

## Choose k folds for selecting tuning parameter
## Since our entropy method is doing hold out cross validation,
## should we try to choose a lower number of folds so that the two methods
## are more comparable?
kfold <- 5

ntaxa <- 150

#### Lasso

# Fitting 100 models to find the optimal tuning parameter lambda
cv.lasso<-cv.glmnet(x, # x and y must be matrices
                    ## if too many response variables, model doesn't converge
                    y[,1:ntaxa], 
                    family = "multinomial",
                    ## alpha=1 means lasso
                    alpha = 1,
                    type.measure= "mse",  
                    nfold= kfold,
                    grouped=FALSE)  #lambda=exp(seq(-2,100,length=100) )
plot(cv.lasso)               
## Comments:
## Function chooses 100 lambdas to perform cross validation from (w/ MSE metric)
## You can hand specify lambdas if you want
## For our microbiome data, we get a monotonic relationship with lambda and MSE
## eg lambda = infinity minimizes MSE
## Convergence of model is sensitive number of y variables and to number of folds (kinda)
## All models can handle the first 150 taxa from microbiome dataset
## Should probably sample the taxa when perform actual data analysis
## With fewer folds in cv, more response variables can lead to convergence
## Since our entropy method is using hold out cv, should we use fewer folds in these models for consistency?
## Can do parallel computing with cv.glmnet, but function doesn't take long enough to require it
## I experimented with parallel, but sometimes models didn't converge when I used PC, 
##    even tho they converged when I didn't use parallel computing
## PC likely more useful when we have more covariates
## glmnet automatically rescales covariates before running and then returns them to original scale


## list of sparse matrices for coefficient of each taxa. "." means 0.
coef.lasso<- coef(cv.lasso)

## Comments:
## For lasso, best model results in only an intercept coefficient

## Predicting y hats and probabilities using model with best tuning parameter lambda
pred.lasso<- predict(cv.lasso, 
                     ## newx is the test data (even tho I used the 5 rows from training data)
                     newx=x[1:5,],
                     ## link outputs yhats
                     type="response",
                     ## s determines which lambda to use
                     ## could use lambda.min, lambda.1se, or random lambda choice
                     s="lambda.min") 
## Comments: 
## pred.lasso reports an array, use pred.lasso[,,1] to get yhats as a matrix

yhats.lasso<- rowSums(y[1:5,1:ntaxa])*pred.lasso[,,1]

mse.yhat.lasso <- mean((y[1:5,1:ntaxa]-yhats.lasso)^2)

mse.prob.lasso<- mean((y[1:5,1:ntaxa]/rowSums(y[1:5,1:ntaxa])-pred.lasso[,,1])^2)


## Comments:
## glmnet can predict probabilities for each response variable (the y hats were a little suspect)
## To compare true y counts to the predicted y probabilities we can either divide true y counts by rowsum
## or multiply predicted y probabilities by true y row sums
## Used both methods to get mses based on yhat and probability

#### Ridge


# Fitting 100 models to find the optimal tuning parameter lambda
cv.ridge<-cv.glmnet(x, y[,1:ntaxa], 
                    family = "multinomial",
                    ## alpha = 0 is for Ridge
                    alpha = 0, 
                    type.measure="mse",
                    nfold= kfold,
                    grouped=FALSE) #lambda=exp(seq(-2,100,length=100) )
plot(cv.ridge)
## Comments:
## See lasso comments

## list of sparse matrices for coefficient of each taxa. "." means 0.
coef.ridge<-coef(cv.ridge)

## Predicting y hats and probabilities using model with best tuning parameter lambda
pred.ridge<- predict(cv.ridge, 
                     newx=x[1:5,],
                     type="response", 
                     s="lambda.min")

yhats.ridge<- rowSums(y[1:5,1:ntaxa])*pred.ridge[,,1]

mse.yhat.ridge <- mean((y[1:5,1:ntaxa]-yhats.ridge)^2)

mse.prob.ridge<- mean((y[1:5,1:ntaxa]/rowSums(y[1:5,1:ntaxa])-pred.ridge[,,1])^2)


## Comments:
## glmnet can predict probabilities for each response variable (the y hats were a little suspect)
## To compare true y counts to the predicted y probabilities we can either divide true y counts by rowsum
## or multiply predicted y probabilities by true y row sums
## Used both methods to get mses based on yhat and probability




### Elastic Net

# Fitting 100 models to find the optimal tuning parameter lambda
cv.elast<-cv.glmnet(x, y[,1:ntaxa], 
                    family = "multinomial", 
                    ## 0 < alpha < 1 
                    alpha = 0.5, 
                    type.measure="mse",
                    nfold=kfold,
                    grouped=FALSE)   #lambda=exp(seq(-2,100,length=100) )
plot(cv.elast)
## Comments:
## Elastic net is a combination of Ridge and Lasso, so choose alpha to indicate how similar
##    you want the penalization to be like Ridge (alpha=0) or Lasso (alpha=1)
## Maybe we want to fit multiple elastic nets, Idk
## See lasso comments

## list of sparse matrices for coefficient of each taxa. "." means 0.
coef.elast<-coef(cv.elast)

## Predicting y hats and probabilities using model with best tuning parameter lambda
pred.elast<- predict(cv.elast, newx=x[1:5,],
                     type="response", s="lambda.min")
yhats.elast<- rowSums(y[1:5,1:ntaxa])*pred.elast[,,1]

mse.yhat.elast <- mean((y[1:5,1:ntaxa]-yhats.elast)^2)

mse.prob.elast<- mean((y[1:5,1:ntaxa]/rowSums(y[1:5,1:ntaxa])-pred.elast[,,1])^2)

## Comments:
## glmnet can predict probabilities for each response variable (the y hats were a little suspect)
## To compare true y counts to the predicted y probabilities we can either divide true y counts by rowsum
## or multiply predicted y probabilities by true y row sums
## Used both methods to get mses based on yhat and probability


## Comparing MSEs of yhat and prob 
## (but plain MLR probably has different number of response variables 
##  than regularization methods)
## mses are different when number of response variables are low
## when number of response variables are high, the lambda that optimizes
## the model shrinks all the covariates to 0 except the intercept
## so all the mses are the same
(mse.yhat<- c(mse.yhat.plain,mse.yhat.lasso,mse.yhat.elast,mse.yhat.ridge))
(mse.prob<- c(mse.prob.plain,mse.prob.lasso,mse.prob.elast,mse.prob.ridge))


## experimented with parallel computing the cv.glmnet function
##    (would just add parallel=TRUE as an option), 
##    but sometimes model wouldn't converge with parallel computing when it did without it
##    Also could be due to the high dimensionality of the data. 
##    Maybe with ~100 response variables, it would work fine.
## Didn't seem to speed up function that much

## Run this first to establish parallel computing
library(doMC)
registerDoMC(cores=3)

test<-cv.glmnet(x, y[,1:100], 
                family = "multinomial", 
                ## 0 < alpha < 1 
                alpha = 0.5, 
                type.measure="mse",
                nfold=kfold,
                grouped=FALSE,
                parallel=TRUE) 

plot(test)




