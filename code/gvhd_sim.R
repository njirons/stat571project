setwd("~/Documents/Classes/Winter 2021/stat 571/final")
library(readr)
library(tidyr)
library(dplyr)
library(nnet)
library(mlogit)
library(reshape2)
library(foreign)
library(stats4)

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

# m <- 50
m <- 138
J <- 100
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
d <- ncol(x)

x <- x[subjects,]
taxa <- taxon[subjects,1:J]

x <- as.matrix(x)
y <- as.matrix(taxa)

ns <- rowSums(y)
n <- sum(ns)

rowSums(y/ns)

# m <- nrow(taxa)
# J <- ncol(taxa)





















# gradient ascent
y <- as.matrix(taxa)
x <- as.matrix(x)
ns <- apply(y,1,sum)
#ns <- matrix(rep(ns,each=J),nrow=m,ncol=J,byrow=TRUE)
n <- sum(ns)

beta <- matrix(0,nrow=dim(x)[2],ncol=J-1)
num <- matrix(0,nrow=m,ncol=J-1)
p <- matrix(0,nrow=m,ncol=J)
t <- 0
# iters <- 1000
# step <- 0.1/n
grad <- matrix(1,nrow=dim(x)[2],ncol=J-1)

alpha <- 0.5
gamma <- 0.8
# for(t in 1:iters){
while(sqrt(sum(grad^2)) > 10^(-4)){
  t <- t+1
  s <- 1/gamma
  # s <- 0.05
  num <- exp(x %*% beta)
  den <- 1+rowSums(num)
  p[,1:(J-1)] <- num/den
  p[,J] <- 1/den
  
  #compute gradient
  yhat <- ns * p
  grad <- t(x) %*% (y[,1:(J-1)]-yhat[,1:(J-1)]) / n
  
  # yhat2 <- matrix(0,nrow=m,ncol=J)
  # for(i in 1:m){
  #   for(j in 1:J){
  #     #yhat2[i,j] <- ns[i,j]*p[i,j]
  #     yhat2[i,j] <- ns[i]*p[i,j]
  #   }
  # }
  
  # for(i in 1:m){
  #   for(j in 1:(J-1)){
  #     grad[,j] <- grad[,j] + x[i,]*(y[i,j]-ns[i]*p[i,j])/n
  #   }
  # }

  #line search
  p_tmp <- p
  while(sum(y*log(p_tmp))/n < sum(y*log(p))/n+alpha*s*sum(grad^2)){
    s <- gamma*s
    num_tmp <- exp(x %*% (beta+s*grad))
    den_tmp <- 1+rowSums(num_tmp)
    p_tmp[,1:(J-1)] <- num_tmp/den_tmp
    p_tmp[,J] <- 1/den_tmp
  }
  beta <- beta + s*grad
  if(t %% 100 == 0){
    print(max(beta))
    print(sqrt(sum(grad^2)))
  }
}


# gradient ascent + entropy regularization
y <- as.matrix(taxa)
x <- as.matrix(x)
ns <- apply(y,1,sum)
n <- sum(ns)

beta <- matrix(0,nrow=dim(x)[2],ncol=J-1)
num <- matrix(0,nrow=m,ncol=J-1)
p <- matrix(0,nrow=m,ncol=J)
t <- 0
# iters <- 1000
# step <- 0.1/n
grad <- matrix(1,nrow=dim(x)[2],ncol=J-1)

alpha <- 0.5
gamma <- 0.8
s <- 1/gamma

# epsilon <- 10000
epsilon <- 100
eps <- epsilon*ns
# eps <- epsilon

# for(t in 1:iters){
while(sqrt(sum(grad^2)) > 10^(-4)){
# while(sqrt(sum((s*grad)^2)) > 10^(-8)){
  t <- t+1
  #s <- 1/gamma
  s <- 0.01
  num <- exp(x %*% beta)
  den <- 1+rowSums(num)
  p[,1:(J-1)] <- num/den
  p[,J] <- 1/den
  
  #compute gradient
  H <- -rowSums(p*log(p))
  #yhat <- (ns*p)*(1+epsilon*(H+log(p)))
  # yhat <- p*(ns+eps*(H+log(p)))
  yhat <- (p*ns) +(eps*(p*log(p))) +(eps*H)*p
  grad <- t(x) %*% (y[,1:(J-1)]-yhat[,1:(J-1)]) / n
  
  # for(i in 1:m){
  #   for(j in 1:(J-1)){
  #     grad[,j] <- grad[,j] + x[i,]*(y[i,j]-ns[i]*p[i,j])/n #ns[i]
  #     for(k in 1:J){
  #       if(k == j){
  #         grad[,j] <- grad[,j] - eps*(1+log(p[i,j]))*
  #           (x[i,]*num[i,j]*(1+sum(num[i,-j])))/(den[i]^2*n) #ns[i])
  #       }
  #       if(k == J){
  #         grad[,j] <- grad[,j] + eps*(1-log(den[i]))*
  #           (x[i,]*num[i,j])/(den[i]^2*n) #ns[i])
  #       }
  #       if((k != j) & (k < J)){
  #         grad[,j] <- grad[,j] + eps*x[i,]*((1+log(p[i,k]))*
  #                                             (num[i,j]*num[i,k])/(den[i]^2*n)) #ns[i]))
  #       }
  #     }
  #   }
  # }
  
  # backtracking line search
  # p_tmp <- p
  # H_tmp <- H
  # # while(sum(y*log(p_tmp))/n-eps*sum(p_tmp*log(p_tmp))/n <
  # #       sum(y*log(p))/n-eps*sum(p*log(p))/n+alpha*s*sum(grad^2)){
  # # while(sum(y*log(p_tmp)+eps*H_tmp)/n <
  # #       sum(y*log(p)+eps*H)/n+alpha*s*sum(grad^2)){
  # while(sum(y*log(p_tmp)+eps*H_tmp)/n <= sum(y*log(p)+eps*H)/n){
  #   s <- gamma*s
  #   num_tmp <- exp(x %*% (beta+s*grad))
  #   den_tmp <- 1+rowSums(num_tmp)
  #   p_tmp[,1:(J-1)] <- num_tmp/den_tmp
  #   p_tmp[,J] <- 1/den_tmp
  #   H_tmp <- -rowSums(p_tmp*log(p_tmp))
  # }
  
  # gradient step
  beta <- beta + s*grad
  if(t %% 100 == 0){
    print(s)
    print(max(beta))
    print(sqrt(sum(grad^2)))
  }
}






# gradient ascent + ridge regression
y <- as.matrix(taxa)
x <- as.matrix(x)
ns <- apply(y,1,sum)
#ns <- matrix(rep(ns,each=J),nrow=m,ncol=J,byrow=TRUE)
n <- sum(ns)

beta_r <- matrix(0,nrow=dim(x)[2],ncol=J-1)
num <- matrix(0,nrow=m,ncol=J-1)
p <- matrix(0,nrow=m,ncol=J)
t <- 0
# iters <- 1000
# step <- 0.1/n
grad <- matrix(1,nrow=dim(x)[2],ncol=J-1)

alpha <- 0.5
gamma <- 0.8

lambda <- 0.01
# for(t in 1:iters){
while(sqrt(sum(grad^2)) > 10^(-4)){
  t <- t+1
  # s <- 1/gamma
  s <- 0.1
  num <- exp(x %*% beta_r)
  den <- 1+rowSums(num)
  p[,1:(J-1)] <- num/den
  p[,J] <- 1/den
  
  #compute gradient
  yhat <- ns * p
  grad <- t(x) %*% (y[,1:(J-1)]-yhat[,1:(J-1)]) / n - lambda*beta_r
  
  # yhat2 <- matrix(0,nrow=m,ncol=J)
  # for(i in 1:m){
  #   for(j in 1:J){
  #     #yhat2[i,j] <- ns[i,j]*p[i,j]
  #     yhat2[i,j] <- ns[i]*p[i,j]
  #   }
  # }
  
  # for(i in 1:m){
  #   for(j in 1:(J-1)){
  #     grad[,j] <- grad[,j] + x[i,]*(y[i,j]-ns[i]*p[i,j])/n
  #   }
  # }
  
  #line search
  # p_tmp <- p
  # beta_tmp <- beta
  # while(sum(y*log(p_tmp))/n-lambda/2*sum(beta_tmp^2) < 
  #       sum(y*log(p))/n-lambda/2*sum(beta^2)+alpha*s*sum(grad^2)){
  #   s <- gamma*s
  #   beta_tmp <- beta+s*grad
  #   num_tmp <- exp(x %*% beta_tmp)
  #   den_tmp <- 1+rowSums(num_tmp)
  #   p_tmp[,1:(J-1)] <- num_tmp/den_tmp
  #   p_tmp[,J] <- 1/den_tmp
  # }
  beta_r <- beta_r + s*grad
  if(t %% 100 == 0){
    print(max(beta_r))
    print(sqrt(sum(grad^2)))
  }
}


sqrt(sum((beta-beta_r)^2))/(sqrt(sum((beta_r)^2)))
abs(beta-beta_r)
