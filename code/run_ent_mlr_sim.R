#### Run simulation studies for entropic MLR

library(caret)
library(glmnet)
library(MASS)

## necessary functions

# gradient ascent MLR + entropy regularization
ent_mlr <- function(x,y,n,ns,m,J,d,epsilon,eps_ns=TRUE,prop = FALSE,
                    cutoff = 10^(-5), iters = 10^5,
                    s=1,step="search",alpha = 0.5,gamma = 0.8){
  t <- 0
  beta <- matrix(0,nrow=d,ncol=J-1)
  theta <- matrix(0,nrow=d,ncol=J-1)
  # num <- matrix(0,nrow=m,ncol=J-1)
  # p <- matrix(0,nrow=m,ncol=J)
  # H <- -rowSums(p*log(p))
  # grad <- matrix(1,nrow=d,ncol=J-1)
  
  if(prop){
    eps <- epsilon
    # const <- 1
    # const <- n^(-1/3)
    # const <- 1/mean(n/ns)
    # const <- mean(ns/n)
    const <- 1/100
    ns <- 1
  }else{
    const <- 1/n
    if(eps_ns){
      eps <- epsilon * ns
    }
    else{
      eps <- epsilon * n
    }
  }
  # print(const)
  
  increase <- Inf
  obj_tmp <- -Inf
  obj_vals <- c()
  
  # while(sqrt(sum(grad^2)) > cutoff & t < iters){
  while((increase > cutoff | increase <0) & (t < iters)){
    t <- t+1
    # num <- exp(x %*% beta)
    num <- exp(x %*% theta)
    den <- 1+rowSums(num)
    # p[,1:(J-1)] <- num/den
    # p[,J] <- 1/den
    p <- cbind(num/den,1/den)
    
    #compute gradient
    H <- -rowSums(p*log(p))
    yhat <- p*(ns+eps*(H+log(p))) 
    grad <- const*(t(x) %*% (y[,1:(J-1)]-yhat[,1:(J-1)]))
    obj_theta <- const*(sum(y*log(p))+sum(eps*H))
    # print(obj_theta)
    if(step=="search"){
      s <- 1/gamma #step size
      # backtracking line search
      p_tmp <- p
      # H_tmp <- H
      obj_step <- obj_theta
      # while(sum(y*log(p_tmp))/n+sum(eps*H_tmp)/n <
      #       obj+alpha*s*sum(grad^2)){
      # print(obj_step)
      # print(obj_step < obj_theta+alpha*s*sum(grad^2))
      # print(grad)
      while(obj_step < obj_theta+alpha*s*sum(grad^2)){
        s <- gamma*s
        # num_tmp <- exp(x %*% (beta+s*grad))
        num_tmp <- exp(x %*% (theta+s*grad))
        den_tmp <- 1+rowSums(num_tmp)
        # p_tmp[,1:(J-1)] <- num_tmp/den_tmp
        # p_tmp[,J] <- 1/den_tmp
        # print(p_tmp)
        p_tmp <- cbind(num_tmp/den_tmp,1/den_tmp)
        H_tmp <- -rowSums(p_tmp*log(p_tmp))
        obj_step <- const*(sum(y*log(p_tmp))+sum(eps*H_tmp))
        # print(obj_step)
      }
    }
    
    # gradient step
    # beta <- beta + s*grad
    
    # fast gradient step
    beta_tmp <- theta + s*grad
    theta <- beta_tmp + (t/(t+3))*(beta_tmp-beta)
    # theta <- beta_tmp
    beta <- beta_tmp
    
    #compute objective value
    num <- exp(x %*% beta)
    den <- 1+rowSums(num)
    # p[,1:(J-1)] <- num/den
    # p[,J] <- 1/den
    p <- cbind(num/den,1/den)
    H <- -rowSums(p*log(p))
    obj <- const*(sum(y*log(p))+sum(eps*H))
    obj_vals <- c(obj_vals,obj)
    increase = obj - obj_tmp
    obj_tmp <- obj
  }
  return(list(beta=beta,obj_vals=obj_vals))
}

# k-fold CV
mlr_cv <- function(x,y,n,ns,m,J,d,epsilons, eps_ns=TRUE, prop = FALSE,
                   cutoff = 10^(-5), iters = 10^5,
                   s=1,step="search",alpha = 0.5,gamma = 0.8, 
                   folds = 2, reps = 100, loss="cMSPE"){
  num_eps <- length(epsilons)
  # losses <- c("cMSPE","pMSPE","TV","KL","ChiSq")
  # cv_error <- matrix(0,ncol=length(losses),nrow=num_eps)
  cv_error <- rep(0,num_eps)
  
  for(i in 1:reps){
    # split data into folds
    flds <- createFolds(1:m, k = folds, list = TRUE, returnTrain = FALSE)
    print(paste0("rep: ",round(100*i/reps),"%"))
    for(k in 1:folds){
      # print(paste0("rep: ",round(100*i/reps),
      #              "%, fold: ",round(100*k/folds),"%"))
      
      # split data into training and test sets
      x_train <- x[-(flds[[k]]),]
      x_test <- x[flds[[k]],]
      y_train <- y[-(flds[[k]]),]
      y_test <- y[flds[[k]],]
      for(j in 1:num_eps){
        # print(paste0("rep: ",round(100*i/reps),
        #              "%, fold: ",round(100*k/folds),
        #              "%, epsilon: ",round(100*j/num_eps),"%"))
        
        # fit to training set
        fit <- ent_mlr(x=x_train,y=y_train,n=sum(y_train),ns=rowSums(y_train),
                       m=dim(x_train)[1],J=J,d=d,
                       epsilon=epsilons[j], eps_ns=eps_ns, prop = prop,
                       cutoff = cutoff, iters = iters,
                       s = s, step=step,alpha = alpha,gamma = gamma)
        beta <- fit$beta
        
        # compute prediction error on test set
        num <- exp(x_test %*% beta)
        den <- 1+rowSums(num)
        phat <- cbind(num/den,1/den)
        ns_test <- rowSums(y_test)
        
        if(loss == "cMSPE"){
          #mean squared prediction error on counts
          # pred <- sum((y_test-ns_test*phat)^2)/(sum(y_test^2)*reps*folds)
          # pred <- sum((y_test-ns_test*phat)^2)/(n*m*J*reps*folds)
          pred <- mean((y_test-ns_test*phat)^2)/n
          # cv_error[j,1] <- cv_error[j,1] + pred 
        }
        if(loss == "pMSPE"){
          #mean squared prediction error on probabilities
          pred <- sum(((y_test/ns_test)-phat)^2)/(reps*folds)
          # cv_error[j,2] <- cv_error[j,2] + pred 
        }
        if(loss=="TV"){
          # TV distance / MAE of probabilities
          pred <- sum(abs((y_test/ns_test)-phat))/(reps*folds)
          # cv_error[j,3] <- cv_error[j,3] + pred 
        }
        if(loss=="KL"){
          # KL divergence/relative entropy
          pred <- -sum((y_test/ns_test)*log(phat))/(reps*folds*100)
          # cv_error[j,4] <- cv_error[j,4] + pred 
        }
        if(loss == "ChiSq"){
          #chi squared prediction error on counts
          pred <- mean((y_test-ns_test*phat)^2/(ns_test*phat))/(m*J*reps*folds)
          # cv_error[j,5] <- cv_error[j,5] + pred 
        }
        cv_error[j] <- cv_error[j] + pred
      }
    }
  }
  
  # find optimal epsilon 
  # loss_num <- which(losses==loss)
  # epsilon_opt <- epsilons[which(cv_error[,loss_num]==min(cv_error[,loss_num]))]
  epsilon_opt <- epsilons[which(cv_error==min(cv_error))]
  
  # fit model with optimal epsilon on full data set
  fit_opt <- ent_mlr(x=x,y=y,n=n,ns=ns,m=m,J=J,d=d, 
                     eps_ns = eps_ns, prop = prop,
                     epsilon=epsilon_opt,cutoff = cutoff, iters = iters,
                     s = s, step=step,alpha = alpha,gamma = gamma)
  
  return(list(beta=fit_opt$beta, obj_vals = fit_opt$obj_vals,
              epsilon = epsilon_opt,
              cv_error=cv_error, epsilons=epsilons))
}

## generate fake data

set.seed(2021)
num_sets <- 50 # number of datasets to generate for each parameter setting
d <- 11; num_cov <- d-1;
Js <- c(25,50,100) # number of response categories
ms <- 2*Js # sample size/number of observations

#for sim 2
m2 <- 100 
J2 <- 50 
rhos <- c(0.1,0.5,0.9)   # correlation parameters

# for sim 3
sep <- 0.1 #fraction of linearly separable response categories

####### simulation 1: correct data generating mechanism

## We use d = 11 covariates (including the intercept)
## 5 of which are continuous (normal) and 5 binary

## We allow the number of observations to vary among m = 20, 50, 100
## We also vary the number of response categories J = 20, 50, 100

# generate data sets of each type
x1 <- list(m1J1 = list(), m1J2 = list(), m1J3 = list(),
           m2J1 = list(), m2J2 = list(), m2J3 = list(),
           m3J1 = list(), m3J2 = list(), m3J3 = list())
y1 <- list(m1J1 = list(), m1J2 = list(), m1J3 = list(),
           m2J1 = list(), m2J2 = list(), m2J3 = list(),
           m3J1 = list(), m3J2 = list(), m3J3 = list())
beta1 <- list(m1J1 = list(), m1J2 = list(), m1J3 = list(),
              m2J1 = list(), m2J2 = list(), m2J3 = list(),
              m3J1 = list(), m3J2 = list(), m3J3 = list())

for(i in 1:length(ms)){
  for(j in 1:length(Js)){
    index <- length(Js)*(i-1)+j
    m <- ms[i]
    J <- Js[j]
    for(n in 1:num_sets){
      # generate covariates independently
      bin_cov <- matrix(rbinom(n=m*num_cov/2,size=1,prob=1/2),nrow=m)
      cts_cov <- matrix(rnorm(n=m*num_cov/2,mean=0,sd=1),nrow=m) 
      x <- cbind(rep(1,m),bin_cov,cts_cov)
      x1[[index]][[n]] <- x
      
      beta <- matrix(runif(n=d*(J-1),min=-1,max=1),nrow=d)
      beta1[[index]][[n]] <- beta
      
      num <- exp(x %*% beta)
      den <- 1+rowSums(num)
      p <- cbind(num/den,1/den)
      y <- matrix(0,nrow=m,ncol=J)
      ns <- 10000*rpois(n=m,lambda=20)
      for(k in 1:m){
        y[k,] <- rmultinom(1,size=ns[k],prob=p[k,])
      }
      y1[[index]][[n]] <- y
    }
  }
}


####### simulation 2: (exhangeable) correlated covariates
# Similar setup as in simulation 1, but with correlated normal covariates

set.seed(2021)

m <- m2
J <- J2

# generate data sets of each type
x2 <- list(rho1 = list(),rho2 = list(),rho3 = list())
y2 <- list(rho1 = list(),rho2 = list(),rho3 = list())
beta2 <- list(rho1 = list(),rho2 = list(),rho3 = list())

for(r in 1:length(rhos)){
  rho <- rhos[r]
  Sigma <- ((1-rho)*diag(num_cov/2)+rho*matrix(1,nrow=num_cov/2,ncol=num_cov/2)) 
  # exchangeable covariance matrix with correlation rho
  for(n in 1:num_sets){
    # generate covariates 
    bin_cov <- matrix(rbinom(n=m*num_cov/2,size=1,prob=1/2),nrow=m)
    cts_cov <- mvrnorm(n=m, mu=rep(0,num_cov/2),Sigma = Sigma)
    x <- cbind(rep(1,m),bin_cov,cts_cov)
    x2[[r]][[n]] <- x
    
    beta <- matrix(runif(n=d*(J-1),min=-1,max=1),nrow=d)
    beta2[[r]][[n]] <- beta
    
    num <- exp(x %*% beta)
    den <- 1+rowSums(num)
    p <- cbind(num/den,1/den)
    y <- matrix(0,nrow=m,ncol=J)
    ns <- 10000*rpois(n=m,lambda=20)
    for(k in 1:m){
      y[k,] <- rmultinom(1,size=ns[k],prob=p[k,])
    }
    y2[[r]][[n]] <- y
  }
}


####### simulation 3: linearly separable data
# Similar setup as in simulation 1, but now we let a fixed fraction of the 
# response categories be linearly separable w.r.t. the binary covariates

set.seed(2021)

## We use d = 11 covariates (including the intercept)
## 5 of which are continuous (normal) and 5 binary

## We allow the number of observations to vary among m = 20, 50, 100
## We also vary the number of response categories J = 20, 50, 100

# generate data sets of each type
x3 <- list(m1J1 = list(), m1J2 = list(), m1J3 = list(),
           m2J1 = list(), m2J2 = list(), m2J3 = list(),
           m3J1 = list(), m3J2 = list(), m3J3 = list())
y3 <- list(m1J1 = list(), m1J2 = list(), m1J3 = list(),
           m2J1 = list(), m2J2 = list(), m2J3 = list(),
           m3J1 = list(), m3J2 = list(), m3J3 = list())
beta3 <- list(m1J1 = list(), m1J2 = list(), m1J3 = list(),
              m2J1 = list(), m2J2 = list(), m2J3 = list(),
              m3J1 = list(), m3J2 = list(), m3J3 = list())

for(i in 1:length(ms)){
  for(j in 1:length(Js)){
    index <- length(Js)*(i-1)+j
    m <- ms[i]
    J <- Js[j]
    for(n in 1:num_sets){
      # generate covariates independently
      bin_cov <- matrix(rbinom(n=m*num_cov/2,size=1,prob=1/2),nrow=m)
      cts_cov <- matrix(rnorm(n=m*num_cov/2,mean=0,sd=1),nrow=m) 
      x <- cbind(rep(1,m),bin_cov,cts_cov)
      x3[[index]][[n]] <- x
      
      beta <- matrix(runif(n=d*(J-1),min=-1,max=1),nrow=d)
      beta3[[index]][[n]] <- beta
      
      num <- exp(x %*% beta)
      num[which(bin_cov[,1]==1),1:(sep*J)] <- 0
      den <- 1+rowSums(num)
      p <- cbind(num/den,1/den)
      y <- matrix(0,nrow=m,ncol=J)
      ns <- 10000*rpois(n=m,lambda=20)
      for(k in 1:m){
        y[k,] <- rmultinom(1,size=ns[k],prob=p[k,])
      }
      y3[[index]][[n]] <- y
    }
  }
}


### save generated data sets

save(x1,file="x1_new.RData")
save(y1,file="y1_new.RData")
save(beta1,file="beta1_new.RData")

save(x2,file="x2_new.RData")
save(y2,file="y2_new.RData")
save(beta2,file="beta2_new.RData")

save(x3,file="x3_new.RData")
save(y3,file="y3_new.RData")
save(beta3,file="beta3_new.RData")



##### run our method for each study
cutoff <- 10^(-5)  # stopping criterion for gradient ascent
L <- 50  # size of epsilon grid to search over
epsilons <- exp(seq(-8,4,length.out=L))
num_folds <- 5  # for k-fold CV

### simulation 1
set.seed(2021)

#### mse of betas
mse1 <- list(m1J1 = rep(0,num_sets), m1J2 = rep(0,num_sets), m1J3 = rep(0,num_sets),
             m2J1 = rep(0,num_sets), m2J2 = rep(0,num_sets), m2J3 = rep(0,num_sets),
             m3J1 = rep(0,num_sets), m3J2 = rep(0,num_sets), m3J3 = rep(0,num_sets))

### mean squared prediction error on counts
mspe1 <- list(m1J1 = rep(0,num_sets), m1J2 = rep(0,num_sets), m1J3 = rep(0,num_sets),
              m2J1 = rep(0,num_sets), m2J2 = rep(0,num_sets), m2J3 = rep(0,num_sets),
              m3J1 = rep(0,num_sets), m3J2 = rep(0,num_sets), m3J3 = rep(0,num_sets))

#### mean absolute error of betas
mae1 <- list(m1J1 = rep(0,num_sets), m1J2 = rep(0,num_sets), m1J3 = rep(0,num_sets),
             m2J1 = rep(0,num_sets), m2J2 = rep(0,num_sets), m2J3 = rep(0,num_sets),
             m3J1 = rep(0,num_sets), m3J2 = rep(0,num_sets), m3J3 = rep(0,num_sets))

### mean absolute prediction error on counts
mape1 <- list(m1J1 = rep(0,num_sets), m1J2 = rep(0,num_sets), m1J3 = rep(0,num_sets),
              m2J1 = rep(0,num_sets), m2J2 = rep(0,num_sets), m2J3 = rep(0,num_sets),
              m3J1 = rep(0,num_sets), m3J2 = rep(0,num_sets), m3J3 = rep(0,num_sets))

for(i in 1:length(ms)){
  for(j in 1:length(Js)){
    index <- length(Js)*(i-1)+j
    m <- ms[i]
    J <- Js[j]
    for(n in 1:num_sets){
      x <- x1[[index]][[n]]
      y <- y1[[index]][[n]]
      beta <- beta1[[index]][[n]]
      
      trainid <- createDataPartition(1:m,times=1,p=0.8)
      
      xtrain <- x[trainid$Resample1,]
      ytrain <- y[trainid$Resample1,]
      ns_train <- rowSums(ytrain)
      n_train <- sum(ns_train)
      
      xtest <- x[-(trainid$Resample1),]
      ytest <- y[-(trainid$Resample1),]
      ns_test <- rowSums(ytest)
      n_test <- sum(ns_test)
      
      cv <- mlr_cv(x=xtrain,y=ytrain,n=n_train,ns=ns_train,m=m*0.8,J=J,d=d,
                   epsilons=epsilons, eps_ns = TRUE,
                   cutoff = cutoff, iters = 5*10^3,
                   s=1,step="search",alpha = 0.5,gamma = 0.8,
                   folds = num_folds, reps = 1, loss="cMSPE")
      
      betahat <- cv$beta
      mse1[[index]][n] <- mean((betahat-beta)^2)
      mae1[[index]][n] <- mean(abs(betahat-beta))
      
      numhat <- exp(xtest %*% beta)
      denhat <- 1+rowSums(numhat)
      phat <- cbind(numhat/denhat,1/denhat)
      yhat <- ns_test * phat
      mspe1[[index]][n] <- mean((ytest-yhat)^2)/(100*m*J)
      mape1[[index]][n] <- mean(abs(ytest-yhat))/(100*m*J)
    }
    print(paste0("Sim 1: ", round(100*index/(length(ms)*length(Js))), "% done."))
  }
}
save(mse1, file="mse1.RData")
save(mae1, file="mae1.RData")
save(mspe1, file="mspe1.RData")
save(mape1, file="mape1.RData")

### simulation 2
set.seed(2021)

#### mse of betas
mse2 <- list(rho1 = rep(0,num_sets),rho2 = rep(0,num_sets),rho3 = rep(0,num_sets))

### mean squared prediction error on counts
mspe2 <- list(rho1 = rep(0,num_sets),rho2 = rep(0,num_sets),rho3 = rep(0,num_sets))

#### mean absolute error of betas
mae2 <- list(rho1 = rep(0,num_sets),rho2 = rep(0,num_sets),rho3 = rep(0,num_sets))

### mean absolute prediction error on counts
mape2 <- list(rho1 = rep(0,num_sets),rho2 = rep(0,num_sets),rho3 = rep(0,num_sets))

for(r in 1:length(rhos)){
  for(n in 1:num_sets){
    x <- x2[[r]][[n]]
    y <- y2[[r]][[n]]
    beta <- beta2[[r]][[n]]
    
    trainid <- createDataPartition(1:m,times=1,p=0.8)
      
    xtrain <- x[trainid$Resample1,]
    ytrain <- y[trainid$Resample1,]
    ns_train <- rowSums(ytrain)
    n_train <- sum(ns_train)
      
    xtest <- x[-(trainid$Resample1),]
    ytest <- y[-(trainid$Resample1),]
    ns_test <- rowSums(ytest)
    n_test <- sum(ns_test)
      
    cv <- mlr_cv(x=xtrain,y=ytrain,n=n_train,ns=ns_train,m=m*0.8,J=J,d=d,
                 epsilons=epsilons, eps_ns = TRUE,
                 cutoff = cutoff, iters = 5*10^3,
                 s=1,step="search",alpha = 0.5,gamma = 0.8,
                 folds = num_folds, reps = 1, loss="cMSPE")
      
    betahat <- cv$beta
    mse2[[r]][n] <- mean((betahat-beta)^2)
    mae2[[r]][n] <- mean(abs(betahat-beta))
      
    numhat <- exp(xtest %*% beta)
    denhat <- 1+rowSums(numhat)
    phat <- cbind(numhat/denhat,1/denhat)
    yhat <- ns_test * phat
    mspe2[[r]][n] <- mean((ytest-yhat)^2)/(100*m*J)
    mape2[[r]][n] <- mean(abs(ytest-yhat))/(100*m*J)
  }
  print(paste0("Sim 2: ", round(100*r/(length(rhos))), "% done."))
}

save(mse2, file="mse2.RData")
save(mae2, file="mae2.RData")
save(mspe2, file="mspe2.RData")
save(mape2, file="mape2.RData")

### simulation 3
set.seed(2021)

#### mse of betas
mse3 <- list(m1J1 = rep(0,num_sets), m1J2 = rep(0,num_sets), m1J3 = rep(0,num_sets),
             m2J1 = rep(0,num_sets), m2J2 = rep(0,num_sets), m2J3 = rep(0,num_sets),
             m3J1 = rep(0,num_sets), m3J2 = rep(0,num_sets), m3J3 = rep(0,num_sets))

### mean squared prediction error on counts
mspe3 <- list(m1J1 = rep(0,num_sets), m1J2 = rep(0,num_sets), m1J3 = rep(0,num_sets),
              m2J1 = rep(0,num_sets), m2J2 = rep(0,num_sets), m2J3 = rep(0,num_sets),
              m3J1 = rep(0,num_sets), m3J2 = rep(0,num_sets), m3J3 = rep(0,num_sets))

#### mean absolute error of betas
mae3 <- list(m1J1 = rep(0,num_sets), m1J2 = rep(0,num_sets), m1J3 = rep(0,num_sets),
             m2J1 = rep(0,num_sets), m2J2 = rep(0,num_sets), m2J3 = rep(0,num_sets),
             m3J1 = rep(0,num_sets), m3J2 = rep(0,num_sets), m3J3 = rep(0,num_sets))

### mean absolute prediction error on counts
mape3 <- list(m1J1 = rep(0,num_sets), m1J2 = rep(0,num_sets), m1J3 = rep(0,num_sets),
              m2J1 = rep(0,num_sets), m2J2 = rep(0,num_sets), m2J3 = rep(0,num_sets),
              m3J1 = rep(0,num_sets), m3J2 = rep(0,num_sets), m3J3 = rep(0,num_sets))

for(i in 1:length(ms)){
  for(j in 1:length(Js)){
    index <- length(Js)*(i-1)+j
    m <- ms[i]
    J <- Js[j]
    for(n in 1:num_sets){
      x <- x3[[index]][[n]]
      y <- y3[[index]][[n]]
      beta <- beta3[[index]][[n]]
      
      trainid <- createDataPartition(1:m,times=1,p=0.8)
      
      xtrain <- x[trainid$Resample1,]
      ytrain <- y[trainid$Resample1,]
      ns_train <- rowSums(ytrain)
      n_train <- sum(ns_train)
      
      xtest <- x[-(trainid$Resample1),]
      ytest <- y[-(trainid$Resample1),]
      ns_test <- rowSums(ytest)
      n_test <- sum(ns_test)
      
      cv <- mlr_cv(x=xtrain,y=ytrain,n=n_train,ns=ns_train,m=m*0.8,J=J,d=d,
                   epsilons=epsilons, eps_ns = TRUE,
                   cutoff = cutoff, iters = 5*10^3,
                   s=1,step="search",alpha = 0.5,gamma = 0.8,
                   folds = num_folds, reps = 1, loss="cMSPE")
      
      betahat <- cv$beta
      mse3[[index]][n] <- mean((betahat-beta)^2)
      mae3[[index]][n] <- mean(abs(betahat-beta))
      
      numhat <- exp(xtest %*% beta)
      denhat <- 1+rowSums(numhat)
      phat <- cbind(numhat/denhat,1/denhat)
      yhat <- ns_test * phat
      mspe3[[index]][n] <- mean((ytest-yhat)^2)/(100*m*J)
      mape3[[index]][n] <- mean(abs(ytest-yhat))/(100*m*J)
    }
    print(paste0("Sim 3: ", round(100*index/(length(ms)*length(Js))), "% done."))
  }
}
save(mse3, file="mse3.RData")
save(mae3, file="mae3.RData")
save(mspe3, file="mspe3.RData")
save(mape3, file="mape3.RData")

