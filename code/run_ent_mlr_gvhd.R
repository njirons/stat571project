### run CV on GVHD data and record MSPE

rm(list=ls())

setwd("~/Documents/Classes/Winter 2021/stat571project")

library(caret)
library(glmnet)
library(MASS)

## necessary functions

# gradient ascent MLR + entropy regularization
ent_mlr <- function(x,y,n,ns,m,J,d,epsilon,eps_ns=TRUE,prop = FALSE,
                    cutoff = 10^(-5), iters = 5*10^3,
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
  return(list(beta=beta,obj_vals=obj_vals,iters = t))
}

# k-fold CV
mlr_cv <- function(x,y,n,ns,m,J,d,epsilons, eps_ns=TRUE, prop = FALSE,
                   cutoff = 10^(-5), iters = 5*10^3,
                   s=1,step="search",alpha = 0.5,gamma = 0.8, 
                   folds = 5, reps = 1, loss="cMSPE"){
  num_eps <- length(epsilons)
  # losses <- c("cMSPE","pMSPE","TV","KL","ChiSq")
  # cv_error <- matrix(0,ncol=length(losses),nrow=num_eps)
  cv_error <- rep(0,num_eps)
  
  for(i in 1:reps){
    # split data into folds
    flds <- createFolds(1:m, k = folds, list = TRUE, returnTrain = FALSE)
    # print(paste0("rep: ",round(100*i/reps),"%"))
    for(k in 1:folds){
      # print(paste0("rep: ",round(100*i/reps),
      #              "%, fold: ",round(100*k/folds),"%"))
      
      # split data into training and test sets
      x_train <- x[-(flds[[k]]),]
      x_test <- x[flds[[k]],]
      y_train <- y[-(flds[[k]]),]
      y_test <- y[flds[[k]],]
      for(j in 1:num_eps){
        print(paste0("rep: ",round(100*i/reps),
                     "%, fold: ",round(100*k/folds),
                     "%, epsilon: ",round(100*j/num_eps),"%"))
        
        # fit to training set
        fit <- ent_mlr(x=x_train,y=y_train,n=sum(y_train),ns=rowSums(y_train),
                       m=dim(x_train)[1],J=J,d=d,
                       epsilon=epsilons[j], eps_ns=eps_ns, prop = prop,
                       cutoff = cutoff, iters = iters,
                       s = s, step=step,alpha = alpha,gamma = gamma)
        beta <- fit$beta
        
        # print(fit$iters)
        plot(fit$obj_vals)
        
        # compute prediction error on test set
        num <- exp(x_test %*% beta)
        den <- 1+rowSums(num)
        phat <- cbind(num/den,1/den)
        ns_test <- rowSums(y_test)
        
        if(loss == "cMSPE"){
          #mean squared prediction error on counts
          # pred <- sum((y_test-ns_test*phat)^2)/(sum(y_test^2)*reps*folds)
          # pred <- sum((y_test-ns_test*phat)^2)/(n*m*J*reps*folds)
          pred <- mean((y_test-ns_test*phat)^2)
          # cv_error[j,1] <- cv_error[j,1] + pred 
        }
        if(loss == "pMSPE"){
          #mean squared prediction error on probabilities
          pred <- mean(((y_test/ns_test)-phat)^2)
          # cv_error[j,2] <- cv_error[j,2] + pred 
        }
        if(loss=="TV"){
          # TV distance / MAE of probabilities
          pred <- mean(abs((y_test/ns_test)-phat))
          # cv_error[j,3] <- cv_error[j,3] + pred 
        }
        if(loss=="KL"){
          # KL divergence/relative entropy
          pred <- -mean((y_test/ns_test)*log(phat))
          # cv_error[j,4] <- cv_error[j,4] + pred 
        }
        if(loss == "ChiSq"){
          #chi squared prediction error on counts
          pred <- mean((y_test-ns_test*phat)^2/(ns_test*phat))
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
              epsilon = epsilon_opt, iters = fit_opt$iters,
              cv_error=cv_error, epsilons=epsilons))
}

##### run our method for each study
cutoff <- 10^(-5)  # stopping criterion for gradient ascent
L <- 30 # size of epsilon grid to search over
epsilons <- exp(seq(-4,4,length.out=L))
iters <- 5*10^3
num_folds <- 3  # for k-fold CV
num_sets <- 50

set.seed(2021)

load("GvHD Data/xtest.RData")
load("GvHD Data/ytest.RData")
load("GvHD Data/xtrain.RData")
load("GvHD Data/ytrain.RData")

xtrains <- xtrain
ytrains <- ytrain
xtests <- xtest
ytests <- ytest

# fitted betahats and train/test splits
betahats <- list()

### mean squared prediction error on counts
mspe <- rep(0,num_sets)

### mean absolute prediction error on counts
mape <- rep(0,num_sets)

for(n in 1:num_sets){
  xtrain <- xtrains[[n]]
  ytrain <- ytrains[[n]]
  ns_train <- rowSums(ytrain)
  n_train <- sum(ns_train)
  m <- nrow(ytrain)
  J <- ncol(ytrain)
  d <- ncol(xtrain)
    
  xtest <- xtests[[n]]
  ytest <- ytests[[n]]
  ns_test <- rowSums(ytest)
  n_test <- sum(ns_test)
  
  cv <- mlr_cv(x=xtrain,y=ytrain,n=n_train,ns=ns_train,m=m,J=J,d=d,
               epsilons=epsilons, eps_ns = TRUE,
               cutoff = cutoff, iters = iters,
               s=1,step="search",alpha = 0.5,gamma = 0.8,
               folds = num_folds, reps = 1, loss="cMSPE")
  
  print(paste0("iters: ",cv$iters))
  print(paste0("log(epsilon_opt): ",log(cv$epsilon)))
  plot(log(cv$epsilons),cv$cv_error)
  
  betahat <- cv$beta
  betahats[[n]] <- betahat
  
  numhat <- exp(xtest %*% betahat)
  denhat <- 1+rowSums(numhat)
  phat <- cbind(numhat/denhat,1/denhat)
  yhat <- ns_test * phat
  mspe[n] <- mean((ytest-yhat)^2)
  mape[n] <- mean(abs(ytest-yhat))
  
  print(paste0("GVHD Sim: ", round(100*n/num_sets), "% done."))
}

save(mspe, file="GvHD data/mspe.RData")
save(mape, file="GvHD data/mape.RData")
save(betahats, file="GvHD data/betahats.RData")
