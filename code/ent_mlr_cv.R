library(caret)
library(glmnet)

# gradient ascent MLR + entropy regularization
ent_mlr <- function(x,y,n,ns,m,J,d,epsilon,eps_ns=TRUE,prop = FALSE,
                    cutoff = 10^(-5), iters = 10^5,
                    s=1,step="search",alpha = 0.5,gamma = 0.8){
  t <- 0
  beta <- matrix(0,nrow=d,ncol=J-1)
  theta <- matrix(0,nrow=d,ncol=J-1)
  # num <- matrix(0,nrow=m,ncol=J-1)
  p <- matrix(0,nrow=m,ncol=J)
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
    p[,1:(J-1)] <- num/den
    p[,J] <- 1/den
    # p <- cbind(num/den,1/den)
    
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
        p_tmp[,1:(J-1)] <- num_tmp/den_tmp
        p_tmp[,J] <- 1/den_tmp
        # print(p_tmp)
        # p_tmp <- cbind(num_tmp/den_tmp,1/den_tmp)
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
    p[,1:(J-1)] <- num/den
    p[,J] <- 1/den
    # p <- cbind(num/den,1/den)
    H <- -rowSums(p*log(p))
    obj <- const*(sum(y*log(p))+sum(eps*H))
    obj_vals <- c(obj_vals,obj)
    increase = obj - obj_tmp
    obj_tmp <- obj
  }
  return(list(beta=beta,obj_vals=obj_vals))
}

## examples
beta <- ent_mlr(x,y/ns,n,ns,m,J,d,epsilon=exp(-4),eps_ns=TRUE, prop=TRUE,
                cutoff = 10^(-5), iters = 5*10^3,
                s = 0.3, step="search",alpha = 0.5,gamma = 0.8)
plot(beta$obj_vals)

beta <- ent_mlr(x,y,n,ns,m,J,d,epsilon=exp(0),eps_ns=TRUE, prop=FALSE,
                cutoff = 10^(-5), iters = 5*10^3,
                s = 0.3, step="search",alpha = 0.5,gamma = 0.8)
plot(beta$obj_vals)

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
    for(k in 1:folds){
      print(paste0("rep: ",round(100*i/reps),
                   "%, fold: ",round(100*k/folds),"%"))
      
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

## examples
cv_c <- mlr_cv(x=x,y=y,n=n,ns=ns,m=m,J=J,d=d,
               epsilons=exp(seq(-8,4,length.out=20)),
               eps_ns=TRUE, prop=FALSE,
               # epsilons = 0,
               cutoff = 10^(-4), iters = 5*10^3,
               s=1,step="search",alpha = 0.5,gamma = 0.8,
               folds = 5, reps = 1, loss="cMSPE")

cv_c <- mlr_cv(x=x,y=y,n=n,ns=ns,m=m,J=J,d=d,
               epsilons=exp(seq(-1,4,length.out=10)), 
               eps_ns=TRUE, prop=FALSE,
               # epsilons = 0,
               cutoff = 10^(-5), iters = 5*10^3,
               s=1,step="search",alpha = 0.5,gamma = 0.8, 
               folds = 2, reps = 1, loss="ChiSq")

cv_p <- mlr_cv(x=x,y=y/ns,n=n,ns=ns,m=m,J=J,d=d,
               epsilons=exp(seq(-10,5,length.out=50)),
               eps_ns=TRUE, prop=TRUE,
               # epsilons = 0,
               cutoff = 10^(-5), iters = 5*10^3,
               s=1,step="search",alpha = 0.5,gamma = 0.8,
               folds = 10, reps = 10, loss="pMSPE")

cv_e <- mlr_cv(x,y,n,ns,m,J,d,
               epsilons=exp(seq(-10,-6,length.out=10)),eps_ns=FALSE, prop=FALSE,
               cutoff = 10^(-5), iters = 5*10^3,
               s=1,step="search",alpha = 0.5,gamma = 0.8,
               folds = 2, reps = 10, loss="cMSPE")

cv <- cv_c
cv <- cv_p
cv <- cv_e

#plot CV errors
plot(log(cv$epsilons),cv$cv_error,
     type="l",col=1, ylab="CV error", xlab="log(epsilon)")
log(cv$epsilon)
plot(cv$obj_vals)
# plot(log(cv$epsilons),
#      (cv$cv_error[,1]-min(cv$cv_error[,1]))/
#        (max(cv$cv_error[,1])-min(cv$cv_error[,1])),
#      type="l",col=1,ylim=c(0,1), ylab="CV error", xlab="log(epsilon)")
# lines(log(cv$epsilons),
#      (cv$cv_error[,2]-min(cv$cv_error[,2]))/
#        (max(cv$cv_error[,2])-min(cv$cv_error[,2])),
#      type="l",col=2)
# lines(log(cv$epsilons),
#      (cv$cv_error[,3]-min(cv$cv_error[,3]))/
#        (max(cv$cv_error[,3])-min(cv$cv_error[,3])),
#      type="l",col=3)
# lines(log(cv$epsilons),
#      (cv$cv_error[,4]-min(cv$cv_error[,4]))/
#        (max(cv$cv_error[,4])-min(cv$cv_error[,4])),
#      type="l",col=4)

# beta0 <- cv$beta

# mean((cv$beta-beta0)^2/abs(beta0))
# mean(abs((cv$beta-beta0)/beta0))
# # mean((cv$beta-beta0)^2/abs(cv$beta))
# mean((cv$beta-beta0)^2/beta0^2)
# mean((cv$beta-beta0)^2/cv$beta^2)
mean((cv$beta-beta0)^2)
mean(abs(cv$beta-beta0))

max(abs(cv$beta))
max(abs(beta0))

s <- 1
plot(beta0[s,],ylim=c(min(c(beta0[s,],cv$beta[s,])),
                      max(c(beta0[s,],cv$beta[s,]))))
points(cv$beta[s,],col=2)
s <- s+1

s <- 1
plot(beta0[,s],ylim=c(min(c(beta0[,s],cv$beta[,s])),
                      max(c(beta0[,s],cv$beta[,s]))))
points(cv$beta[,s],col=2)
s <- s+1


### glmnet
fit <- cv.glmnet(x,y,nfolds=5,family="multinomial", #nlambda=100,
                 maxit=10^6,type.measure = "dev",alpha=0.5,
                 intercept=FALSE)

length(fit$lambda)
dim(fit$glmnet.fit$beta[[100]])

fit$lambda.min
which(fit$lambda==fit$lambda.min)
beta_glm <- fit$glmnet.fit$beta[[fit$index[1]]]
sum((beta_glm[,2:dim(beta_glm)[2]]-beta0)^2)

fit <- cv.glmnet(x,y,nfolds=5,family="multinomial", #nlambda=10,
                 maxit=10^6,type.measure = "mse", alpha=0.5)

length(fit$lambda)
dim(fit$glmnet.fit$beta[[1]])


















# coordinate descent with fixed Hessian newton-raphson iterations
ent_mlr_coord <- function(x,y,n,ns,m,J,d,epsilon,eps_ns=TRUE,
                    cutoff = 10^(-5), iters = 10^5,inv_hess){
  t <- 0
  beta <- matrix(0,nrow=d,ncol=J-1)
  num <- matrix(1,nrow=m,ncol=J-1)
  den <- rep(J,m)
  p <- matrix(1/J,nrow=m,ncol=J)
  H <- -rep(log(1/J),m)
  
  if(eps_ns){
    eps <- epsilon * ns
  }
  else{
    eps <- epsilon * n
  }
  
  increase <- 1
  obj_tmp <- -999999
  obj_vals <- c()
  
  while((increase > cutoff | increase <0) & (t < iters)){
    t <- t+1
    for(j in 1:(J-1)){
      #compute gradient
      # yhat <- p*(ns+eps*(H+log(p)))
      yhat <- p[,j]*(ns+eps*(H+log(p[,j])))
      # grad <- t(x) %*% (y[,j]-yhat[,j])
      grad <- t(x) %*% (y[,j]-yhat)
      
      # newton-raphson step
      beta[,j] <- beta[,j] - inv_hess %*% grad
      
      # update
      den <- den - num[,j]
      num[,j] <- exp(x %*% beta[,j])
      # den <- 1+rowSums(num)
      den <- den + num[,j]
      p[,1:(J-1)] <- num/den
      p[,J] <- 1/den
      H <- -rowSums(p*log(p))
    }
    
    #compute objective value
    obj <- sum(y*log(p))/n+sum(eps*H)/n
    obj_vals <- c(obj_vals,obj)
    increase = obj - obj_tmp
    obj_tmp <- obj
  }
  return(list(beta=beta,obj_vals=obj_vals))
}

epsilon <- exp(1)
grad <- ent_mlr(x,y,n,ns,m,J,d,epsilon=epsilon,eps_ns=TRUE,
                cutoff = 10^(-5), iters = 5*10^3,
                s = 0.3, step="search",alpha = 0.5,gamma = 0.8)
plot(grad$obj_vals)

A <- sqrt(ns*(1+epsilon))*x
# hess <- -(t(A) %*% A)*((1/J)*(1-(1/J)))
hess <- -(t(A) %*% A)*((1/2)*(1-(1/J)))
inv_hess <- solve(hess)
coord <- ent_mlr_coord(x,y,n,ns,m,J,d,epsilon=epsilon,eps_ns=TRUE,
              cutoff = 10^(-5), iters = 10^5,inv_hess=inv_hess)
plot(coord$obj_vals)

