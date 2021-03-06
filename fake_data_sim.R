m <- 20
J <- 50
p <- 5
ns <- rep(100000,m)
n <- sum(ns)

# separable data
x0 <- matrix(0,nrow=m/2,ncol=p)
for(i in 1:(m/2)){
  x0[i,] <- c(1,0,rnorm(p-2))
}

x1 <- matrix(0,nrow=m/2,ncol=p)
for(i in 1:(m/2)){
  x1[i,] <- c(1,1,rnorm(p-2))
}

beta0 <- matrix(rnorm(p*(J-1),mean=4,sd=1),nrow=p,ncol=J-1)

p0 <- matrix(0,nrow=m/2,ncol=J)
num0 <- exp(x0 %*% beta0)
den0 <- 1+rowSums(num0)
p0[,1:(J-1)] <- num0/den0
p0[,J] <- 1/den0

p1 <- matrix(0,nrow=m/2,ncol=J)
num1 <- exp(x1 %*% beta0)
num1[,1] <- 0
den1 <- 1+rowSums(num1)
p1[,1:(J-1)] <- num1/den1
p1[,J] <- 1/den1

y0 <- matrix(0,nrow=m/2,ncol=J)
for(i in 1:(m/2)){
  y0[i,] <- rmultinom(1,size=ns[i],prob=p0[i,])
}

y1 <- matrix(0,nrow=m/2,ncol=J)
for(i in 1:(m/2)){
  y1[i,] <- rmultinom(1,size=ns[i],prob=p1[i,])
}

y <- rbind(y0,y1)
x <- rbind(x0,x1)

# correct data generating mechanism







# run simulation

# gradient ascent MLR
beta <- matrix(0,nrow=dim(x)[2],ncol=J-1)
num <- matrix(0,nrow=m,ncol=J-1)
p <- matrix(0,nrow=m,ncol=J)
t <- 0
grad <- matrix(1,nrow=dim(x)[2],ncol=J-1)

alpha <- 0.5
gamma <- 0.8
while(sqrt(sum(grad^2)) > 10^(-5)){
  t <- t+1
  s <- 1/gamma
  num <- exp(x %*% beta)
  den <- 1+rowSums(num)
  p[,1:(J-1)] <- num/den
  p[,J] <- 1/den
  
  #compute gradient
  yhat <- ns * p
  grad <- t(x) %*% (y[,1:(J-1)]-yhat[,1:(J-1)]) / n
  
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
  if(t %% 10000 == 0){
    print(max(beta))
    print(sqrt(sum(grad^2)))
  }
}

beta0[2,1]
beta[2,1]

# gradient ascent MLR + entropy regularization
beta_e <- matrix(0,nrow=dim(x)[2],ncol=J-1)
num <- matrix(0,nrow=m,ncol=J-1)
p <- matrix(0,nrow=m,ncol=J)
t <- 0
grad <- matrix(1,nrow=dim(x)[2],ncol=J-1)

alpha <- 0.5
gamma <- 0.8
s <- 1/gamma

epsilon <- exp(log(n)/2)
eps <- epsilon

# epsilon <- 0.2
# eps <- epsilon*ns

while(sqrt(sum(grad^2)) > 10^(-5)){
  t <- t+1
  s <- 1/gamma
  num <- exp(x %*% beta_e)
  den <- 1+rowSums(num)
  p[,1:(J-1)] <- num/den
  p[,J] <- 1/den
  
  #compute gradient
  H <- -rowSums(p*log(p))
  #yhat <- (ns*p)*(1+epsilon*(H+log(p)))
  # yhat <- p*(ns+eps*(H+log(p)))
  yhat <- (p*ns) +(eps*(p*log(p))) +(eps*H)*p
  grad <- t(x) %*% (y[,1:(J-1)]-yhat[,1:(J-1)]) / n
  
  # backtracking line search
  p_tmp <- p
  H_tmp <- H
  while(sum(y*log(p_tmp))/n+sum(eps*H_tmp)/n <
        sum(y*log(p))/n+sum(eps*H)/n+alpha*s*sum(grad^2)){
    s <- gamma*s
    num_tmp <- exp(x %*% (beta_e+s*grad))
    den_tmp <- 1+rowSums(num_tmp)
    p_tmp[,1:(J-1)] <- num_tmp/den_tmp
    p_tmp[,J] <- 1/den_tmp
    H_tmp <- -rowSums(p_tmp*log(p_tmp))
  }
  
  # gradient step
  beta_e <- beta_e + s*grad
  if(t %% 10000 == 0){
    print(max(beta_e))
    print(sqrt(sum(grad^2)))
  }
}

beta0[2,1]
beta_e[2,1]




# gradient ascent MLR + ridge regression
beta_r <- matrix(0,nrow=dim(x)[2],ncol=J-1)
num <- matrix(0,nrow=m,ncol=J-1)
p <- matrix(0,nrow=m,ncol=J)
t <- 0
grad <- matrix(1,nrow=dim(x)[2],ncol=J-1)

alpha <- 0.5
gamma <- 0.8

lambda <- 1/n*100
while(sqrt(sum(grad^2)) > 10^(-5)){
  t <- t+1
  s <- 1/gamma
  num <- exp(x %*% beta_r)
  den <- 1+rowSums(num)
  p[,1:(J-1)] <- num/den
  p[,J] <- 1/den
  
  #compute gradient
  yhat <- ns * p
  grad <- t(x) %*% (y[,1:(J-1)]-yhat[,1:(J-1)]) / n - lambda*beta_r
  
  #line search
  p_tmp <- p
  beta_tmp <- beta_r
  while(sum(y*log(p_tmp))/n-lambda/2*sum(beta_tmp^2) <
        sum(y*log(p))/n-lambda/2*sum(beta_r^2)+alpha*s*sum(grad^2)){
    s <- gamma*s
    beta_tmp <- beta_r+s*grad
    num_tmp <- exp(x %*% beta_tmp)
    den_tmp <- 1+rowSums(num_tmp)
    p_tmp[,1:(J-1)] <- num_tmp/den_tmp
    p_tmp[,J] <- 1/den_tmp
  }
  beta_r <- beta_r + s*grad
  if(t %% 10000 == 0){
    print(max(beta_r))
    print(sqrt(sum(grad^2)))
  }
}

beta0[2,1]
beta_r[2,1]


sqrt(sum((beta-beta_r)^2))/(sqrt(sum((beta_r)^2)))
abs(beta-beta_r)

sqrt(sum((beta0-beta_e)^2))/(sqrt(sum((beta0)^2)))
abs(beta-beta_r)

