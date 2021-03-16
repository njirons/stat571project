########################## just "futzing around"
m <- 200; J <- 10;
d <- 10
ns <- rpois(n=m,lambda=(1:m)*10000)
# ns <- rpois(n=m,lambda=100000)
n <- sum(ns)

# separable data
x0 <- matrix(0,nrow=m/2,ncol=d)
for(i in 1:(m/2)){
  x0[i,] <- c(1,0,rnorm(d-2))
}

x1 <- matrix(0,nrow=m/2,ncol=d)
for(i in 1:(m/2)){
  x1[i,] <- c(1,1,rnorm(d-2))
}

beta0 <- matrix(rnorm(d*(J-1),mean=0,sd=2),nrow=d,ncol=J-1)

p0 <- matrix(0,nrow=m/2,ncol=J)
num0 <- exp(x0 %*% beta0)
den0 <- 1+rowSums(num0)
p0[,1:(J-1)] <- num0/den0
p0[,J] <- 1/den0

p1 <- matrix(0,nrow=m/2,ncol=J)
num1 <- exp(x1 %*% beta0)
# num1[,1] <- 0
den1 <- 1+rowSums(num1)
p1[,1:(J-1)] <- num1/den1
p1[,J] <- 1/den1

y0 <- matrix(0,nrow=m/2,ncol=J)
for(i in 1:(m/2)){
  y0[i,] <- rmultinom(1,size=ns[i],prob=p0[i,])
}

y1 <- matrix(0,nrow=m/2,ncol=J)
for(i in 1:(m/2)){
  y1[i,] <- rmultinom(1,size=ns[i+m/2],prob=p1[i,])
}

y <- rbind(y0,y1)
x <- rbind(x0,x1)

rowSums(y)/ns
###########################################

### Generate and save data sets for simulation studies

library(MASS)

####### simulation 1: correct data generating mechanism

set.seed(2021)
num_sets <- 50 # number of datasets to generate for each parameter setting

## We use d = 11 covariates (including the intercept)
## 5 of which are continuous (normal) and 5 binary

## We allow the number of observations to vary among m = 20, 50, 100
## We also vary the number of response categories J = 20, 50, 100
d <- 11; num_cov <- d-1;
ms <- 5*c(20,50,100)
Js <- c(20,50,100)

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
    index <- length(ms)*(i-1)+j
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
save(x1,file="x1.RData")
save(y1,file="y1.RData")
save(beta1,file="beta1.RData")


####### simulation 2: (exhangeable) correlated covariates
# Similar setup as in simulation 1, but with correlated normal covariates

set.seed(2021)
num_sets <- 50 # number of datasets to generate for each parameter setting

d <- 11; num_cov <- d-1;
m <- 5*50
J <- 50
rhos <- c(0.1,0.5,0.9)   # correlation parameters

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
save(x2,file="x2.RData")
save(y2,file="y2.RData")
save(beta2,file="beta2.RData")


####### simulation 3: linearly separable data
# Similar setup as in simulation 1, but now we let a fixed fraction of the 
# response categories be linearly separable w.r.t. the binary covariates

set.seed(2021)
num_sets <- 50 # number of datasets to generate for each parameter setting

## We use d = 11 covariates (including the intercept)
## 5 of which are continuous (normal) and 5 binary

## We allow the number of observations to vary among m = 20, 50, 100
## We also vary the number of response categories J = 20, 50, 100
d <- 11; num_cov <- d-1;
ms <- 5*c(20,50,100)
Js <- c(20,50,100)
sep <- 0.1 #fraction of linearly separable response categories

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
    index <- length(ms)*(i-1)+j
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
save(x3,file="x3.RData")
save(y3,file="y3.RData")
save(beta3,file="beta3.RData")




