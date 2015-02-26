library(mvtnorm)
library(MCMCpack)
devtools::load_all()
##data generation
mu = list()
#cluster means
mu[[1]] = c(5,1,2);mu[[2]]=c(-5,-1,-2);mu[[3]] = c(0,2,0)
n = 1000
pi = c(0.6,0.3,0.1)
x = matrix(,n,3);index=c()
for (i in 1:n){
  index[i] = sample(1:3,1,prob=pi)
  x[i,] = mvrnorm(1,mu[[index[i]]],diag(3))  #identity covariance matrices
}

K=5
res = MCMC(x,K,100,50)
plot(x)
points(x,col=res$z)


y=matrix(,n,3)
for (i in 1:n) {
  index = sample(1:K,1,prob=res$pires)
  y[i,] = mvrnorm(1,res$mu[[index]],res$Sigma[[index]])
}

par(mfrow=c(2,3))
hist(x[,1],breaks=100)
hist(x[,2],breaks=100)
hist(x[,3],breaks=100)
hist(y[,1],breaks=52,xlim=c(-10,10))
hist(y[,2],breaks=100,xlim=c(-5,5))
hist(y[,3],breaks=100,xlim=c(-5,5))
