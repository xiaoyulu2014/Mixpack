library(mvtnorm)
library(MCMCpack)
##data generation
mu = list()
#cluster means
mu[[1]] = c(5,1);mu[[2]]=c(-5,-1)
n = 1000
pi = c(0.6,0.4)
x = matrix(,n,2);index=c()
for (i in 1:n){
  index[i] = sample(1:2,1,prob=pi)
  x[i,] = mvrnorm(1,mu[[index[i]]],diag(2))  #identity covariance matrices
}

K=5
res = MCMC(x,K,100,50)
plot(x,col=res$z)

# x1 = x[which(index==1),]
# x_mean = colMeans(x1)
# x1cov = t(x1-x_mean)%*%(x1-x_mean)
# riwish(length(which(index==1)),x1cov)
# cov(x1,x1)
# 
# t(x1)%*%(diag(length(which(index==1)))-matrix(1,length(which(index==1)),length(which(index==1)))/length(which(index==1)))%*%(x1)



y=matrix(,n,2)
for (i in 1:n) {
  index = sample(1:K,1,prob=res$pires)
  y[i,] = mvrnorm(1,res$mu[[index]],res$Sigma[[index]])
}

par(mfrow=c(2,2))
hist(x[,1],breaks=100)
hist(x[,2],breaks=100)
hist(y[,1],breaks=52)
hist(y[,2],breaks=100)
# 
# n = 1000
# y=matrix(,n,2)
# for (i in 1:n) {
#   index = sample(1:2,1,prob=pi)
#   y[i,] = mvrnorm(1,mu[[index]],Sigma[[index]])
# }