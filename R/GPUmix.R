library(mvtnorm)
library(MCMCpack)

MCMC = function(x,K,N,burnin){
  nr = nrow(x) ; nl = ncol(x); gamma=10; nu=1; Phi = diag(nl);mu0=rep(0,nl) ;alpha=1
  kappa0 = 1/gamma; nu0 = nu + 2; Lambda0 = solve(nu*Phi)
  #initialization
  Sigma = list();mu=list();
  for (j in 1:K) { 
    Sigma[[j]] = riwish(nu0,Lambda0)
    mu[[j]] = mvrnorm(1,mu0,gamma*Sigma[[j]])
  }
  
  pi = rep(1/K,K)
  
  z = c();n=rep(0,K);out = c()
  pdf_func = function(x,pi,mu,Sigma) {
    for (j in 1:K){
      out[j] = pi[j] * dmvnorm(x,mu[[j]],Sigma[[j]]) 
    }
    out = out/sum(out)   
    return(out)
  }
  
  #iterate
  for (step in 1:N) {
    
    for (i in 1:nr) {
      pdf = pdf_func(x[i,],pi,mu,Sigma)
      #resample indices
      #z[i] = sample(1:K,size=1, prob=pdf)
      z[i] = which.max(pdf)
    }
    
    for (k in 1:K) {
      n[k] = length(which(z==k))
    }
    
    #resample means and covariance matrices
    #prior
    
    
    #posterior parameters
    for (j in 1:K) {
      X = x[which(z==j),]
      if (n[j]==0) {
        Sigma[[j]] = riwish(nu0,Lambda0)
        mu[[j]] = mvrnorm(1,mu0,gamma*Sigma[[j]]) }
      else {
        if (n[j]==1) {X = t(X)}
        X_mean = colMeans(X)
        mun = kappa0/(kappa0+n[j])*mu0 + n[j]/(kappa0+n[j])* X_mean
        kappan = kappa0 + n[j]
        nun = nu0 + n[j]
        # Lambdan = Lambda0 + t(X-X_mean)%*%(X-X_mean) + kappa0*n[j]/(kappa0+n[j])*((X_mean-mu0)%*% t(X_mean-mu0))
        Lambdan = Lambda0 + t(X)%*%(diag(n[j])-matrix(1,n[j],n[j])/n[j])%*%(X) + kappa0*n[j]/(kappa0+n[j])*((X_mean-mu0)%*% t(X_mean-mu0))
        
        #simulate
        Sigma[[j]] = riwish(nun,Lambdan)
        mu[[j]] = mvrnorm(1,mun,1/kappan*Sigma[[j]])
      }
    }
    
    #resample mixture weights
    v = c() ;     
    for (j in 1:(K-1)){
      a = alpha + sum(n[(j+1):K])
      v[j] = rbeta(1,1+n[j],a)
      if (j ==1)  pi[j] = v[1]
      else pi[j] = v[j] * prod((1-v)[1:(j-1)])
      pi[K] = 1-sum(pi[1:(K-1)])
    }
    
    #resample DP precision
    #alpha?
    
    #  # take average of the parameters
    if (step == burnin){
      Sigmares = unlist(Sigma)
      mures = unlist(mu)
      pires = pi
    }
    if (step > burnin) {
      pires = (pires + pi)/2
      mures = (mures+unlist(mu))/2
      Sigmares = (Sigmares + unlist(Sigma))/2
    }
  }
  
  for (i in 1:K) {
    mu[[i]] = mures[((i-1)*nl+1):(nl*i)]
    Sigma[[i]] = matrix(Sigmares[((i-1)*nl^2+1):(nl^2*i)],nl,nl)   
  }
  return(list(pires = pi,mu = mu, Sigma = Sigma, z=z))
}
