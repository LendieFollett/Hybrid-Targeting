library(mvtnorm)
library(MASS)
library(dplyr)
library(ggplot2)
library(truncnorm)



community <- rep(1:K, (N0+N1)/K) #each training + testing household assigned to a community

CovMat=diag(P) ## covariance matrix for the covariates
for(i in 1:(P-1)){
  for(j in (i+1):P){
    CovMat[i,j]=rho^(abs(i-j))
    CovMat[j,i]=rho^(abs(i-j))
  }
}
# X_MICRO - both training and testing has same micro covariate information
X_micro = rmvnorm(N0 + N1, mean = rep(0, P), sigma = CovMat) ## covariate matrix for ALL sampled individuals
X_micro <- cbind(1, X_micro)
X_micro1 = X_micro[1:N1,] ## covariate matrix for ALL sampled individuals
X_micro0 = X_micro[-c(1:N1),]

X_comm <- aggregate(X_micro, list(community), mean)[,-1] %>%
  as.matrix()

Y_comm <- array(NA, dim = c(K, A)) 
Y_micro <- array(NA, dim = c(N0, M)) #only training has micro response (e.g., consumption)
Z <- array(NA, dim = c(N1, R)) #only testing has latent ranks (e.g., consumption)

#parameter values
omega_comm_true <- rep(1, A)
omega_micro_true <- rep(2, M)
omega_rank_true <- rep(.5, R)
beta_true = c(0,rep(1, P)) #first column is intercept

#random effect variances
sigma2_comm <- sigma2_micro <- sigma2_rank <- 1
#Fill "responses"
gamma_comm_true <- rnorm(K, 0, sigma2_comm^.5) 
for (a in 1:A){ #fill community measures
  Y_comm[,a] <-  rnorm(K, X_comm %*% beta_true, sqrt(1/omega_comm_true[a])) +gamma_comm_true
}


gamma_micro_true <-  rnorm(N0, 0, sigma2_micro^.5)
for (m in 1:M){ #fill micro-data
  Y_micro[,m] <-  rnorm(N0, X_micro0 %*% beta_true, sqrt(1/omega_micro_true[m])) + gamma_micro_true
}

gamma_rank_true <-  rnorm(N1, 0, sigma2_rank^.5)
for (r in 1:R){ #fill latent Z scores
  Z[,r] <-  rnorm(N1, X_micro1 %*% beta_true, sqrt(1/omega_rank_true[r])) +gamma_rank_true
}



Tau <- apply(Z, 2, rank) #R rankings (what we actually observe)

#real life observe: Tau, Y_micro ('test' only), Y_comm 
#do NOT observe: Z


pair.comp.ten = array(NA, dim = c(N1, N1, R)) ## get pairwise comparison matrices from the ranking lists
for(r in 1:R){
  pair.comp.ten[,,r] = FullRankToPairComp( Tau[,r] )
}


iter.max = 5000   ## Gibbs sampler total iterations
iter.burn =1000   ## Gibbs sampler burn-in iterations
print.opt = 100  ## print a message every print.opt steps

#remove parameters we wouldn't have
rm(A)
rm(N1)
rm(N0)
rm(R)
rm(M)
rm(P)
