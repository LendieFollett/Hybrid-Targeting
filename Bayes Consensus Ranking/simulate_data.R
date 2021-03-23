
library(mvtnorm)
library(dplyr)
R = 8  ## number of rankers
A = 2   ## number of aggregate/community-level variables captured
K = 10 ## number of communities
M = 2   ## number of micro-level variables captured
N0 = 100## number of unranked/training items
N1 = 50 ## number of ranked/test items
P = 3   ## number of covariates
rho=0.5 ## correlation for covariates

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
omega_comm_true <- rep(.5, A)
omega_micro_true <- rep(2, M)
omega_rank_true <- rep(1, R)
beta_true = c(0,rep(1, P)) #first column is intercept

#Fill "responses"
for (a in 1:A){ #fill community measures
  Y_comm[,a] <-  rnorm(K, X_comm %*% beta_true, omega_comm_true[a])
}

for (m in 1:M){ #fill micro-data
  Y_micro[,m] <-  rnorm(N0, X_micro0 %*% beta_true, omega_micro_true[m])
}

for (r in 1:R){ #fill latent Z scores
  Z[,r] <-  rnorm(N1, X_micro1 %*% beta_true, omega_rank_true[r])  
}

Tau <- apply(Z, 2, rank) #R rankings (what we actually observe)

