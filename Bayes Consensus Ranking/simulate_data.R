
library(mvtnorm)

R = 10  ## number of rankers
A = 2   ## number of aggregate/community-level variables captured
K = 3 ## number of communities
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
X.mat = rmvnorm(N0 + N1, mean = rep(0, P), sigma = CovMat) ## covariate matrix for ALL sampled individuals
X.mat1 = X.mat[1:N1,] ## covariate matrix for ALL sampled individuals
X.mat0 = X.mat[-c(1:N1),]

X_comm <- aggregate(X.mat, list(community), mean)

Y_comm <- array(NA, dim = c(K, A))
Y_micro <- array(NA, dim = c(N0, A)) #only training has micro response (e.g., consumption)

beta.true = c(3,2,1)
mu.true0 =  as.vector(  X.mat0 %*% beta.true )  ## true evaluation score
mu.true1 =  as.vector(  X.mat1 %*% beta.true )  ## true evaluation score
rank.true0 = rank(mu.true0)  ## true ranking list
rank.true1 = rank(mu.true1)  ## true ranking list
sigma.true = 5  ## noise level
Z.real1 = t( rmvnorm(M, mean = mu.true1, sigma = sigma.true^2 * diag(N1) ) ) ## scores for test rankers (wouldn't observe in real life)
#Z_obs1 = rmvnorm(1, mean = mu.true0, sigma = sigma.true^2 * diag(N0) ) #OBSERVED scores. e.g., consumption for 'training' sample
fullrank.real0 = apply(Z.real1, 2, rank)  ## OBSERVED ranking lists for 'testing' individuals

