library(truncnorm)
library(mvtnorm)
M = 10  ## number of rankers
N0 = 50  ## number of ranked items (number of entities in 'test' sample)
N1 = 70 ## number of items involved in PMT sample (in 'training' sample)
N = N0 + N1 ## total number of sampled individuals
L = 3   ## number of covariates
rho=0.5   ## correlation for covariates
CovMat=diag(L) ## covariance matrix for the covariates
for(i in 1:(L-1)){
  for(j in (i+1):L){
    CovMat[i,j]=rho^(abs(i-j))
    CovMat[j,i]=rho^(abs(i-j))
  }
}
X.mat = rmvnorm(N, mean = rep(0, L), sigma = CovMat) ## covariate matrix for ALL sampled individuals
X.mat1 = X.mat[1:N1,] ## covariate matrix for ALL sampled individuals
Xmat0 = X.mat[-c(1:N1),]
beta.true = c(3,2,1)
mu.true0 = rowSums( X.mat0^2 ) + as.vector(  X.mat0 %*% beta.true )  ## true evaluation score
mu.true1 = rowSums( X.mat1^2 ) + as.vector(  X.mat1 %*% beta.true )  ## true evaluation score
rank.true0 = rank(mu.true0)  ## true ranking list
rank.true1 = rank(mu.true1)  ## true ranking list
sigma.true = 5  ## noise level
Z.real0 = t( rmvnorm(M, mean = mu.true0, sigma = sigma.true^2 * diag(N0) ) ) ## scores for test rankers
Z.real1 = t( rmvnorm(M, mean = mu.true1, sigma = sigma.true^2 * diag(N1) ) ) ## scores for training rankers
Z_obs1 = Z.real0[1,] #OBSERVED scores. e.g., consumption for 'training' sample
fullrank.real0 = apply(Z.real0, 2, rank)  ## OBSERVED ranking lists for 'testing' individuals

