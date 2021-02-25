library(truncnorm)
library(mvtnorm)
M = 10  ## number of rankers
N = 50  ## number of ranked items
L = 3   ## number of covariates
rho=0.5   ## correlation for covariates
CovMat=diag(L) ## covariance matrix for the covariates
for(i in 1:(L-1)){
  for(j in (i+1):L){
    CovMat[i,j]=rho^(abs(i-j))
    CovMat[j,i]=rho^(abs(i-j))
  }
}
X.mat = rmvnorm(N, mean = rep(0, L), sigma = CovMat) ## covariate matrix 
beta.true = c(3,2,1)
mu.true = rowSums( X.mat^2 ) + as.vector(  X.mat %*% beta.true )  ## true evaluation score
rank.true = rank(mu.true)  ## true ranking list
sigma.true = 5  ## noise level
Z.real = t( rmvnorm(M, mean = mu.true, sigma = sigma.true^2 * diag(N) ) ) ## scores for all rankers
fullrank.real = apply(Z.real, 2, rank)  ## OBSERVED ranking lists