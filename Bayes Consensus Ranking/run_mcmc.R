library(truncnorm)
library(mvtnorm)

source("Bayes Consensus Ranking/simulate_data.R")

#Set up
pair.comp.ten = array(NA, dim = c(N, N, M)) ## get pairwise comparison matrices from the ranking lists
for(j in 1:M){
  pair.comp.ten[,,j] = FullRankToPairComp( fullrank.real[,j] )
}
X.mat.sd = t( (t( X.mat ) - colMeans(X.mat)) / apply(X.mat, 2, sd) )  ## standardized covariates
iter.max = 1000   ## Gibbs sampler total iterations
iter.burn = 200   ## Gibbs sampler burn-in iterations
print.opt = 100  ## print a message every print.opt steps

