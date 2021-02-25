rm(list = ls())
library(truncnorm)
library(mvtnorm)

source("Bayes Consensus Ranking/simulate_data.R")
source("Bayes Consensus Ranking/functions.R")
#Set up
pair.comp.ten = array(NA, dim = c(N0, N0, M)) ## get pairwise comparison matrices from the ranking lists
for(j in 1:M){
  pair.comp.ten[,,j] = FullRankToPairComp( fullrank.real0[,j] )
}
X.mat.sd0 = t( (t( X.mat0 ) - colMeans(X.mat)) / apply(X.mat, 2, sd) )  ## standardized covariates
X.mat.sd1 = t( (t( X.mat1 ) - colMeans(X.mat)) / apply(X.mat, 2, sd) )  ## standardized covariates
iter.max = 1000   ## Gibbs sampler total iterations
iter.burn = 200   ## Gibbs sampler burn-in iterations
print.opt = 100  ## print a message every print.opt steps

BARCW.fit = BayesRankCovWeight(pair.comp.ten = pair.comp.ten, 
                               X.mat = X.mat.sd, 
                               Z_obs = Z_obs,
                               tau2.alpha = 1^2, nu.alpha = 3,
                               tau2.beta = 10^2, nu.beta = 3,
                               iter.max = iter.max, print.opt = print.opt)

BARCW.fit$agg.rank = apply(BARCW.fit$mu[, -c(1:iter.burn)], 1, mean)  ## aggregated ranking list
RankDist(BARCW.fit$agg.rank, rank.true)   ## Kendall tau distance between estimated and true ranking lists

rowMeans( BARCW.fit$weight.vec )  ## posterior means of weights for all rankers
