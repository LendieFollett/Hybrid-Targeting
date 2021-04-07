rm(list = ls())
library(mvtnorm)
library(MASS)
library(dplyr)
library(ggplot2)
library(truncnorm)
source("Bayes Consensus Ranking/functions.R")
R = 8  ## number of rankers
A = 2   ## number of aggregate/community-level variables captured
K = 20 ## number of communities
M = 2   ## number of micro-level variables captured
N0 = 100## number of unranked/training items
N1 = 60 ## number of ranked/test items
P = 13  ## number of covariates
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
omega_comm_true <- rep(1, A)
omega_micro_true <- rep(2, M)
omega_rank_true <- rep(.5, R)
beta_true = c(0,rep(1, P)) #first column is intercept

#Fill "responses"
gamma_comm_true <- rnorm(K, 0, 1) 
for (a in 1:A){ #fill community measures
  Y_comm[,a] <-  rnorm(K, X_comm %*% beta_true, sqrt(1/omega_comm_true[a])) +gamma_comm_true
}


gamma_micro_true <-  rnorm(N0, 0, 1)
for (m in 1:M){ #fill micro-data
  Y_micro[,m] <-  rnorm(N0, X_micro0 %*% beta_true, sqrt(1/omega_micro_true[m])) + gamma_micro_true
}

gamma_rank_true <-  rnorm(N1, 0, 1)
for (r in 1:R){ #fill latent Z scores
  Z[,r] <-  rnorm(N1, X_micro1 %*% beta_true, sqrt(1/omega_rank_true[r])) +gamma_rank_true
}


Tau <- apply(Z, 2, rank) #R rankings (what we actually observe)

pair.comp.ten = array(NA, dim = c(N1, N1, R)) ## get pairwise comparison matrices from the ranking lists
for(r in 1:R){
  pair.comp.ten[,,r] = FullRankToPairComp( Tau[,r] )
}


iter.max = 5000   ## Gibbs sampler total iterations
iter.burn =1000   ## Gibbs sampler burn-in iterations
print.opt = 100  ## print a message every print.opt steps

rm(A)
rm(N1)
rm(N0)
rm(R)
rm(M)
rm(P)

#Run MCMC for Bayesian Consensus Targeting
temp <- BCTarget(pair.comp.ten=pair.comp.ten, X_comm = X_comm, X_micro0 = X_micro0, X_micro1 = X_micro1,
                           Y_comm = Y_comm, Y_micro = Y_micro,
                               sigma_beta = 2.5,
                               weight.prior.value = c(0.5, 1, 2), 
                               weight.prior.prob = rep(1/length(weight.prior.value), length(weight.prior.value)),
                               N1 = dim(pair.comp.ten)[1], 
                               R = dim(pair.comp.ten)[3], 
                               iter.max = iter.max, para.expan = TRUE, print.opt = 100,
                               initial.list = NULL)

mu_mean <- apply(temp$mu, 2, mean)

#posterior summaries of ranks
tau_post <-apply(temp$mu, 1, rank)
tau_post_summary <- data.frame(
  mean =  rank(apply(temp$mu, 2, mean)),#Rank of posterior means of xbeta + gamma
  min = (apply(tau_post, 1, min)), #minimum rank seen in MCMC draws
  max = (apply(tau_post, 1, max)),#maximum rank seen in MCMC draws
  quantile = apply(tau_post, 1, quantile, .75)
)
tau_post_summary$naive_agg <- apply(Tau, 1, mean)
tau_post_summary$PMT <- rank(X_micro1%*%solve(t(X_micro0)%*%X_micro0)%*%t(X_micro0)%*%apply(Y_micro, 1, mean))



ggplot(data = tau_post_summary) +
  geom_pointrange(aes(x = naive_agg, y = mean,ymin = min, ymax = max)) +
  geom_abline(aes(slope = 1, intercept = 0)) +
  labs(x = "Mean Aggregation of R Ranks", y = "Posterior Summaries of Tau(alpha + X*Beta)")

ggplot(data = tau_post_summary) +
  geom_pointrange(aes(x = (PMT), y = mean,ymin = min, ymax = max)) +
  geom_abline(aes(slope = 1, intercept = 0)) +
  labs(x = "PMT-based ranks", y = "Posterior Summaries of Tau(alpha + X*Beta)")

ggplot(data = tau_post_summary[order(tau_post_summary$mean),]) +
  geom_point(aes(x = 1:nrow(tau_post_summary), y = mean)) +
  geom_point(aes(x = 1:nrow(tau_post_summary), y = naive_agg), colour = "tomato") +
  geom_point(aes(x = 1:nrow(tau_post_summary), y = PMT), colour = "steelblue") +
  geom_abline(aes(slope = 1, intercept = 0)) +
  labs(x = "ID", y = "Different rankings")



data.frame(postmean =  (apply(tau_post, 1, median)), rank(apply(Tau, 1, mean)))

#posteriors of quality weights - compare to truths
apply(temp$omega_comm, 2, mean)
apply(temp$omega_micro, 2, mean)
apply(temp$omega_rank, 2, mean) 


qplot(gamma_rank_true,apply(temp$gamma_rank, 2, mean)) +geom_abline(aes(intercept = 0, slope = 1))
qplot(gamma_micro_true,apply(temp$gamma_micro, 2, mean)) +geom_abline(aes(intercept = 0, slope = 1))
qplot(gamma_comm_true,apply(temp$gamma_comm, 2, mean)) +geom_abline(aes(intercept = 0, slope = 1))


plot(temp$sigma2_comm%>%sqrt)
plot(temp$sigma2_micro%>%sqrt)
plot(temp$sigma2_rank%>%sqrt)


