rm(list = ls())
library(truncnorm)
library(mvtnorm)
library(LaplacesDemon)
source("Bayes Consensus Ranking/functions.R")

#parameters for simulation
R = 5  ## number of rankers
A = 2   ## number of aggregate/community-level variables captured
K = 20 ## number of communities
M = 2   ## number of micro-level variables captured
N0 = 100## number of unranked/training items
N1 = 60 ## number of ranked/test items
P = 13  ## number of covariates
rho=0.5 ## correlation for covariates

iter.keep = 10000   ## Gibbs sampler kept iterations (post burn-in)
iter.burn =5000   ## Gibbs sampler burn-in iterations 
print.opt = 100  ## print a message every print.opt steps

#simulate data based on parameters
source("Bayes Consensus Ranking/simulate_data.R")

#Run MCMC for Bayesian Consensus Targeting
temp <- BCTarget(pair.comp.ten=pair.comp.ten, X_comm = X_comm, X_micro0 = X_micro0, X_micro1 = X_micro1,
                 Y_comm = Y_comm, Y_micro = Y_micro,
                 sigma_beta = 2.5,
                 weight.prior.value = c(0.5, 1, 2), 
                 weight.prior.prob = rep(1/length(weight.prior.value), length(weight.prior.value)),
                 N1 = dim(pair.comp.ten)[1], 
                 R = dim(pair.comp.ten)[3], 
                 iter.keep = iter.keep,
                 iter.burn = iter.burn,
                 para.expan = TRUE, print.opt = 100,
                 initial.list = NULL)

mu_mean <- apply(temp$mu, 2, mean) # mu = X_micro1 %*% beta + Xr_rank[1:N1,]%*%gamma_rank

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
  labs(x = "Mean Aggregation of R Ranks", y = "Posterior Summaries of Tau(alpha + X*Beta + gamma)")


ggplot(data = tau_post_summary) +
  geom_pointrange(aes(x = PMT, y = mean,ymin = min, ymax = max, colour = apply(temp$gamma_rank, 2, mean))) +
  scale_colour_gradient2("Gamma\nPosterior\nMean")+
  geom_abline(aes(slope = 1, intercept = 0)) +
  labs(x = "Tau(alpha + X*Beta)", y = " Tau(alpha + X*Beta + gamma)") +
  ggtitle("Effect of Random Effects on Ultimate Ranks")

ggplot(data = tau_post_summary) +
  geom_point(aes(x = X_micro1%*%solve(t(X_micro0)%*%X_micro0)%*%t(X_micro0)%*%apply(Y_micro, 1, mean),
                      y = X_micro1%*%solve(t(X_micro0)%*%X_micro0)%*%t(X_micro0)%*%apply(Y_micro, 1, mean) +apply(temp$gamma_rank, 2, mean),
                      colour = apply(temp$gamma_rank, 2, mean))) +
  scale_colour_gradient2("Gamma\nPosterior\nMean")+
  geom_abline(aes(slope = 1, intercept = 0)) +
  labs(x = "alpha + X*Beta", y = "alpha + X*Beta + gamma")+
  ggtitle("Effect of Random Effects on Posterior Means")

ggplot(data = tau_post_summary[order(tau_post_summary$mean),]) +
  geom_point(aes(x = 1:nrow(tau_post_summary), y = mean)) +
  geom_point(aes(x = 1:nrow(tau_post_summary), y = naive_agg), colour = "tomato") +
  geom_point(aes(x = 1:nrow(tau_post_summary), y = PMT), colour = "steelblue") +
  geom_abline(aes(slope = 1, intercept = 0)) +
  labs(x = "ID", y = "Different rankings")



data.frame(postmean =  (apply(tau_post, 1, median)), rank(apply(Tau, 1, mean)))

#posteriors of quality weights - compare to truths
apply(temp$omega_comm, 2, mean); omega_comm_true
apply(temp$omega_micro, 2, mean);omega_micro_true
apply(temp$omega_rank, 2, mean) ;omega_rank_true


###convergence diagnostics------
qplot(gamma_rank_true,apply(temp$gamma_rank, 2, mean), 
      colour = apply(temp$gamma_rank, 2, ESS)) +
  scale_colour_gradient2("ESS", midpoint = 100)+#color by effective sample size - want > 100
  geom_abline(aes(intercept = 0, slope = 1))

qplot(gamma_micro_true,apply(temp$gamma_micro, 2, mean),
      colour = apply(temp$gamma_micro, 2, ESS)) +
  scale_colour_gradient2("ESS", midpoint = 100)+#color by effective sample size - want > 100
  geom_abline(aes(intercept = 0, slope = 1))
qplot(gamma_comm_true,apply(temp$gamma_comm, 2, mean) ,
      colour = apply(temp$gamma_comm, 2, ESS)) +
  scale_colour_gradient2("ESS", midpoint = 100)+#color by effective sample size - want > 100
  geom_abline(aes(intercept = 0, slope = 1))

plot(temp$sigma2_comm%>%sqrt)
plot(temp$sigma2_micro%>%sqrt)
plot(temp$sigma2_rank%>%sqrt)

plot(temp$gamma_rank[,4])


#do all ESS checks

doESS <- function(x){
  if(!is.null(dim(x))){ #if it's a data frame
    return(apply(x, 2, ESS))
  }else{
    return(ESS(x))
  }
}

lapply(temp[-c(1,7)], doESS) %>% str()
