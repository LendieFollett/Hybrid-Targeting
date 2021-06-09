rm(list = ls())
library(truncnorm)
library(mvtnorm)
library(LaplacesDemon)
library(lme4)
library(Matrix) #for sparse matrices
library(MASS)
library(dplyr)
library(ggplot2)
library(Rcpp)

#devtools::install_github("adzemski/rtnorm")
sourceCpp("functions.cpp")
source("Bayes Consensus Ranking/functions.R")

#parameters for simulation

A = 2   ## number of aggregate/community-level variables captured
K = R = 200 ## number of communities equal to number of rankers
M = 2   ## number of micro-level variables captured
N0 = 1000## number of unranked/training items
N1 = 1000 ## number of ranked/test items
P = 6  ## number of covariates
rho=0 ## correlation for covariates

iter.keep = 100   ## Gibbs sampler kept iterations (post burn-in)
iter.burn =100   ## Gibbs sampler burn-in iterations 
print.opt = 1  ## print a message every print.opt steps

#simulate data based on parameters
source("Bayes Consensus Ranking/simulate_data.R")


#starting values for random effects
temp_data <- data.frame(y = as.vector(apply(apply(Tau, 2, function(x){(x - mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)}), 1, function(x){sum(x, na.rm=TRUE)})), kronecker(rep(1, ncol(Z)), X_micro1))
form <- formula(paste0("y~-1+", paste0("X", 1:ncol(X_micro1), collapse = "+")))
#gamma_start <- ranef(lmer(form, data = temp_data))[[1]]$`(Intercept)` 
beta_start <-coef(lm(form, data = temp_data))%>%as.vector()
initial_list <- list(#gamma_rank = gamma_start,
                     beta = beta_start)

#Run MCMC for Bayesian Consensus Targeting
temp <- BCTarget(pair.comp.ten=pair.comp.ten, X_comm = X_comm, X_micro0 = X_micro0, X_micro1 = X_micro1,
                 X_elite = 7,#7th position (including intercept)
                 Y_comm = Y_comm, Y_micro = Y_micro,
                 weight.prior.value = c(0.5, 1, 2), 
                 weight.prior.prob = rep(1/length(weight.prior.value), length(weight.prior.value)),
                 iter.keep = iter.keep,
                 iter.burn = iter.burn,
                 para.expan = TRUE, print.opt = 100,
                 initial.list = initial_list)

mu_mean <- apply(temp$mu, 2, mean) # mu = X_micro1 %*% beta + Xr_rank[1:N1,]%*%gamma_rank

#posterior summaries of ranks
tau_post <-apply(temp$mu, 1, rank)
tau_post_summary <- data.frame(
  mean =  rank(apply(temp$mu, 2, mean)),#Rank of posterior means of xbeta + gamma
  min = (apply(tau_post, 1, min)), #minimum rank seen in MCMC draws
  max = (apply(tau_post, 1, max)),#maximum rank seen in MCMC draws
  quantile = apply(tau_post, 1, quantile, .75)
)
tau_post_summary$naive_agg <- apply(Tau, 1, function(x){mean(x, na.rm=TRUE)})
tau_post_summary$PMT <- rank(X_micro1%*%solve(t(X_micro0)%*%X_micro0)%*%t(X_micro0)%*%apply(Y_micro, 1, mean))



ggplot(data = tau_post_summary) +
  geom_violin(aes(x = naive_agg, y = mean, group = naive_agg)) +
  geom_jitter(aes(x = naive_agg, y = mean)) +
  #geom_abline(aes(slope = 1, intercept = 0)) +
  labs(x = "Rank within community", 
       y = "Posterior Mean of Tau(X*Beta )")


#ggplot(data = tau_post_summary) +
#  geom_pointrange(aes(x = PMT, y = mean,ymin = min, ymax = max, colour = apply(temp$gamma_rank, 2, mean))) +
#  scale_colour_gradient2("Gamma\nPosterior\nMean")+
#  geom_abline(aes(slope = 1, intercept = 0)) +
#  labs(x = "Tau(alpha + X*Beta)", y = " Tau(alpha + X*Beta + gamma)") +
#  ggtitle("Effect of Random Effects on Ultimate Ranks")

#ggplot(data = tau_post_summary) +
#  geom_point(aes(x = X_micro1%*%solve(t(X_micro0)%*%X_micro0)%*%t(X_micro0)%*%apply(Y_micro, 1, mean),
#                      y = X_micro1%*%solve(t(X_micro0)%*%X_micro0)%*%t(X_micro0)%*%apply(Y_micro, 1, mean) +apply(temp$gamma_rank, 2, mean),
#                      colour = apply(temp$gamma_rank, 2, mean))) +
#  scale_colour_gradient2("Gamma\nPosterior\nMean")+
#  geom_abline(aes(slope = 1, intercept = 0)) +
#  labs(x = "alpha + X*Beta", y = "alpha + X*Beta + gamma")+
#  ggtitle("Effect of Random Effects on Posterior Means")

ggplot(data = tau_post_summary[order(tau_post_summary$mean),]) +
  geom_point(aes(x = 1:nrow(tau_post_summary), y = mean-PMT,colour = naive_agg)) +
  #geom_point(aes(x = 1:nrow(tau_post_summary), y = naive_agg), colour = "tomato") +
  #geom_point(aes(x = 1:nrow(tau_post_summary), y = PMT,colour = naive_agg)) +
  #geom_abline(aes(slope = 1, intercept = 0)) +
  labs(x = "ID", y = "Different rankings")



data.frame(postmean =  (apply(tau_post, 1, median)), rank(apply(Tau, 1, mean)))

#posteriors of quality weights - compare to truths
apply(temp$omega_comm, 2, mean); omega_comm_true
apply(temp$omega_micro, 2, mean);omega_micro_true
apply(temp$omega_rank, 2, mean) ;omega_rank_true

apply(temp$beta, 2, mean) ;beta_true

###convergence diagnostics------
#qplot(gamma_rank_true,apply(temp$gamma_rank, 2, mean), 
#      colour = apply(temp$gamma_rank, 2, ESS)) +
#  scale_colour_gradient2("ESS", midpoint = 100)+#color by effective sample size - want > 100
#  geom_abline(aes(intercept = 0, slope = 1))

qplot(gamma_micro_true,apply(temp$gamma_micro, 2, mean),
      colour = apply(temp$gamma_micro, 2, ESS)) +
  scale_colour_gradient2("ESS", midpoint = 100)+#color by effective sample size - want > 100
  geom_abline(aes(intercept = 0, slope = 1))
qplot(gamma_comm_true,apply(temp$gamma_comm, 2, mean) ,
      colour = apply(temp$gamma_comm, 2, ESS)) +
  scale_colour_gradient2("ESS", midpoint = 100)+#color by effective sample size - want > 100
  geom_abline(aes(intercept = 0, slope = 1))

plot(temp$sigma2_comm%>%sqrt); sigma2_comm %>%sqrt
plot(temp$sigma2_micro%>%sqrt); sigma2_micro%>%sqrt
#plot(temp$sigma2_rank%>%sqrt); sigma2_rank%>%sqrt 
#these are consistently overestimated - reconsider the inverse chi-squared prior

plot(temp$gamma_rank[,4])


#do all ESS checks

doESS <- function(x){
  if(!is.null(dim(x))){ #if it's a data frame
    return(apply(x, 2, ESS))
  }else{
    return(ESS(x))
  }
}

effectiv_ss <- lapply(temp[-c(1)], doESS) 

effectiv_ss%>% str()

effectiv_ss
