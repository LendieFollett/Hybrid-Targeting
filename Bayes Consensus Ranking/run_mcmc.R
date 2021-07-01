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
#sourceCpp("functions.cpp")
source("Bayes Consensus Ranking/functions.R")

#parameters for simulation

A = 2   ## number of aggregate/community-level variables captured
K = R = 200 ## number of communities equal to number of rankers
M = 2   ## number of micro-level variables captured
N0 = 1000## number of unranked/training items
N1 = 1000 ## number of ranked/test items
P = 6  ## number of covariates
rho=0 ## correlation for covariates

iter.keep = 1000   ## Gibbs sampler kept iterations (post burn-in)
iter.burn =500   ## Gibbs sampler burn-in iterations 
print.opt = 1  ## print a message every print.opt steps

#simulate data based on parameters
source("Bayes Consensus Ranking/simulate_data.R")


#starting values for random effects
temp_data <- data.frame(y = as.vector(apply(apply(Y_micro, 2, function(x){(x - mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)}), 1, function(x){mean(x, na.rm=TRUE)})),
                        kronecker(rep(1, ncol(Y_micro)), X_micro0))
form <- formula(paste0("y~-1+", paste0("X", 1:ncol(X_micro1), collapse = "+")))
#gamma_start <- ranef(lmer(form, data = temp_data))[[1]]$`(Intercept)` 
beta_start <-coef(lm(form, data = temp_data))%>%as.vector()
initial_list <- list(#gamma_rank = gamma_start,
                     beta_rank = beta_start,
                     beta_comm = beta_start,
                     beta_micro = beta_start)

#Run MCMC for Bayesian Consensus Targeting
temp <- BCTarget(Tau=Tau, 
                 X_micro0 = X_micro0, 
                 X_micro1 = X_micro1,
                 X_elite = 7,#7th position (including intercept)
                 Y_micro = Y_micro,
                 iter.keep = iter.keep,
                 iter.burn = iter.burn,
                 para.expan = TRUE, print.opt = 100)

mu_mean <- apply(temp$mu, 2, mean) # mu = X_micro1 %*% mu_beta

#posterior summaries of ranks
tau_post <-apply(temp$mu, 1, rank)
tau_post_summary <- data.frame(
  mean_rank =  rank(apply(temp$mu, 2, mean)),#Rank of posterior means of xbeta + gamma
  min_rank = (apply(tau_post, 1, min)), #minimum rank seen in MCMC draws
  max_rank = (apply(tau_post, 1, max)),#maximum rank seen in MCMC draws
  mean_score = mu_mean
)
tau_post_summary$naive_rank <- apply(Tau, 1, function(x){mean(x, na.rm=TRUE)})
tau_post_summary$PMT_rank <- rank(X_micro1%*%solve(t(X_micro0)%*%X_micro0)%*%t(X_micro0)%*%apply(Y_micro, 1, mean))

tau_post_summary$PMT <- (X_micro1%*%solve(t(X_micro0)%*%X_micro0)%*%t(X_micro0)%*%apply(Y_micro, 1, mean))

ggplot(data = tau_post_summary) +
  geom_violin(aes(x = naive_rank, y = mean_rank, group = naive_rank)) +
  geom_jitter(aes(x = naive_rank, y = mean_rank)) +
  #geom_abline(aes(slope = 1, intercept = 0)) +
  labs(x = "Rank within community", 
       y = "Posterior Mean of Tau(X*Beta )")


ggplot(data = tau_post_summary) +
  geom_point(aes(x = PMT_rank, y = mean_rank))

ggplot(data = tau_post_summary) +
  geom_point(aes(x = PMT, y = mu_mean))

ggplot(data = tau_post_summary) +
  geom_point(aes(x = PMT, y = mu_mean)) +
  geom_abline(aes(slope = 1, intercept = 0))


ggplot(data = tau_post_summary[order(tau_post_summary$mean_rank),]) +
  geom_point(aes(x = 1:nrow(tau_post_summary), y= mean_rank)) +
  geom_point(aes(x = 1:nrow(tau_post_summary), y= PMT_rank, colour = community[1:nrow(tau_post_summary)])) 


ggplot(data = tau_post_summary[order(tau_post_summary$mean_rank),]) +
  geom_point(aes(x = PMT_rank, y = mean_rank-PMT_rank,colour = naive_rank)) +
  #geom_point(aes(x = PMT_rank, y = mean_rank-PMT_rank,colour = naive_rank)) +
  #geom_point(aes(x = 1:nrow(tau_post_summary), y = naive_agg), colour = "tomato") +
  #geom_point(aes(x = 1:nrow(tau_post_summary), y = PMT,colour = naive_agg)) +
  #geom_abline(aes(slope = 1, intercept = 0)) +
  labs(x = "Pure PMT rank", y = "Hybrid ranks lower<--->PMT ranks lower \n lower == more need")+
  scale_colour_distiller("Pure CBT Rank\n(w/in village)")



data.frame(postmean =  (apply(tau_post, 1, median)), rank(apply(Tau, 1, mean)))

#posteriors of quality weights - compare to truths
mean(temp$omega_comm); omega_comm_true
mean(temp$omega_micro);omega_micro_true
mean(temp$omega_rank) ;omega_rank_true

apply(temp$beta_rank, 2, mean) ;beta_rank_true
apply(temp$beta_comm, 2, mean) ;beta_comm_true
apply(temp$beta_micro, 2, mean) ;beta_micro_true

###convergence diagnostics------


#do all ESS checks

doESS <- function(x){
  if(!is.null(dim(x))){ #if it's a data frame
    return(apply(x, 2, ESS))
  }else{
    return(ESS(x))
  }
}

effectiv_ss <- lapply(temp, doESS) 

effectiv_ss%>% str()

effectiv_ss
