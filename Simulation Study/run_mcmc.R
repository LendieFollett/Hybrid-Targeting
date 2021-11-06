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
source("Bayes Consensus Ranking/HybridTarget.R")
source("Bayes Consensus Ranking/CBTarget.R")
#parameters for simulation

doESS <- function(x){
  
  if(!is.null(dim(x))){ #if it's a data frame
    return(apply(x, 2, ESS))
  }else{
    return(ESS(x))
  }
}
K = 30 # 50 communities
R = 150 #200/50 rankers per community (may or may not be crossed)
N_CBT = 10*K ## number of ranked/test items (10 per community here)
N_PMT = 300## number of unranked/training items
N_Program = 300 ## number of ranked/test items
P = 6  ## number of covariates


iter_keep = 1000   ## Gibbs sampler kept iterations (post burn-in)
iter_burn =3000   ## Gibbs sampler burn-in iterations 
print_opt = 100  ## print a message every print.opt steps

#simulate data based on parameters
source("Simulation Study/simulate_data.R")


#starting values for random effects
temp_data <- data.frame(y = Y_micro,X_PMT)
form <- formula(paste0("y~-1+", paste0("X", 1:ncol(X_PMT), collapse = "+")))
#gamma_start <- ranef(lmer(form, data = temp_data))[[1]]$`(Intercept)` 
beta_start <-coef(lm(form, data = temp_data))%>%as.vector()
initial_list <- list(beta_rank = beta_start,
                     beta_micro = beta_start)
groups <- rep(c(1:(R/K)), K)
weight_prior_value = c(0.5, 1, 2)
prior_prob_rank=list(rep(1/length(weight_prior_value), length(weight_prior_value)),
                     rep(1/length(weight_prior_value), length(weight_prior_value)),
                     rep(1/length(weight_prior_value), length(weight_prior_value)),
                     rep(1/length(weight_prior_value), length(weight_prior_value)),
                     rep(1/length(weight_prior_value), length(weight_prior_value)))
#Run MCMC for Bayesian Consensus Targeting
temp <- HybridTarget(Tau=Tau, 
                     X_PMT = X_PMT, 
                     X_CBT = X_CBT,
                     X_program = X_Program,
                     X_elite = 7,
                     groups = groups,
                     prior_prob_rank = prior_prob_rank,
                     weight_prior_value =weight_prior_value,
                     Y_micro = Y_micro, #needs to be a matrix, not vector
                     iter_keep = iter_keep,
                     iter_burn = iter_burn,
                     print_opt = print_opt,
                     initial.list = initial_list)
X_program = X_Program
X_elite = 7
initial.list = initial_list
#weight_prior_value = c(0.5, 1, 2)
#prior_prob_rank = list(rep(1/length(weight_prior_value), length(weight_prior_value))) #override if heterogeneous
#groups = rep(1, ncol(Tau))

lapply(temp[c("mu_beta", "beta_rank", "beta_micro")], doESS) 

mu_beta_mean <- apply(temp$mu_beta, 2, mean)
beta_rank_mean <- apply(temp$beta_rank, 2, mean)
beta_micro_mean <- apply(temp$beta_micro, 2, mean)

alpha_mean <- apply(temp$alpha, 2, mean)
plot(alpha_true,alpha_mean)
plot(temp$sigma2_alpha)
mean(sqrt(temp$sigma2_alpha))

plot(temp$omega_micro)

#true coefficients:
beta_rank_true[-1] #intercept is irrelevant
beta_rank_mean[-1]
beta_micro_true
beta_micro_mean
apply(temp$omega_rank, 2, mean)
omega_rank_true

data.frame(parameter = paste("Parameter", c(1:length(mu_beta_mean))),
           mu_beta = mu_beta_mean,
           beta_rank = beta_rank_mean[-1],
           beta_micro=beta_micro_mean[-1])%>%
  melt(id.var = "parameter") %>%
  ggplot() +
  geom_hline(aes(yintercept = 0))+
  geom_line(aes(x = parameter, y = value, colour = variable, group = variable)) +
  geom_point(aes(x = parameter, y = value, colour = variable, group = variable)) +
  coord_flip() +
  labs(x = "Coefficient ", y = "Estimate") +
  scale_colour_brewer("Parameter", palette = "Set1") +
  theme_bw()


poverty_rate <- .3

test_data <- data.frame(X_micro1, consumption = Y_micro1[,1], community=community[1:nrow(X_micro1)]) %>% 
  mutate(hybrid_prediction =         apply(temp$mu, 2, mean),
         hybrid_prediction_noelite =         apply(temp$mu_noelite, 2, mean),
         micro_prediction = X_micro1%*%beta_micro_mean)%>%
  group_by(community) %>%
  mutate(hybrid_rank =rank(hybrid_prediction)/length(community),
         hybrid_noelite_rank =rank(hybrid_prediction_noelite)/length(community),
         pmt_rank =rank(micro_prediction)/length(community),
         consumption_rank = rank(consumption)/length(community)) %>%
  mutate(hybrid_inclusion = hybrid_rank < poverty_rate,
         hybrid_noelite_inclusion = hybrid_noelite_rank < poverty_rate,
         pmt_inclusion = pmt_rank < poverty_rate,
         consumption_inclusion = consumption_rank<poverty_rate) %>%ungroup() %>%
  mutate_at(vars(matches("inclusion")), as.factor)

rbind(confusionMatrix(test_data$hybrid_inclusion, test_data$consumption_inclusion,positive = "TRUE")$byClass,
      confusionMatrix(test_data$hybrid_noelite_inclusion, test_data$consumption_inclusion,positive = "TRUE")$byClass,
      confusionMatrix(test_data$pmt_inclusion, test_data$consumption_inclusion,positive = "TRUE")$byClass) %>%as.data.frame%>%
  mutate(Method = c("Hybrid", "Hybrid Connection Corrected", "PMT"),
         TD = Sensitivity - (1-Specificity)) %>%
  dplyr::select(c(Method,Sensitivity, Specificity, TD))

#posterior summaries of ranks
tau_post <-apply(temp$mu, 1, rank)
tau_post_summary <- data.frame(
  mean_rank =  rank(apply(temp$mu, 2, mean)),#Rank of posterior means of xbeta
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
