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
library(reshape2)
library(gridExtra)
library(LaplacesDemon)

doESS <- function(x){
  
  if(!is.null(dim(x))){ #if it's a data frame
    return(apply(x, 2, ESS))
  }else{
    return(ESS(x))
  }
}

source("Bayes Consensus Ranking/functions.R")
#parameters for simulation

full_data <- read.csv("Empirical Study/alatas.csv") %>%
  select(-c("hhsize_ae"))
#add a community id
#rank is a decimal; needs to be an integer rank
full_data <- full_data%>% 
  mutate(community_id = as.numeric(factor(interaction(village, province, district, subdistrict))))%>%
  group_by(village, province, district, subdistrict)%>%
         mutate(rank = ifelse(is.na(rank), NA, floor(rank(rank)))) %>%ungroup %>%
  subset(community == 1 | hybrid == 1)
set.seed(572319852)
train_idx <-which(full_data$community_id %in% sample(unique(full_data$community_id), replace=FALSE, length(unique(full_data$community_id))*.7))

#groups of x variables
m1 <- c("connected","hhsize","hhage","hhmale","hhmarried",
        "hheduc2","hheduc3","hheduc4",
        "age04","higheduc2","higheduc3","higheduc4","depratio")
m2.1 <- c("pcfloor", "tfloor","twall", "toilet","water","lighting", "troof",
          "fcook","house", "ac","computer","radio","tv", "dvd","satellite", 
          "gas", "refrigerator", "bicycle", "motorcycle", "auto", "hp", 
          "jewelry","chicken","cow")
m2 <- c(m1, m2.1)
m3.1 <- c("credit","hhsector1", "hhsector2","hhsector3",
          "formal","informal", "eschild","jschild","sschild")
m3 <- c(m2,m3.1) #full collection
test_data <- full_data[-train_idx,] #%>% subset(community == 1)
train_data <- full_data[train_idx,] #%>% subset(community == 0)

X_micro0 <- cbind(1, train_data[,m3]%>%apply(2, function(x){(x - mean(x))/sd(x)})) 
X_micro1 <- cbind(1, test_data[,m3]%>%apply(2, function(x){(x - mean(x))/sd(x)}))

Y_micro <- as.matrix(log(train_data$consumption))
Y_micro <- apply(Y_micro, 2, function(x){(x - mean(x))/sd(x)})

R = test_data %>% group_by(village, province, district, subdistrict) %>% summarise(n = length(cow))%>%ungroup() %>%nrow
M = 1  ## just consumption
N0 = train_data %>%nrow
N1 = test_data %>%nrow
P = ncol(X_micro0)-1 #(-1 since i've included the intercept)


#starting values for random effects
temp_data <- data.frame(y = Y_micro,
                         X_micro0)
form <- formula(paste0("y~", paste0(colnames(X_micro0), collapse = "+")))
#gamma_start <- ranef(lmer(form, data = temp_data))[[1]]$`(Intercept)` 
beta_start <-coef(lm(form, data = temp_data))%>%as.vector()
initial_list <- list(#gamma_rank = gamma_start,
  beta_rank = c(0,beta_start[-1]),
  beta_micro = beta_start,
  mu_beta = beta_start[-1])

#create rank matrix: one column per 'ranker' (community)
Tau <- array(NA, dim = c(nrow(test_data), R))
j = 0
for ( idx in unique(test_data$community_id)){ #loop over columns
  j = j + 1
  Tau[test_data$community_id == idx,j] <- test_data$rank[test_data$community_id == idx]
}


iter_keep = 1000   ## Gibbs sampler kept iterations (post burn-in)
iter_burn =1000   ## Gibbs sampler burn-in iterations 
print_opt = 100  ## print a message every print.opt steps


#Run MCMC for Bayesian Consensus Targeting
temp <- BCTarget(Tau=Tau, 
                 X_micro0 = X_micro0, 
                 X_micro1 = X_micro1,
                 X_elite = "connected",
                 Y_micro = Y_micro, #needs to be a matrix, not vector
                 prior_prob_rank = c(.2,.2,.6),
                 prior_prob_micro = c(.6,.2,.2),
                 iter_keep = iter_keep,
                 iter_burn = iter_burn,
                  print_opt = print_opt,
                 initial.list = initial_list)


lapply(temp[c("mu_beta", "beta_rank", "beta_micro")], doESS) 

mu_beta_mean <- apply(temp$mu_beta, 2, mean)
beta_rank_mean <- apply(temp$beta_rank, 2, mean)
beta_micro_mean <- apply(temp$beta_micro, 2, mean)

data.frame(parameter = colnames(X_micro0)[-1],
           mu_beta = mu_beta_mean,
           beta_rank = beta_rank_mean[-1],
           beta_micro=beta_micro_mean[-1])%>%
  melt(id.var = "parameter") %>%
  mutate(parameter = factor(parameter, levels = colnames(X_micro0)[-1][order(mu_beta_mean)]))%>%
  ggplot() +
  geom_hline(aes(yintercept = 0))+
  geom_line(aes(x = parameter, y = value, colour = variable, group = variable)) +
  geom_point(aes(x = parameter, y = value, colour = variable, group = variable)) +
  coord_flip() +
  labs(x = "Coefficient ", y = "Estimate") +
  scale_colour_brewer("Parameter", palette = "Set1") +
  theme_bw()
ggsave("coefficients.pdf", width = 6, height = 10)

#get back on log(consumption scale) --->sigma*predicted + mu
test_data$hybrid_prediction <-         apply(temp$mu, 2, mean)#apply(temp$mu*sd(log(train_data$consumption)) + mean(log(train_data$consumption)),2,mean)
test_data$hybrid_prediction_noelite <- apply(temp$mu_noelite, 2, mean)#apply(temp$mu_noelite*sd(log(train_data$consumption)) + mean(log(train_data$consumption)),2,mean)
#beta_start is the OLS estimate of beta
test_data$micro_prediction <- (X_micro1[,-1]%*%beta_start[-1])#*sd(log(train_data$consumption)) + mean(log(train_data$consumption))


poverty_rate <- .3

test_data <- test_data %>% group_by(village, province, district, subdistrict) %>%
  mutate(hybrid_rank =rank(hybrid_prediction)/length(village),
         hybrid_noelite_rank =rank(hybrid_prediction_noelite)/length(village),
         pmt_rank =rank(micro_prediction)/length(village),
         consumption_rank = rank(consumption)/length(village),
         cbt_rank = rank/length(village)) %>%
  mutate(hybrid_inclusion = hybrid_rank < poverty_rate,
         hybrid_noelite_inclusion = hybrid_noelite_rank < poverty_rate,
         pmt_inclusion = pmt_rank < poverty_rate,
         consumption_inclusion = consumption_rank<poverty_rate) %>%ungroup() %>%
  mutate(Z_mean = apply(temp$Z, 2, mean))%>%
  mutate_at(vars(matches("inclusion")), as.factor)

ggplot(data = test_data) + 
  geom_boxplot(aes(x = connected, y = consumption_rank,group=connected, colour = elite))+
  geom_jitter(aes(x = connected, y = consumption_rank,group=elite,colour = elite))


qplot(Z_mean, cbt_rank, colour = pmt_rank,data = test_data)

ggplot(data = test_data) +
  geom_point(aes(cbt_rank,Z_mean, colour = pmt_rank)) +
               geom_line(aes(x = cbt_rank, y = Z_mean, group = community_id), alpha = I(.3))


library(caret)

rbind(confusionMatrix(test_data$hybrid_inclusion, test_data$consumption_inclusion,positive = "TRUE")$byClass,
confusionMatrix(test_data$hybrid_noelite_inclusion, test_data$consumption_inclusion,positive = "TRUE")$byClass,
confusionMatrix(test_data$pmt_inclusion, test_data$consumption_inclusion,positive = "TRUE")$byClass) %>%as.data.frame%>%
  mutate(Method = c("Hybrid", "Hybrid Connection Corrected", "PMT"),
         TD = Sensitivity - (1-Specificity)) %>%
  dplyr::select(c(Method,Sensitivity, Specificity, TD))

ggplot(data = test_data) + 
 geom_point(aes(log(consumption),cbt_rank, colour = cbt_rank),width = .025, height = .025) #+
 # geom_hline(aes(yintercept = poverty_rate)) +
 # geom_vline(aes(xintercept = poverty_rate))

p1 <- ggplot(test_data) +
  geom_point(aes(x = hybrid_prediction, y = consumption)) +
  geom_abline(aes(intercept = 0, slope = 1))+
  ggtitle("Hybrid") +scale_y_continuous(limits = c(-2,2))








