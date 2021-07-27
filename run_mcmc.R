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
library(caret)
library(parallel)
detectCores(logical=FALSE)

doESS <- function(x){
  
  if(!is.null(dim(x))){ #if it's a data frame
    return(apply(x, 2, ESS))
  }else{
    return(ESS(x))
  }
}

source("Bayes Consensus Ranking/functions.R")
#parameters for simulation


poverty_rate <- .3
iter_keep = 1000   ## Gibbs sampler kept iterations (post burn-in)
iter_burn =1000   ## Gibbs sampler burn-in iterations 
print_opt = 100  ## print a message every print.opt steps


full_data <- read.csv("Empirical Study/alatas.csv") %>%
  select(-c("hhsize_ae")) %>% arrange(village, province, district, subdistrict)%>% 
  mutate(community_id = as.numeric(factor(interaction(village, province, district, subdistrict))))%>%
  group_by(village, province, district, subdistrict)%>%
  mutate(prop_rank = rank) %>%
  mutate(rank = ifelse(is.na(rank), NA, floor(rank(rank)))) %>%ungroup
#x variables to include in model
m3 <- c("connected","hhsize","hhage","hhmale","hhmarried",
        "hheduc2","hheduc3","hheduc4",
        "age04","higheduc2","higheduc3","higheduc4","depratio","pcfloor",
        "tfloor","twall", "toilet","water","lighting", "troof",
        "fcook","house", "computer","radio","tv", "dvd","satellite", #LRF REMOVED AC FOR NOW
        "gas", "refrigerator", "bicycle", "motorcycle", "auto", "hp", 
        "jewelry","chicken","cow", "credit","hhsector1", "hhsector2","hhsector3",
        "formal","informal", "eschild","jschild","sschild")

#50% of the full data is surveyed for PMT. get both X and y=consumption
#set.seed(572319852)
PMT_idx <-which(full_data$pmt == 1)#which(full_data$community_id %in% sample(unique(full_data$community_id), replace=FALSE,  length(unique(full_data$community_id))*.5))
#Note: the hh index for program is everything else, e.g. , full_data[-PMT_idx,]
#a subset of the program data is CBT






#parallelized across CBT proportions via mcapply
CBT_prop_list <- c(.05,.1, .25)  
 results <-  mclapply(CBT_prop_list, function(CBT_prop){
   i <- 0
   r <- list()
  for(rep in c(1:5)){
    print(paste("***********Rep ", rep," of CBT proportion ", CBT_prop, "**************"))
    i = i + 1
    
Program_idx <- which(full_data$community_id %in% sample(unique(full_data$community_id[- PMT_idx]), 
                                                    replace=FALSE, 
                                                    length(unique(full_data$community_id[-PMT_idx]))*0.5))
CBT_idx <-which(full_data$community_id %in% sample(unique(full_data$community_id[-c(PMT_idx, Program_idx)]), 
                                                   replace=FALSE, 
                                                   length(unique(full_data$community_id[-c(PMT_idx, Program_idx)]))*CBT_prop))

CBT_data <- full_data[CBT_idx,] #this is a subset of the program data!
PMT_data <- full_data[PMT_idx,] #%>% subset(community == 0)
Program_data <- full_data[Program_idx,]


X_PMT <- cbind(1, PMT_data[,m3]%>%apply(2, function(x){(x - mean(x))/sd(x)})) 
X_CBT <- cbind(1, CBT_data[,m3]%>%apply(2, function(x){(x - mean(x))/sd(x)}))
X_program <- cbind(1, Program_data[,m3]%>%apply(2, function(x){(x - mean(x))/sd(x)}))

Y_micro <- as.matrix(log(PMT_data$consumption))
Y_micro <- apply(Y_micro, 2, function(x){(x - mean(x))/sd(x)})

R = CBT_data %>% group_by(village, province, district, subdistrict) %>% summarise(n = length(cow))%>%ungroup() %>%nrow
M = 1  ## just consumption
N0 = PMT_data %>%nrow
N1 = CBT_data %>%nrow
P = ncol(PMT_data)-1 #(-1 since i've included the intercept)


#starting values for random effects
temp_data <- data.frame(y = Y_micro,
                        X_PMT)
form <- formula(paste0("y~", paste0(colnames(X_PMT), collapse = "+")))
#gamma_start <- ranef(lmer(form, data = temp_data))[[1]]$`(Intercept)` 
beta_start <-coef(lm(form, data = temp_data))%>%as.vector()
initial_list <- list(#gamma_rank = gamma_start,
  beta_rank = c(0,beta_start[-1]),
  beta_micro = beta_start,
  mu_beta = beta_start[-1])

#create rank matrix: one column per 'ranker' (community)
Tau <- array(NA, dim = c(nrow(CBT_data), R))
j = 0
for ( idx in unique(CBT_data$community_id)){ #loop over columns
  j = j + 1
  Tau[CBT_data$community_id == idx,j] <- CBT_data$rank[CBT_data$community_id == idx]
}

#Run MCMC for Bayesian Consensus Targeting

Hybridtemp <- HybridTarget(Tau=Tau, 
                 X_PMT = X_PMT, 
                 X_CBT = X_CBT,
                 X_program = X_program,
                 X_elite = "connected",
                 Y_micro = Y_micro, #needs to be a matrix, not vector
                 prior_prob_rank = c(.025, .025,.95),
                 prior_prob_micro = c(.95,.025, .025),
                 iter_keep = iter_keep,
                 iter_burn = iter_burn,
                  print_opt = print_opt,
                 initial.list = initial_list)

CBtemp <- CBTarget(Tau=Tau, 
                 X_CBT = X_CBT,
                 X_program = X_program,
                 X_elite = "connected",
                 prior_prob_rank = c(.025, .025,.95),
                 iter_keep =iter_keep,
                 iter_burn = iter_burn,
                 print_opt = print_opt,
                 initial.list = initial_list)


#lapply(Hybridtemp[c("mu_beta", "beta_rank", "beta_micro")], doESS) 

Hybrid_mu_beta_mean <- apply(Hybridtemp$mu_beta, 2, mean)
Hybrid_beta_rank_mean <- apply(Hybridtemp$beta_rank, 2, mean)
Hybrid_beta_micro_mean <- apply(Hybridtemp$beta_micro, 2, mean)

CB_beta_rank_mean <- apply(CBtemp$beta_rank, 2, mean)

data.frame(parameter = m3,
           mu_beta = Hybrid_mu_beta_mean,
           beta_rank = Hybrid_beta_rank_mean[-1],
           beta_micro=Hybrid_beta_micro_mean[-1])%>%
  melt(id.var = "parameter") %>%
  mutate(parameter = factor(parameter, levels = colnames(X_CBT)[-1][order(Hybrid_mu_beta_mean)]))%>%
  ggplot() +
  geom_hline(aes(yintercept = 0))+
  geom_line(aes(x = parameter, y = value, colour = variable, group = variable)) +
  geom_point(aes(x = parameter, y = value, colour = variable, group = variable)) +
  coord_flip() +
  labs(x = "Coefficient ", y = "Estimate") +
  scale_colour_brewer("Parameter", palette = "Set1") +
  theme_bw()
ggsave(paste0("coefficients_",CBT_prop, "_", i, ".pdf"), width = 6, height = 10)

m <- lm(prop_rank ~ ., data = CBT_data[,c("prop_rank", m3)])

#get back on log(consumption scale) --->sigma*predicted + mu
#HYBRID-BASED PREDICTION
Program_data$hybrid_prediction <-        (X_program[,-c(1)]%*%Hybrid_beta_rank_mean[-c(1)]) #apply(Hybridtemp$mu, 2, mean)
Program_data$hybrid_prediction_noelite <-(cbind(0,X_program[,-c(1,2)])%*%Hybrid_beta_rank_mean[-c(1)]) #apply(Hybridtemp$mu_noelite, 2, mean)

#CBT SCORE-BASED PREDICTION
Program_data$cbt_model_prediction <- (X_program[,-c(1)]%*%CB_beta_rank_mean[-c(1)])

#BAYESIAN LOGISTIC REGRESSION CBT-BASED PREDICTION
Program_data$CBT_noshrink_prediction<- predict(m, as.data.frame(X_program[,-1])) #(LRF NEEDS TO CHANGE TO) logistic regression

#OLS-BASED PMT PREDICTION
Program_data$PMT_prediction <- (X_program[,-1]%*%beta_start[-1])#beta_start is the OLS estimate of beta

Program_data <- Program_data%>%group_by(village, province, district, subdistrict) %>%
  mutate(#hybrid_rank =rank(hybrid_prediction)/length(village),
         hybrid_noelite_rank =rank(hybrid_prediction_noelite)/length(village),
         pmt_rank =rank(PMT_prediction)/length(village),
         cbt_model_rank = rank(cbt_model_prediction)/length(village),
         consumption_rank = rank(consumption)/length(village),
         cbt_rank = rank/length(village),
         CBT_noshrink_rank = rank(CBT_noshrink_prediction)/length(village)) %>%
  mutate(#hybrid_inclusion = hybrid_rank <= poverty_rate,
         hybrid_noelite_inclusion = hybrid_noelite_rank <= poverty_rate,
         pmt_inclusion = pmt_rank <= poverty_rate,
         consumption_inclusion = consumption_rank<=poverty_rate,
         cbt_model_inclusion = cbt_model_rank<=poverty_rate,
         CBT_noshrink_inclusion = CBT_noshrink_rank<=poverty_rate,
         cbt_inclusion = cbt_rank <= poverty_rate) %>%ungroup() %>%
  mutate_at(vars(matches("inclusion")), as.factor)



r[[i]] <- rbind(#confusionMatrix(Program_data$hybrid_inclusion,   Program_data$cbt_inclusion,positive = "TRUE")$byClass,
confusionMatrix(Program_data$hybrid_noelite_inclusion, Program_data$cbt_inclusion,positive = "TRUE")$byClass,
confusionMatrix(Program_data$cbt_model_inclusion,      Program_data$cbt_inclusion,positive = "TRUE")$byClass,
confusionMatrix(Program_data$pmt_inclusion,            Program_data$cbt_inclusion,positive = "TRUE")$byClass,
confusionMatrix(Program_data$CBT_noshrink_inclusion,            Program_data$cbt_inclusion,positive = "TRUE")$byClass) %>%as.data.frame%>%
  mutate(Method = c( "Hybrid Model","CBT Score Bayesian", "PMT OLS", "CBT OLS"),
         CBT_prop = CBT_prop,
         TD = Sensitivity - (1-Specificity)) %>%
  dplyr::select(c(Method,CBT_prop,Sensitivity, Specificity, TD))


}
   return(r)
  }, mc.cores = length(CBT_prop_list))

all_results_1 <- unlist(results, recursive = FALSE)
all_results <- do.call("rbind", all_results_1)


all_results %>%melt(id.var = c("Method", "CBT_prop")) %>%
  ggplot() +geom_boxplot(aes(x = Method, y = value,colour = Method, group = interaction(Method, CBT_prop))) + 
  geom_point(aes(x = Method, y = value,colour = Method, group = interaction(Method, CBT_prop))) + 
  facet_grid(variable~CBT_prop, scales = "free") +theme(axis.text.x = element_text(angle = 45))
ggsave("results_with_con.pdf")

qplot(1:1000,temp$beta_rank[,21]) + 
  geom_point(aes(1:1000,temp$beta_micro[,21]), colour = "red")+
  geom_line(aes(1:1000,temp$mu_beta[,20]), alpha = I(.4))



p1 <-ggplot(data = test_data) + 
  geom_boxplot(aes(x = connected, y = consumption_rank,group=connected, colour = elite))+
  geom_jitter(aes(x = connected, y = consumption_rank,group=elite,colour = elite))


p2 <-qplot(Z_mean, cbt_rank, colour = pmt_rank,data = test_data)

p3 <-ggplot(data = test_data) +
  geom_point(aes(cbt_rank,Z_mean, colour = pmt_rank)) +
  geom_line(aes(x = cbt_rank, y = Z_mean, group = community_id), alpha = I(.3))

ggplot(data = test_data) + 
 geom_point(aes(log(consumption),cbt_rank, colour = cbt_rank),width = .025, height = .025) #+
 # geom_hline(aes(yintercept = poverty_rate)) +
 # geom_vline(aes(xintercept = poverty_rate))

p1 <- ggplot(test_data) +
  geom_point(aes(x = hybrid_prediction, y = consumption)) +
  geom_abline(aes(intercept = 0, slope = 1))+
  ggtitle("Hybrid") +scale_y_continuous(limits = c(-2,2))

