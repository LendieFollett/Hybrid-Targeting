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
iter_keep = 1500   ## Gibbs sampler kept iterations (post burn-in)
iter_burn =1500   ## Gibbs sampler burn-in iterations 
print_opt = 100  ## print a message every print.opt steps


full_data <- read.csv("Empirical Study/alatas.csv") %>%
  dplyr::select(-c("hhsize_ae")) %>% arrange(village, province, district, subdistrict)%>% 
  mutate(community_id = as.numeric(factor(interaction(village, province, district, subdistrict))))%>%
  group_by(village, province, district, subdistrict)%>%
  mutate(prop_rank = rank) %>%
  mutate(rank = ifelse(is.na(rank), NA, floor(rank(rank)))) %>%ungroup

#x variables to include in model
m_num <- c("hhsize","hhage",
           "age04","depratio","pcfloor",
           "eschild","jschild","sschild")

m_bin <- c("connected","hhmale","hhmarried",
           "hheduc2","hheduc3","hheduc4",
           "higheduc2","higheduc3","higheduc4",
           "tfloor","twall", "toilet","water","lighting", "troof",
           "fcook","house", "computer","radio","tv", "dvd","satellite", "ac",
           "gas", "refrigerator", "bicycle", "motorcycle", "auto", "hp", 
           "jewelry","chicken","cow", "credit","hhsector1", "hhsector2","hhsector3",
           "formal", "informal")# lrf removed informal for now
m3 <- c(m_num, m_bin, "hhage2", "hhsize2")

#50% of the full data is surveyed for PMT. get both X and y=consumption
#set.seed(572319852)
PMT_idx <-which(full_data$pmt == 1)#which(full_data$community_id %in% sample(unique(full_data$community_id), replace=FALSE,  length(unique(full_data$community_id))*.5))
#Note: the hh index for program is everything else, e.g. , full_data[-PMT_idx,]
#a subset of the program data is CBT

full_data <- full_data %>%mutate_at(m_num, function(x){(x - mean(x))/(2*sd(x))}) %>%
  mutate(hhage2 = hhage^2,
         hhsize2 = hhsize^2)


CBT_prop <- 0.1


whats_left <- unique(full_data$community_id[- PMT_idx])
Program_idx <- which(full_data$community_id %in% sample(whats_left, 
                                                        replace=FALSE, 
                                                        length(whats_left)*0.5))
whats_left <- unique(full_data$community_id[-c(PMT_idx, Program_idx)])
CBT_idx <- which(full_data$community_id %in% sample(whats_left, 
                                                    replace=FALSE, 
                                                    length(whats_left)*CBT_prop))


CBT_data <- full_data[CBT_idx,] #this is a subset of the program data!
PMT_data <- full_data[PMT_idx,] #%>% subset(community == 0)

while(any(apply(CBT_data[,m3], 2, var) == 0)){ #have to do to deal with complete separation and ML estimation of logistic regression (#shouldadonebayes)
  Program_idx <- which(full_data$community_id %in% sample(unique(full_data$community_id[- PMT_idx]), 
                                                          replace=FALSE, 
                                                          length(unique(full_data$community_id[-PMT_idx]))*0.5))
  CBT_idx <-which(full_data$community_id %in% sample(unique(full_data$community_id[-c(PMT_idx, Program_idx)]), 
                                                     replace=FALSE, 
                                                     length(unique(full_data$community_id[-c(PMT_idx, Program_idx)]))*CBT_prop))
  
  
  CBT_data <- full_data[CBT_idx,] #this is a subset of the program data!
  PMT_data <- full_data[PMT_idx,] #%>% subset(community == 0)
}

Program_data <- full_data[Program_idx,]


X_PMT <-     cbind(1,PMT_data[,m3]) %>%as.matrix()#cbind(1, PMT_data[,m3]%>%apply(2, function(x){(x - mean(x))/sd(x)})) 
X_CBT <-     cbind(1,CBT_data[,m3]) %>%as.matrix()
X_program <- cbind(1,Program_data[,m3]) %>%as.matrix()
X_program_noelite <- X_program
X_program_noelite[,"connected"] <- 0
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

temp <- HybridTarget(Tau=Tau, 
                           X_PMT = X_PMT, 
                           X_CBT = X_CBT,
                           X_program = X_program,
                           X_elite = "connected",
                           Y_micro = Y_micro, #needs to be a matrix, not vector
                           prior_prob_rank = c(1,1,1)/3,
                           prior_prob_micro = c(1,1,1)/3,
                           iter_keep = iter_keep,
                           iter_burn = iter_burn,
                           print_opt = print_opt,
                           initial.list = initial_list)


apply(temp$omega_rank, 2, mean); apply(Hybridtemp$omega_rank, 2, var)

plot(temp$con)
plot(temp$Z[,2])

i = i + 1
plot(temp$beta_micro[,i])

lapply(temp, doESS)


CBtemp <- CBTarget(Tau=Tau, 
                   X_CBT = X_CBT,
                   X_program = X_program,
                   X_elite = "connected",
                   prior_prob_rank = c(1,1,1)/3,
                   iter_keep =iter_keep,
                   iter_burn = iter_burn,
                   print_opt = print_opt,
                   initial.list = initial_list)


lapply(CBtemp, doESS)


