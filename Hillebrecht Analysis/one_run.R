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
source("Bayes Consensus Ranking/HybridTarget.R")
source("Bayes Consensus Ranking/CBTarget.R")
#parameters for simulation


poverty_rate <- .3
iter_keep = 1500   ## Gibbs sampler kept iterations (post burn-in)
iter_burn =1500   ## Gibbs sampler burn-in iterations 
print_opt = 100  ## print a message every print.opt steps


full_data <- read.csv("Hillebrecht Analysis/hillebrecht.csv") %>%
  group_by(community)%>%
  mutate(informant1 = ifelse(is.na(informant1), NA, floor(rank(-informant1))),
         informant2 = ifelse(is.na(informant2), NA, floor(rank(-informant2))),
         informant3 = ifelse(is.na(informant3), NA, floor(rank(-informant3))))%>%ungroup%>%  arrange(community)
#x variables to include in model
m_num <- c("rooms")

m_bin <- colnames(full_data)[which(colnames(full_data)=="floor"):which(colnames(full_data)=="pig")]
m3 <- c(m_num, m_bin)
full_data <- full_data %>%mutate_at(m_num, function(x){(x - mean(x))/(2*sd(x))}) 

#50% of the full data is surveyed for PMT. get both X and y=consumption

PMT_idx <-which(full_data$year == 2008) #training data is all of 2008 data
PMT_data <- full_data[PMT_idx,] #%>% subset(community == 0)

full_data_left <- full_data[-PMT_idx,]

CBT_ncomm <- 25#c(5,10,15,25)

whats_left <- unique(full_data_left$community) #all 2009 data
Program_idx <- which(full_data_left$community %in% sample(whats_left, 
                                                        replace=FALSE, 
                                                        length(whats_left)*0.5)) #this will leave 26 communities for the CBT
whats_left <- unique(full_data_left$community[-c( Program_idx)])
length(whats_left)
whats_left
CBT_idx <- which(full_data_left$community %in% sample(whats_left, 
                                                    replace=FALSE, 
                                                    CBT_ncomm))


CBT_data <- full_data_left[CBT_idx,] #this is a subset of the program data!

while(any(apply(CBT_data[,m3], 2, var) == 0)){ #have to do to deal with complete separation and ML estimation of logistic regression (#shouldadonebayes)
  whats_left <- unique(full_data_left$community)
  Program_idx <- which(full_data_left$community %in% sample(whats_left, 
                                                       replace=FALSE, 
                                                       length(whats_left)*0.5))
  whats_left <- unique(full_data_left$community[-c(Program_idx)])
  CBT_idx <- which(full_data_left$community %in% sample(whats_left, 
                                                   replace=FALSE, 
                                                   CBT_ncomm))
  CBT_data <- full_data_left[CBT_idx,] #this (can be) a subset of the program data!

}
Program_data <- full_data[Program_idx,]


X_PMT <-     cbind(1,PMT_data[,m3]) %>%as.matrix()#cbind(1, PMT_data[,m3]%>%apply(2, function(x){(x - mean(x))/sd(x)})) 
X_CBT <-     cbind(1,CBT_data[,m3]) %>%as.matrix()
X_program <- cbind(1,Program_data[,m3]) %>%as.matrix()
Y_micro <- as.matrix(log(PMT_data$consumption + .1)) #ROUGH FIX - USE OTHER TRANSFORMATION
Y_micro <- apply(Y_micro, 2, function(x){(x - mean(x))/sd(x)})

R = 3*length(unique(CBT_data$community))#(CBT_data %>% group_by(community) %>% summarise(n = length(floor))%>%ungroup() %>%nrow)
M = 1  ## just consumption
N_PMT = PMT_data %>%nrow
N_CBT = CBT_data %>%nrow
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

#create rank matrix: one column per 'ranker' (within community)
Tau <- array(NA, dim = c(N_CBT, R))
j = 0
for ( idx in unique(CBT_data$community)){ #loop over columns
  for (infmt in c("informant1", "informant2", "informant3")){
  j = j + 1
  Tau[CBT_data$community == idx,j] <- floor(rank(pull(CBT_data[CBT_data$community == idx,], infmt)))
  }
}

#Run MCMC for Bayesian Consensus Targeting
groups <- rep(1, R) #there are 3 rankers per comm. but they aren't crossed
weight_prior_value = c(0.5, 1, 2)
prior_prob_rank=list(rep(1/length(weight_prior_value), length(weight_prior_value)))

temp <- HybridTarget(Tau=Tau, 
                     X_PMT = X_PMT, 
                     X_CBT = X_CBT,
                     X_program = X_program,
                     X_elite = NULL,
                     groups = groups,
                     weight_prior_value = weight_prior_value,
                     prior_prob_rank=prior_prob_rank,
                     Y_micro = Y_micro, #needs to be a matrix, not vector
                     iter_keep = iter_keep,
                     iter_burn = iter_burn,
                     print_opt = print_opt,
                     initial.list = initial_list)
X_program = X_program
X_elite = NULL
initial.list = initial_list

CBtemp <- CBTarget(Tau=Tau, 
                   X_CBT = X_CBT,
                   X_program = X_program,
                   X_elite = NULL,
                   iter_keep =iter_keep,
                   iter_burn = iter_burn,
                   print_opt = print_opt,
                   initial.list = initial_list)


#HYBRID-BASED PREDICTION
Program_data$hybrid_prediction <-        apply(temp$mu, 2, mean)
#CBT SCORE-BASED PREDICTION
Program_data$cbt_model_prediction <- apply(CBtemp$mu, 2, mean)#(X_program%*%CB_beta_rank_mean)

#BAYESIAN LOGISTIC REGRESSION CBT-BASED PREDICTION
#how to define for multiple rankers?
#Program_data$CBT_LR_prediction<- -predict(lr, as.data.frame(X_program[,-1])) #(LRF NEEDS TO CHANGE TO) logistic regression

#OLS-BASED PMT PREDICTION
Program_data$PMT_prediction <- (X_program[,-1]%*%beta_start[-1])#beta_start is the OLS estimate of beta

Program_data <- Program_data%>%group_by(community) %>%
  mutate(
    hybrid_rank =rank(hybrid_prediction)/length(community),
    pmt_rank =rank(PMT_prediction)/length(community),
    cbt_model_rank = rank(cbt_model_prediction)/length(community),
    consumption_rank = rank(consumption)/length(community)
    #CBT_LR_rank = rank(CBT_LR_prediction)/length(village)
    ) %>%
  mutate(#hybrid_inclusion = hybrid_rank <= poverty_rate,
    hybrid_inclusion = hybrid_rank <= poverty_rate,
    pmt_inclusion = pmt_rank <= poverty_rate,
    consumption_inclusion = consumption_rank<=poverty_rate,
    cbt_model_inclusion = cbt_model_rank<=poverty_rate,
    #CBT_LR_inclusion = CBT_LR_rank<=poverty_rate,
    cbt_inclusion = ifelse(treated == 1, TRUE, FALSE)) %>%ungroup() %>%
  mutate_at(vars(matches("inclusion")), as.factor)

CBT_prop <- .3

rbind(
  confusionMatrix(Program_data$hybrid_inclusion, Program_data$cbt_inclusion,positive = "TRUE")$byClass,
  confusionMatrix(Program_data$cbt_model_inclusion,      Program_data$cbt_inclusion,positive = "TRUE")$byClass,
  confusionMatrix(Program_data$pmt_inclusion,            Program_data$cbt_inclusion,positive = "TRUE")$byClass) %>%as.data.frame%>%
  mutate(Method = c( "Hybrid Score","CBT Score", "PMT OLS"),
         CBT_prop = CBT_prop,
         TD = Sensitivity - (1-Specificity)) %>%
  dplyr::select(c(Method,CBT_prop,Sensitivity, Specificity, TD))


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


