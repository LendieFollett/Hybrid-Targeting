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
sourceCpp("functions.cpp")
source("Bayes Consensus Ranking/functions.R")
source("Bayes Consensus Ranking/HybridTarget.R")
source("Bayes Consensus Ranking/CBTarget.R")

iter_keep = 2000   ## Gibbs sampler kept iterations (post burn-in)
iter_burn = 2000   ## Gibbs sampler burn-in iterations 
print_opt = 100  ## print a message every print.opt steps


full_data <- read.csv("Hillebrecht Analysis/hillebrecht.csv") %>%
  group_by(community, year)%>%
  mutate(informant1 = ifelse(is.na(informant1), NA, floor(rank(-informant1))),
         informant2 = ifelse(is.na(informant2), NA, floor(rank(-informant2))),
         informant3 = ifelse(is.na(informant3), NA, floor(rank(-informant3))),
         treat_rate = sum(treated)/length(treated)) %>% ungroup%>%  arrange(community)
#x variables to include in model
m_num <- c("rooms")

m_bin <- colnames(full_data)[which(colnames(full_data)=="floor"):which(colnames(full_data)=="pig")]
m3 <- c(m_num, m_bin)

#50% of the full data is surveyed for PMT. get both X and y=consumption
#set.seed(572319852)
full_data <- full_data %>%mutate_at(m_num, function(x){(x - mean(x))/(2*sd(x))})

PMT_idx <-which(full_data$year == 2008) #training data is all of 2008 data
full_data_left <- full_data[-PMT_idx,]

CBT_ncomm = 26    
    
    #RANDOM SAMPLES OF CBT DATA, PMT DATA, AND PROGRAM (testing) DATA
    whats_left <- unique(full_data_left$community) #communities not in PMT
    samps <- data.frame(community = whats_left,
                        samp =     rep(c("CBT2", "Program", "NA"), c( CBT_ncomm, 25, max(51-10-CBT_ncomm - 25, 0)))[sample.int(length(whats_left))])
    
    CBT2_data <- full_data_left %>%subset(community %in% samps$community[samps$samp == "CBT2"])
    Program_data <- full_data_left %>%subset(community %in% samps$community[samps$samp == "Program"])    
    PMT_data <- full_data[PMT_idx,] #%>% subset(community == 0)
    
    X_PMT <-     cbind(1,PMT_data[,m3]) %>%as.matrix()#cbind(1, PMT_data[,m3]%>%apply(2, function(x){(x - mean(x))/sd(x)})) 
    X_CBT2 <-     cbind(1,CBT2_data[,m3]) %>%as.matrix()
    X_program <- cbind(1,Program_data[,m3]) %>%as.matrix()
    Y_micro <- as.matrix(log(PMT_data$consumption + .1))
    Y_micro <- apply(Y_micro, 2, function(x){(x - mean(x))/sd(x)})
    
    R2 = 3*(CBT2_data %>% group_by(community) %>% summarise(n = length(floor))%>%ungroup() %>%nrow)
    
    which_noelite <- which(colnames(X_CBT2) == "minority") #NOTE THIS INDEX INCLUDES THE FIRST POSITION OF INTERCEPT
    
    #starting values 
    temp_data <- data.frame(Y_micro = Y_micro,
                            X_PMT)
    form <- formula(paste0("Y_micro~", paste0(colnames(X_PMT), collapse = "+")))
    
    PMT_beta_start <-coef(lm(form, data = temp_data))%>%as.vector()
    PMT_beta_start[is.na(PMT_beta_start)] <- 0
    
    temp_data <- data.frame(rank = apply(CBT2_data[,c("informant1", "informant2", "informant3")], 1, mean) ,
                            X_CBT2) %>%
      mutate(rank = (rank - mean(rank))/sd(rank))
    
    form <- formula(paste0("rank~", paste0(colnames(X_CBT2), collapse = "+")))
    
    CBT_beta_start <-coef(lm(form, data = temp_data))%>%as.vector()
    CBT_beta_start[is.na(CBT_beta_start)] <- 0
    
    mu_beta_start <- apply(cbind(PMT_beta_start, CBT_beta_start), 1, mean) %>%c()
    

    
    initial_list_noelite <- list(
      beta_rank = c(0,CBT_beta_start[-1]),
      beta_micro = PMT_beta_start,
      mu_beta = mu_beta_start[-1])
    
    
    #create rank matrix: one column per 'ranker' (community)
    
    Tau2 <- array(NA, dim = c(nrow(CBT2_data), R2))
    j = 0
    for ( idx in unique(CBT2_data$community)){ #loop over columns
      for (infmt in c("informant1", "informant2", "informant3")){
        j = j + 1
        temp = pull(CBT2_data[CBT2_data$community == idx,], infmt)  
        if (infmt == "informant3"){
        Tau2[CBT2_data$community == idx,j] <- temp[sample(1:length(temp), size = length(temp))]  
        }else{
        Tau2[CBT2_data$community == idx,j] <- temp
        }
      }
    }
    
    groups = rep(c(1:3), R2/3)
    prior_prob_rank = list(c(1,1,1)/3,
                           c(1,1,1)/3,
                           c(1,1,1)/3)
    #Run MCMC for Bayesian Consensus Targeting
    
    #Run MCMC for Bayesian Consensus Targeting - WITH CORRECTION
    Hybridtemp_noelite <- HybridTarget(Tau=Tau2, 
                                       X_PMT = X_PMT, 
                                       X_CBT = X_CBT2,
                                       X_program = X_program,
                                       X_elite = "minority",
                                       Y_micro = Y_micro, #needs to be a matrix, not vector
                                       iter_keep = iter_keep,
                                       iter_burn = iter_burn,
                                       print_opt = print_opt,
                                       groups = groups,
                                       prior_prob_rank = prior_prob_rank,
                                       initial.list = initial_list_noelite)

    
apply(Hybridtemp_noelite$omega_rank, 2, mean)  
    
        