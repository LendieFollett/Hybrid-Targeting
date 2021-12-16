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
ihs_trans <- function(x){log(x + sqrt(x^2 + 1))}

iter_keep = 4000   ## Gibbs sampler kept iterations (post burn-in)
iter_burn = 1   ## Gibbs sampler burn-in iterations 
print_opt = 100  ## print a message every print.opt steps


full_data <- read.csv("Data/Burkina Faso/Cleaning/hillebrecht.csv") %>%
  group_by(community, year)%>%
  mutate(informant1 = ifelse(is.na(informant1), NA, floor(rank(-informant1))),
         informant2 = ifelse(is.na(informant2), NA, floor(rank(-informant2))),
         informant3 = ifelse(is.na(informant3), NA, floor(rank(-informant3))),
         treat_rate = sum(treated)/length(treated)) %>% ungroup%>%  arrange(community)
#x variables to include in model
m_num <- c("rooms", "hhsize","age1660","age60")

m_bin <- colnames(full_data)[which(colnames(full_data)=="floor"):which(colnames(full_data)=="pig")]
m_bin <- m_bin[!m_bin %in% m_num]
m3 <- c(m_num, m_bin)

full_data <- full_data %>%mutate_at(m_num, function(x){(x - mean(x))/(2*sd(x))})

PMT_idx <-which(full_data$year == 2008) #training data is all of 2008 data
full_data_left <- full_data[-PMT_idx,]

    
    #RANDOM SAMPLES OF CBT DATA, PMT DATA, AND PROGRAM (testing) DATA

    CBT2_data <- Program_data <- full_data_left

    X_CBT2 <-     cbind(1,CBT2_data[,m3]) %>%as.matrix()
    X_program <- cbind(1,Program_data[,m3]) %>%as.matrix()

    R2 = 3*(CBT2_data %>% group_by(community) %>% summarise(n = length(floor))%>%ungroup() %>%nrow)
    
    which_noelite <- which(colnames(X_CBT2) == "minority") #NOTE THIS INDEX INCLUDES THE FIRST POSITION OF INTERCEPT
    
    temp_data <- data.frame(rank = apply(CBT2_data[,c("informant1", "informant2", "informant3")], 1, mean) ,
                            X_CBT2) %>%
      mutate(rank = (rank - mean(rank))/sd(rank))
    
    form <- formula(paste0("rank~", paste0(colnames(X_CBT2), collapse = "+")))
    
    CBT_beta_start <-coef(lm(form, data = temp_data))%>%as.vector()
    CBT_beta_start[is.na(CBT_beta_start)] <- 0
    
    initial_list_noelite <- list(
      beta_rank = c(0,CBT_beta_start[-1]))
    
    #create rank matrix: one column per 'ranker' (community)
    
    Tau2 <- array(NA, dim = c(nrow(CBT2_data), R2))
    j = 0
    for ( idx in unique(CBT2_data$community)){ #loop over columns
      for (infmt in c("informant1", "informant2", "informant3")){
        j = j + 1
        temp = pull(CBT2_data[CBT2_data$community == idx,], infmt)  
        if (infmt == "informant3"){ #SCRAMBLED RANKS FOR THIRD RANKER
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
    
    #Run MCMC for Bayesian Community Based Targeting -  WITH CORRECTION
    CBtemp_noelite <- CBTarget(Tau=Tau2, 
                               X_CBT = X_CBT2,
                               X_program = X_program,
                               X_elite = "minority",
                               iter_keep =iter_keep,
                               iter_burn = iter_burn,
                               print_opt = print_opt,
                               initial.list = initial_list_noelite)

    
apply(CBtemp_noelite$omega_rank, 2, mean)  
    
        