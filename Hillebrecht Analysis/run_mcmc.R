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

CBT_ncomm_list <- c(5,10,15,25)
nrep <- 10
results <-  mclapply(CBT_ncomm_list, function(CBT_ncomm){
  i <- 0
  r <- list()
  c <- list()
  for(rep in c(1:nrep)){
    print(paste("***********Rep ", rep," of CBT proportion ", CBT_ncomm, "**************"))
    i = i + 1
    
    
    
    whats_left <- unique(full_data_left$community)
    Program_idx <- which(full_data_left$community %in% sample(whats_left, 
                                                         replace=FALSE, 
                                                         length(whats_left)*0.5)) #this will leave 26 communities for the CBT 'training'
    whats_left <- unique(full_data_left$community[-c( Program_idx)])
    #length(whats_left)
    #whats_left
    CBT_idx <- which(full_data_left$community %in% sample(whats_left, 
                                                     replace=FALSE, 
                                                     CBT_ncomm))
    
    CBT_data <- full_data_left[CBT_idx,] #this is a subset of the program data!
    PMT_data <- full_data[PMT_idx,] 
    
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
      PMT_data <- full_data[PMT_idx,] 
      
    }
    
    Program_data <- full_data_left[Program_idx,]
    
    
    X_PMT <-     cbind(1,PMT_data[,m3]) %>%as.matrix()#cbind(1, PMT_data[,m3]%>%apply(2, function(x){(x - mean(x))/sd(x)})) 
    X_CBT <-     cbind(1,CBT_data[,m3]) %>%as.matrix()
    X_program <- cbind(1,Program_data[,m3]) %>%as.matrix()
    Y_micro <- as.matrix(log(PMT_data$consumption + .1)) #ROUGH FIX - USE OTHER TRANSFORMATION
    Y_micro <- apply(Y_micro, 2, function(x){(x - mean(x))/sd(x)})
    
    R = 3*(CBT_data %>% group_by(community) %>% summarise(n = length(floor))%>%ungroup() %>%nrow)

    #starting values 
    temp_data <- data.frame(Y_micro = Y_micro,
                            X_PMT)
    form <- formula(paste0("Y_micro~", paste0(colnames(X_PMT), collapse = "+")))
    
    PMT_beta_start <-coef(lm(form, data = temp_data))%>%as.vector()
    PMT_beta_start[is.na(PMT_beta_start)] <- 0
    temp_data <- data.frame(rank = apply(CBT_data[,c("informant1", "informant2", "informant3")], 1, mean) ,
                            X_CBT) %>%
      mutate(rank = (rank - mean(rank))/sd(rank))
    
    form <- formula(paste0("rank~", paste0(colnames(X_CBT), collapse = "+")))
    
    CBT_beta_start <-coef(lm(form, data = temp_data))%>%as.vector()
    CBT_beta_start[is.na(CBT_beta_start)] <- 0
    
    mu_beta_start <- apply(cbind(PMT_beta_start, CBT_beta_start), 1, mean) %>%c()
    
    
    initial_list<- list(
      beta_rank = c(0,CBT_beta_start[-c(1)]),
      beta_micro = PMT_beta_start,
      mu_beta = mu_beta_start[-c(1)])
    
    
    #create rank matrix: one column per 'ranker' (community)
    Tau <- array(NA, dim = c(nrow(CBT_data), R))
    j = 0
    for ( idx in unique(CBT_data$community)){ #loop over columns
      for (infmt in c("informant1", "informant2", "informant3")){
        j = j + 1
        Tau[CBT_data$community == idx,j] <- pull(CBT_data[CBT_data$community == idx,], infmt)
      }
    }
    
    #Run MCMC for Bayesian Consensus Targeting
    
    Hybridtemp <- HybridTarget(Tau=Tau, 
                         X_PMT = X_PMT, 
                         X_CBT = X_CBT,
                         X_program = X_program,
                         X_elite = NULL,
                         Y_micro = Y_micro, #needs to be a matrix, not vector
                         prior_prob_rank = c(1,1,1)/3,
                         iter_keep = iter_keep,
                         iter_burn = iter_burn,
                         print_opt = print_opt,
                         initial.list = initial_list)
    
    CBtemp <- CBTarget(Tau=Tau, 
                       X_CBT = X_CBT,
                       X_program = X_program,
                       X_elite = NULL,
                       prior_prob_rank = c(1,1,1)/3,
                       iter_keep =iter_keep,
                       iter_burn = iter_burn,
                       print_opt = print_opt,
                       initial.list = initial_list)
    
    #---Save coefficients from models with/without elite connection accounted for
    
    Hybrid_mu_beta_mean<-  apply(Hybridtemp$mu_beta, 2, mean)
    Hybrid_beta_rank_mean <- apply(Hybridtemp$beta_rank, 2, mean)
    Hybrid_beta_micro_mean <- apply(Hybridtemp$beta_micro, 2, mean)

    CB_beta_rank_mean <- apply(CBtemp$beta_rank, 2, mean)
    
    c[[i]] <- data.frame(parameter = m3,
                         rep = rep,
                         CBT_ncomm = CBT_ncomm,
                         Hybrid_mu_beta_mean = Hybrid_mu_beta_mean,
                         Hybrid_beta_rank_mean = Hybrid_beta_rank_mean[-1],
                         Hybrid_beta_micro_mean=Hybrid_beta_micro_mean[-1],

                         CB_beta_rank_mean = CB_beta_rank_mean[-1],
                         
                         PMT_beta = PMT_beta_start[-1])
    
    lr <- glm(treated~ ., data = CBT_data[,c("treated", m3)], family =binomial(link = "probit"))
    
    #HYBRID-BASED PREDICTION
    Program_data$hybrid_prediction <-        apply(Hybridtemp$mu, 2, mean)
    #CBT SCORE-BASED PREDICTION
    Program_data$cbt_model_prediction <- apply(CBtemp$mu, 2, mean)
    
    #PROBIT REGRESSION CBT-BASED PREDICTION
    Program_data$CBT_LR_prediction<- -predict(lr, as.data.frame(X_program[,-1]))
    #OLS-BASED PMT PREDICTION
    Program_data$PMT_prediction <- (X_program[,-1]%*%PMT_beta_start[-1])#PMT_beta_start is the OLS estimate of beta
    
    Program_data <- Program_data%>%group_by(community) %>%
      mutate(
        hybrid_rank =rank(hybrid_prediction)/length(community),
        pmt_rank =rank(PMT_prediction)/length(community),
        cbt_model_rank = rank(cbt_model_prediction)/length(community),
        #consumption_rank = rank(consumption)/length(community),
        CBT_LR_rank = rank(CBT_LR_prediction)/length(community)
      ) %>%
      mutate(#hybrid_inclusion = hybrid_rank <= poverty_rate,
        hybrid_inclusion = hybrid_rank <= treat_rate,
        pmt_inclusion = pmt_rank <= treat_rate,
        #consumption_inclusion = consumption_rank<=treat_rate,
        cbt_model_inclusion = cbt_model_rank<=treat_rate,
        CBT_LR_inclusion = CBT_LR_rank<=treat_rate,
        cbt_inclusion = ifelse(treated == 1, TRUE, FALSE)) %>%ungroup() %>%
      mutate_at(vars(matches("inclusion")), as.factor)
    
    
    
    r[[i]] <- rbind(
      confusionMatrix(Program_data$hybrid_inclusion, Program_data$cbt_inclusion,positive = "TRUE")$byClass,
      confusionMatrix(Program_data$cbt_model_inclusion,      Program_data$cbt_inclusion,positive = "TRUE")$byClass,
      confusionMatrix(Program_data$pmt_inclusion,            Program_data$cbt_inclusion,positive = "TRUE")$byClass,
      confusionMatrix(Program_data$CBT_LR_inclusion,            Program_data$cbt_inclusion,positive = "TRUE")$byClass) %>%as.data.frame%>%
      mutate(Method = c( "Hybrid Score","CBT Score", "PMT OLS", "CBT Probit"),
             CBT_ncomm = CBT_ncomm,
             TD = Sensitivity - (1-Specificity),
             rep = rep) #%>%
      #dplyr::select(c(Method,CBT_ncomm,Sensitivity, Specificity, TD))
    
    
  }
  
  return(list(r, c))
}, mc.cores = length(CBT_ncomm_list))

all_results_1 <- unlist(results, recursive = FALSE)

all_results_2 <- list()
for( i in seq(1,4*2-1, by = 2)){
all_results_2[[i]] <- do.call("rbind", all_results_1[[i]])
}

all_results_1 <- unlist(results, recursive = FALSE)
#collect accuracies
all_results_2 <- list()
for( i in seq(1,4*2-1, by = 2)){
  all_results_2[[i]] <- do.call("rbind", all_results_1[[i]])
}
all_results <- do.call("rbind", all_results_2)


#collect coefficients 
all_coef_2 <- list()
for( i in seq(2,4*2, by = 2)){
  all_coef_2[[i]] <- do.call("rbind", all_results_1[[i]])
}
all_coef <- do.call("rbind", all_coef_2)


write.csv(all_results, "Hillebrecht Analysis/all_results.csv")

write.csv(all_coef, "Hillebrecht Analysis/all_coef.csv")
