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
#parameters for simulation


iter_keep = 2000   ## Gibbs sampler kept iterations (post burn-in)
iter_burn =2000   ## Gibbs sampler burn-in iterations 
print_opt = 1000  ## print a message every print.opt steps


full_data <- read.csv("Alatas Analysis/alatas.csv") %>%
  dplyr::select(-c("hhsize_ae")) %>% arrange(village, province, district, subdistrict)%>% 
  mutate(community_id = as.numeric(factor(interaction(village, province, district, subdistrict))))%>%
  group_by(village, province, district, subdistrict)%>%
  mutate(prop_rank = rank,
         poverty_rate = mean(treated)) %>%
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
           "formal", "informal")
m3 <- c(m_num, m_bin)

#50% of the full data is surveyed for PMT. get both X and y=consumption
#set.seed(572319852)
PMT_idx <-which(full_data$pmt == 1)#which(full_data$community_id %in% sample(unique(full_data$community_id), replace=FALSE,  length(unique(full_data$community_id))*.5))
#Note: the hh index for program is everything else, e.g. , full_data[-PMT_idx,]
#a subset of the program data is CBT

full_data <- full_data[-PMT_idx,] %>%mutate_at(m_num, function(x){(x - mean(x))/(2*sd(x))}) #%>%
#mutate(hhage2 = hhage^2,
#       hhsize2 = hhsize^2)



  i <- 0
  r <- list()
  c <- list()
  for(rep in c(1:10)){
    print(paste("***********Rep ", rep, "**************"))
    i = i + 1
    
    #SAMPLE PROGRAM DATA - FOR OUT OF SAMPLE TESTING
    whats_left <- unique(full_data$community_id)
    samps <- data.frame(community_id = whats_left,
                        samp = sample(c(1:3), prob = c(.76,.12, .12), size = length(whats_left),replace=TRUE))
    
    CBT1_data <- full_data %>%subset(community_id %in% samps$community_id[samps$samp == 2])
    CBT2_data <- full_data %>%subset(community_id %in% samps$community_id[samps$samp == 3])
    Program_data <- full_data %>%subset(community_id %in% samps$community_id[samps$samp == 1])
    
    X_CBT1 <-     cbind(1,CBT1_data[,m3]) %>%as.matrix()
    X_CBT2 <-     cbind(1,CBT2_data[,m3]) %>%as.matrix()
    X_program <- cbind(1,Program_data[,m3]) %>%as.matrix()
    X_program_noelite <- X_program
    X_program_noelite[,"connected"] <- 0
    
    R1 = CBT1_data %>% group_by(village, province, district, subdistrict) %>% summarise(n = length(cow))%>%ungroup() %>%nrow
    R2 = CBT2_data %>% group_by(village, province, district, subdistrict) %>% summarise(n = length(cow))%>%ungroup() %>%nrow
    M = 1  ## just consumption
    
    #starting values for random effects
    
    which_noelite <- which(colnames(X_CBT1) == "connected") #NOTE THIS INDEX INCLUDES THE FIRST POSITION OF INTERCEPT
    
    temp_data <- data.frame(rank = CBT1_data$rank,
                            X_CBT1[,-1]) %>%
      mutate(rank = (rank - mean(rank))/sd(rank))
    form <- formula(paste0("rank~", paste0(colnames(X_CBT1), collapse = "+")))
    
    CBT_beta_start <-coef(lm(form, data = temp_data))%>%as.vector()
    CBT_beta_start[is.na(CBT_beta_start)] <- 0
    
    
    initial_list_noelite <- list(
      beta_rank = c(0,CBT_beta_start[-1]))
    
    
    #create rank matrix: one column per 'ranker' (community)
    Tau1 <- array(NA, dim = c(nrow(CBT1_data), R1))
    j = 0
    for ( idx in unique(CBT1_data$community_id)){ #loop over columns
      j = j + 1
      Tau1[CBT1_data$community_id == idx,j] <- CBT1_data$rank[CBT1_data$community_id == idx]
    }
    Tau2 <- array(NA, dim = c(nrow(CBT2_data), R2))
    j = 0
    for ( idx in unique(CBT2_data$community_id)){ #loop over columns
      j = j + 1
      Tau2[CBT2_data$community_id == idx,j] <- CBT2_data$rank[CBT2_data$community_id == idx]
    }
    
    

    #PERIOD 1
    CB1_results <- CBTarget(Tau=Tau1, 
                               X_CBT = X_CBT1,
                               X_program = X_program,
                               X_elite = "connected",
                               iter_keep =iter_keep,
                               iter_burn = iter_burn,
                               print_opt = print_opt,
                               initial.list = initial_list_noelite)
    
    CB_beta_rank_mean <- apply(CB1_results$beta_rank, 2, mean)
    CB_beta_rank_vcov <- cov(CB1_results$beta_rank)
    CB_omega_rank_probs <- c(mean(CB1_results$omega_rank[,1]==0.5),mean(CB1_results$omega_rank[,1]==1),mean(CB1_results$omega_rank[,1]==2))
    CB_omega_rank_probs <- apply(cbind(c(1/3,1/3,1/3), CB_omega_rank_probs), 1, mean)
    
    
    #PERIOD 2 - NO DYNAMIC UPDATING
    CB2_results <- CBTarget(Tau=Tau2, 
                            X_CBT = X_CBT2,
                            X_program = X_program,
                            X_elite = "connected",
                            iter_keep =iter_keep,
                            iter_burn = iter_burn,
                            print_opt = print_opt,
                            initial.list = initial_list_noelite) #USING DEFAULTS FOR delta prior means, vars, etc.
    
    
    #PERIOD 2 - DYNAMIC UPDATING
    CB2_results_DU <- CBTarget(Tau=Tau2, 
                            X_CBT = X_CBT2,
                            X_program = X_program,
                            X_elite = "connected",
                            iter_keep =iter_keep,
                            iter_burn = iter_burn,
                            print_opt = print_opt,
                            initial.list = initial_list_noelite,
                            delta_prior_mean = CB_beta_rank_mean[-1],
                            delta_prior_var = 2*mean(diag(CB_beta_rank_vcov)),
                            prior_prob_rank = CB_omega_rank_probs)
    

    
    c[[i]] <- data.frame(parameter = m3,
                         rep = rep,
                         CB_DU = apply(CB2_results_DU$beta_rank, 2, mean)[-1],
                         CB = apply(CB2_results$beta_rank, 2, mean)[-1],
                         P1 = CB_beta_rank_mean[-1])
    c[[i]] %>%
      melt(id.vars = c("parameter", "rep")) %>%
    ggplot() +
      geom_line(aes(x = parameter, y = value, colour = variable, group = variable)) +
      coord_flip()
    
  
    
    #HYBRID-BASED PREDICTION - WITH DYNAMIC UPDATING
    Program_data$DU_prediction <-apply(CB2_results_DU$mu_noelite, 2, mean)
    
    #HYBRID-BASED PREDICTION -  WEAKLY INFORMATIVE PRIORS (NO DU)
    Program_data$prediction <- apply(CB2_results$mu_noelite, 2, mean)
Program_data$p1_prediction <- apply(CB1_results$mu_noelite, 2, mean)
    

    Program_data <- Program_data%>%group_by(village, province, district, subdistrict, poverty_rate) %>%
      mutate(DU_rank =rank(DU_prediction)/length(village),
             stand_rank =rank(prediction)/length(village),
             p1_rank = rank(p1_prediction)/length(village),
             cbt_rank = rank/length(village)) %>%
      mutate(DU_inclusion = DU_rank <= poverty_rate,
             inclusion = stand_rank <= poverty_rate,
             p1_inclusion = p1_rank <= poverty_rate,
             cbt_inclusion = cbt_rank <= poverty_rate) %>%ungroup() %>%
      mutate_at(vars(matches("inclusion")), as.factor)
    
    
    
    r[[i]] <- rbind(
      confusionMatrix(Program_data$DU_inclusion,  Program_data$cbt_inclusion,positive = "TRUE")$byClass,
      confusionMatrix(Program_data$inclusion,     Program_data$cbt_inclusion,positive = "TRUE")$byClass,
      confusionMatrix(Program_data$p1_inclusion,  Program_data$cbt_inclusion,positive = "TRUE")$byClass) %>%as.data.frame%>%
      mutate(Method = c("DU", "Non-DU", "P1"),
             TD = Sensitivity - (1-Specificity),
             rep = rep) 
  }
  

all_results <-  do.call(rbind, r)
all_results %>%  mutate(IER = 1-Precision,
                                     EER = 1-Sensitivity) %>%
  dplyr::select(Method, IER, rep) %>%
  ggplot() + geom_boxplot(aes(x = Method, y = IER))


all_coef <-  do.call(rbind, c)

write.csv(all_results, "Alatas Analysis/Dynamic Updating/all_results.csv")

write.csv(all_coef, "Alatas Analysis/Dynamic Updating/all_coef.csv")
