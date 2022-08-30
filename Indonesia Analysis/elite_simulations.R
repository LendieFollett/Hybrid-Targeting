rm(list = ls())
library(truncnorm)
library(mvtnorm)
library(LaplacesDemon)
library(MASS)
library(dplyr)
library(ggplot2)
library(Rcpp)
library(reshape2)
library(caret)
library(parallel)
library(RcppTN)
detectCores(logical=FALSE)


sourceCpp("Hybrid Targeting/functions.cpp")
source("Hybrid Targeting/functions.R")
source("Hybrid Targeting/HybridTarget.R")
source("Hybrid Targeting/CBTarget.R")
detectCores(logical=FALSE)
#parameters for simulation
iter_keep = 2000   ## Gibbs sampler kept iterations (post burn-in)
iter_burn =2000   ## Gibbs sampler burn-in iterations 
print_opt = 500  ## print a message every print.opt steps


full_data <- read.csv("Data/Indonesia/Cleaning/alatas.csv") %>%
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
           "formal", "informal")# lrf removed informal for now
m3 <- c(m_num, m_bin)


#parallelized across CBT proportions via mcapply
CBT_ncomm_list <- c(10, 25, 75, 100) 
nrep <- 3

do_ploops <- function(full_data){
results <-  mclapply(CBT_ncomm_list, function(CBT_ncomm){
  i <- 0
  r <- list()
  for(rep in c(1:nrep)){
    print(paste("***********Rep ", rep," of CBT proportion ", CBT_ncomm, "**************"))
    i = i + 1
    
    
    #RANDOM SAMPLES OF CBT DATA, PMT DATA, AND PROGRAM (testing) DATA
    whats_left <- unique(full_data$community_id) #communities not in PMT
    samps <- data.frame(community_id = whats_left,
                        samp =     rep(c("CBT2", "Program", "NA"), 
                                       c(CBT_ncomm, 
                                         100, 
                                         length(whats_left) -CBT_ncomm- 100))[sample.int(length(whats_left))])
    
    
    CBT2_data <- full_data %>%subset(community_id %in% samps$community_id[samps$samp == "CBT2"])
    Program_data <- full_data %>%subset(community_id %in% samps$community_id[samps$samp == "Program"])    
    
    while(any(apply(CBT2_data[,m3], 2, var) == 0)){ #have to do to deal with complete separation and ML estimation of logistic regression (#shouldadonebayes)
      whats_left <- unique(full_data$community_id) #communities not in PMT
      samps <- data.frame(community_id = whats_left,
                          samp =     rep(c("CBT2", "Program", "NA"), 
                                         c(CBT_ncomm, 
                                           100, 
                                           length(whats_left) -CBT_ncomm- 100))[sample.int(length(whats_left))])
      
      
      CBT2_data <- full_data %>%subset(community_id %in% samps$community_id[samps$samp == "CBT2"])
      Program_data <- full_data %>%subset(community_id %in% samps$community_id[samps$samp == "Program"])    
    }

    X_CBT2 <-     cbind(1,CBT2_data[,m3]) %>%as.matrix()
    X_program <- cbind(1,Program_data[,m3]) %>%as.matrix()

    R2 = CBT2_data %>% group_by(village, province, district, subdistrict) %>% summarise(n = length(cow))%>%ungroup() %>%nrow

    #Create Rank matrix in format required for CBTarget()
    Tau2 <- array(NA, dim = c(nrow(CBT2_data), R2))
    j = 0
    for ( idx in unique(CBT2_data$community_id)){ #loop over columns
      j = j + 1
      Tau2[CBT2_data$community_id == idx,j] <- CBT2_data$rank[CBT2_data$community_id == idx]
    }
    
    #Run MCMC for Bayesian Community Based Targeting -  WITHOUT CORRECTION
    CBtemp <- CBTarget(Tau=Tau2, 
                       X_CBT = X_CBT2[,-which(colnames(X_CBT2) == "connected")],
                       X_program = X_program[,-which(colnames(X_program) == "connected")],
                       X_elite =NULL,
                       iter_keep =iter_keep,
                       iter_burn = iter_burn,
                       print_opt = print_opt)

  
    #CBT SCORE-BASED PREDICTION -  WITHOUT CORRECTION
    Program_data$cbt_model_prediction <- apply(CBtemp$mu, 2, mean)
    
    Program_data <- Program_data%>%group_by(village, province, district, subdistrict, poverty_rate) %>%
      mutate(cbt_model_rank = rank(cbt_model_prediction)/length(village),
             cbt_rank = rank/length(village)) %>%
      ungroup()
    
    #ith list element contains one row per household, all rank predictions, rep id, and number CBT communities
    r[[i]] <- Program_data %>% dplyr::select(c(hhid, village, province, district, subdistrict,poverty_rate, 
                                               cbt_model_rank, cbt_rank)) %>%
      mutate(rep = rep,
             CBT_ncomm = CBT_ncomm)
    
  }
  
  return(r)
}, mc.cores = length(CBT_ncomm_list))
}

PMT_idx <-which(full_data$pmt == 1)
full_data <- full_data[-PMT_idx,]

full_data_elite1 <- full_data %>% subset(elite == 1)
full_data_elite0 <- full_data %>% subset(elite == 0)
full_data_day1 <- full_data %>% subset(daymeeting == 1)
full_data_day0 <- full_data %>% subset(daymeeting == 0)

r_elite1 <- do_ploops(full_data_elite1)
r_elite0 <- do_ploops(full_data_elite0)
r_day1 <- do_ploops(full_data_day1)
r_day0 <- do_ploops(full_data_day0)


all_results_1 <- unlist(results, recursive = FALSE)
#collect accuracies
all_results_2 <- list()
for( i in seq(1,4*2-1, by = 2)){
  all_results_2[[i]] <- do.call("rbind", all_results_1[[i]])
}
all_results <- do.call("rbind", all_results_2)




