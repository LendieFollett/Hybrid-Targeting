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
           "formal", "informal")# lrf removed informal for now
m3 <- c(m_num, m_bin)

#50% of the full data is surveyed for PMT. get both X and y=consumption
#set.seed(572319852)
PMT_idx <-which(full_data$pmt == 1)#which(full_data$community_id %in% sample(unique(full_data$community_id), replace=FALSE,  length(unique(full_data$community_id))*.5))
#Note: the hh index for program is everything else, e.g. , full_data[-PMT_idx,]
#a subset of the program data is CBT

full_data <- full_data %>%mutate_at(m_num, function(x){(x - mean(x))/(2*sd(x))}) #%>%
#mutate(hhage2 = hhage^2,
#       hhsize2 = hhsize^2)


#parallelized across CBT proportions via mcapply
r <- list()
CBT_ncomm_list <- c(10, 50, 100, 200) 
nrep <- 1000
i=0
for (CBT_ncomm in CBT_ncomm_list){
  print(paste("***********Number of communities =  ", CBT_ncomm, "**************"))
  for(rep in c(1:nrep)){
    i = i + 1
    
    
    #RANDOM SAMPLES OF CBT DATA, PMT DATA, AND PROGRAM (testing) DATA
    whats_left <- unique(full_data$community_id[-PMT_idx]) #communities not in PMT
    samps <- data.frame(community_id = whats_left,
                        samp =     rep(c("CBT2", "Program", "NA"), 
                                       c(CBT_ncomm, 
                                         200, 
                                         length(whats_left) -CBT_ncomm- 200))[sample.int(length(whats_left))])
    

    
    CBT2_data <- full_data %>%subset(community_id %in% samps$community_id[samps$samp == "CBT2"])
    Program_data <- full_data %>%subset(community_id %in% samps$community_id[samps$samp == "Program"])    
    PMT_data <- full_data[PMT_idx,] 
    
    
    while(any(apply(CBT2_data[,m3], 2, var) == 0)){ #have to do to deal with complete separation and ML estimation of logistic regression (#shouldadonebayes)
      whats_left <- unique(full_data$community_id[-PMT_idx]) #communities not in PMT
      samps <- data.frame(community_id = whats_left,
                          samp =     rep(c("CBT2", "Program", "NA"), 
                                         c(CBT_ncomm, 
                                           200, 
                                           length(whats_left) -CBT_ncomm- 200))[sample.int(length(whats_left))])
      
      CBT2_data <- full_data %>%subset(community_id %in% samps$community_id[samps$samp == "CBT2"])
      Program_data <- full_data %>%subset(community_id %in% samps$community_id[samps$samp == "Program"])    
      PMT_data <- full_data[PMT_idx,] 
    }
    
    X_PMT <-     cbind(1,PMT_data[,m3]) %>%as.matrix()#cbind(1, PMT_data[,m3]%>%apply(2, function(x){(x - mean(x))/sd(x)})) 
    X_CBT2 <-     cbind(1,CBT2_data[,m3]) %>%as.matrix()
    X_program <- cbind(1,Program_data[,m3]) %>%as.matrix()
    which_noelite <- which(colnames(X_CBT2) == "connected") #NOTE THIS INDEX INCLUDES THE FIRST POSITION OF INTERCEPT
    
    #PMT - no elite bias correction - NOTE: SAMPLE SIZE VARIES
    #whats_left <- unique(full_data$community_id[PMT_idx])
    #samps <- sample(whats_left, size = CBT_ncomm, replace=FALSE)
    #PMT_data_sub <- subset(PMT_data, community_id %in% samps)
    Y_micro_sub <- as.matrix(log(CBT2_data$consumption))
    Y_micro_sub <- apply(Y_micro_sub, 2, function(x){(x - mean(x))/sd(x)})
    X_PMT_sub <-     cbind(1,CBT2_data[,m3]) %>%as.matrix()
    temp_data <- data.frame(Y_micro = Y_micro_sub,
                            X_PMT_sub)
    form <- formula(paste0("Y_micro~", paste0(colnames(X_PMT)[-which_noelite], collapse = "+")))
    PMT_model <- lm(form, data = temp_data)
    PMT_beta <-coef(PMT_model)%>%as.vector()
   

    #OLS-BASED PMT PREDICTION -  WITHOUT CORRECTION
    Program_data$PMT_prediction <- predict(PMT_model, Program_data)#(X_program[,-c(1, which_noelite)]%*%PMT_beta[-1])#beta_start is the OLS estimate of beta
    
    
    Program_data <- Program_data%>%group_by(village, province, district, subdistrict, poverty_rate) %>%
      mutate(pmt_rank =rank(PMT_prediction)/length(village),
             consumption_rank = rank(consumption)/length(village),
             cbt_rank = rank/length(village)) %>%
      mutate(pmt_inclusion = pmt_rank <= poverty_rate,
             consumption_inclusion = consumption_rank<=poverty_rate,
             cbt_inclusion = cbt_rank <= poverty_rate) %>%ungroup() %>%
      mutate_at(vars(matches("inclusion")), as.factor)
    
    r[[i]] <-rbind(
      confusionMatrix(Program_data$pmt_inclusion, Program_data$consumption_inclusion,positive = "TRUE")$byClass,
      confusionMatrix(Program_data$pmt_inclusion, Program_data$cbt_inclusion,positive = "TRUE")$byClass) %>%as.data.frame%>%
      mutate(Method = c( "PMT OLS", "PMT OLS (cbt target)"),
             CBT_ncomm = CBT_ncomm,
             nonconverge =  any(is.na(PMT_beta)),
             TD = Sensitivity - (1-Specificity),
             rep = rep)
  }
  
}
pmt_results<- do.call(rbind, r) %>%
  mutate(EER = 1-Sensitivity)
write.csv(pmt_results, "Alatas Analysis/pmt_experimental_results.csv", row.names=FALSE)

#---CONSUMPTION TRUTH STANDARD-----------

#error rates
pmt_results %>% subset(Method == "PMT OLS") %>%
  group_by(CBT_ncomm) %>%
  summarise(mean_EER = mean(EER))%>%
  write.csv("Alatas Analysis/PMT_nonconverge_corrected_consumption.csv", row.names=FALSE)
#proportion of experiments discarded
pmt_results %>% subset(Method == "PMT OLS") %>%
  group_by(CBT_ncomm) %>%
  summarise(prop_nonconverge = sum(nonconverge)/length(EER))


#---COMMUNITY RANKING TRUTH STANDARD-----------
#error rates
pmt_results %>% subset(Method == "PMT OLS (cbt target)") %>%
  group_by(CBT_ncomm) %>%
  summarise(mean_EER = mean(EER)) %>%
  mutate(Method = "PMT OLS (corrected)")%>%
  write.csv("Alatas Analysis/PMT_nonconverge_corrected_cbt.csv", row.names=FALSE)
#proportion of experiments discarded
pmt_results %>% subset(Method == "PMT OLS (cbt target)") %>%
  group_by(CBT_ncomm) %>%
  summarise(prop_nonconverge = sum(nonconverge)/length(EER))


