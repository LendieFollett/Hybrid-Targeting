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



ihs_trans <- function(x){log(x + sqrt(x^2 + 1))}

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

#50% of the full data is surveyed for PMT. get both X and y=consumption
#set.seed(572319852)
full_data <- full_data %>%mutate_at(m_num, function(x){(x - mean(x))/(2*sd(x))})

PMT_idx <-which(full_data$year == 2008) #training data is all of 2008 data
full_data_left <- full_data[-PMT_idx,]
PMT_data <- full_data[PMT_idx,] #%>% subset(community == 0)



CBT_ncomm_list <- c(5,10,15,25)
nrep <- 1000
r <- list()
i = 0
for(CBT_ncomm in CBT_ncomm_list){
  print(paste("***********Number of communities =  ", CBT_ncomm, "**************"))
  for(rep in c(1:nrep)){
    i = i + 1
    
    
    #RANDOM SAMPLES OF CBT DATA, PMT DATA, AND PROGRAM (testing) DATA
    whats_left <- unique(full_data_left$community) #communities not in PMT
    samps <- data.frame(community = whats_left,
                        samp =     rep(c("CBT2", "Program", "NA"), 
                                       c(CBT_ncomm, 
                                         25, max(51-CBT_ncomm - 25, 0)))[sample.int(length(whats_left))])
    
    
    CBT2_data <- full_data_left %>%subset(community %in% samps$community[samps$samp == "CBT2"])
    Program_data <- full_data_left %>%subset(community %in% samps$community[samps$samp == "Program"])    
    

    
    while(any(apply(CBT2_data[,m3], 2, var) == 0)){ #have to do to deal with complete separation and ML estimation of logistic regression (#shouldadonebayes)
      whats_left <- unique(full_data_left$community) #communities not in PMT
      samps <- data.frame(community = whats_left,
                          samp =     rep(c("CBT2", "Program", "NA"), 
                                         c(CBT_ncomm, 
                                           25, max(51-CBT_ncomm - 25, 0)))[sample.int(length(whats_left))])
      CBT2_data <- full_data_left %>%subset(community %in% samps$community[samps$samp == "CBT2"])
      Program_data <- full_data_left %>%subset(community %in% samps$community[samps$samp == "Program"])    
      
    }
    
    X_PMT <-     cbind(1,PMT_data[,m3]) %>%as.matrix()#cbind(1, PMT_data[,m3]%>%apply(2, function(x){(x - mean(x))/sd(x)})) 
    
    X_CBT2 <-     cbind(1,CBT2_data[,m3]) %>%as.matrix()
    X_program <- cbind(1,Program_data[,m3]) %>%as.matrix()
    Y_micro <- as.matrix(ihs_trans(PMT_data$consumption))
    Y_micro <- apply(Y_micro, 2, function(x){(x - mean(x))/sd(x)})
    
  
    which_noelite <- which(colnames(X_CBT2) == "minority") #NOTE THIS INDEX INCLUDES THE FIRST POSITION OF INTERCEPT

    
    
    #PMT - no elite bias correction
    Y_micro_sub <- as.matrix(ihs_trans(CBT2_data$consumption))
    Y_micro_sub <- apply(Y_micro_sub, 2, function(x){(x - mean(x))/sd(x)})
    X_PMT_sub <-     cbind(1,CBT2_data[,m3]) %>%as.matrix()
    temp_data <- data.frame(Y_micro = Y_micro_sub,
                            X_PMT_sub)
    form <- formula(paste0("Y_micro~", paste0(colnames(X_PMT_sub)[-which_noelite], collapse = "+")))
    PMT_model <- lm(form, data = temp_data)
    PMT_beta <-coef(PMT_model)%>%as.vector()
    #OLS-BASED PMT PREDICTION -  WITHOUT CORRECTION
    Program_data$PMT_prediction <- predict(PMT_model, Program_data)#(X_program[,-c(1, which_noelite)]%*%PMT_beta[-1])#beta_start is the OLS estimate of beta
    
    Program_data <- Program_data%>%group_by(community) %>%
      mutate(pmt_rank =rank(PMT_prediction)/length(consumption),
             consumption_rank = rank(consumption)/length(consumption)) %>% #consumption = poverty standard
      mutate(pmt_inclusion = pmt_rank <= treat_rate,
             consumption_inclusion = consumption_rank<=treat_rate,
             cbt_inclusion = ifelse(treated == 1, TRUE, FALSE)) %>%ungroup() %>%
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
write.csv(pmt_results, "Hillebrecht Analysis/pmt_experimental_results.csv", row.names=FALSE)

#pmt_results <- read.csv("Hillebrecht Analysis/pmt_experimental_results.csv")

#---CONSUMPTION TRUTH STANDARD-----------

#error rates
pmt_results %>% subset(Method == "PMT OLS") %>%
  group_by(CBT_ncomm) %>%
  summarise(mean_EER = mean(EER))%>%
  write.csv("Hillebrecht Analysis/PMT_nonconverge_corrected_consumption.csv", row.names=FALSE)
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
  write.csv("Hillebrecht Analysis/PMT_nonconverge_corrected_cbt.csv", row.names=FALSE)
#proportion of experiments discarded
pmt_results %>% subset(Method == "PMT OLS (cbt target)") %>%
  group_by(CBT_ncomm) %>%
  summarise(prop_nonconverge = sum(nonconverge)/length(EER))


