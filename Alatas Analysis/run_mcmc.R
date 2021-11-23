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


iter_keep = 20   ## Gibbs sampler kept iterations (post burn-in)
iter_burn =20   ## Gibbs sampler burn-in iterations 
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


#parallelized across CBT proportions via mcapply
CBT_ncomm_list <- c(10, 25, 50, 200) 
 results <-  mclapply(CBT_ncomm_list, function(CBT_ncomm){
   i <- 0
   r <- list()
   c <- list()
  for(rep in c(1:10)){
    print(paste("***********Rep ", rep," of CBT proportion ", CBT_ncomm, "**************"))
    i = i + 1

#SAMPLE PROGRAM DATA - FOR OUT OF SAMPLE TESTING
whats_left <- unique(full_data$community_id[- PMT_idx])
Program_idx <- which(full_data$community_id %in% sample(whats_left, 
                                                    replace=FALSE, 
                                                    length(whats_left)*0.5))
#SAMPLE CBT DATA - FOR TRAINING
whats_left <- unique(full_data$community_id[-c(PMT_idx, Program_idx)])
CBT_idx <- which(full_data$community_id %in% sample(whats_left, 
                                                   replace=FALSE, 
                                                   CBT_ncomm))


CBT_data <- full_data[CBT_idx,] #this is a subset of the program data!
PMT_data <- full_data[PMT_idx,] #%>% subset(community == 0)

while(any(apply(CBT_data[,m3], 2, var) == 0)){ #have to do to deal with complete separation and ML estimation of logistic regression (#shouldadonebayes)
  Program_idx <- which(full_data$community_id %in% sample(unique(full_data$community_id[- PMT_idx]), 
                                                          replace=FALSE, 
                                                          length(unique(full_data$community_id[-PMT_idx]))*0.5))
  CBT_idx <-which(full_data$community_id %in% sample(unique(full_data$community_id[-c(PMT_idx, Program_idx)]), 
                                                     replace=FALSE, 
                                                     CBT_ncomm))
  
  
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

#starting values for random effects

which_noelite <- which(colnames(X_CBT) == "connected") #NOTE THIS INDEX INCLUDES THE FIRST POSITION OF INTERCEPT

temp_data <- data.frame(Y_micro = Y_micro,
                        X_PMT)
form <- formula(paste0("Y_micro~", paste0(colnames(X_PMT), collapse = "+")))

PMT_beta_start <-coef(lm(form, data = temp_data))%>%as.vector()
PMT_beta_start[is.na(PMT_beta_start)] <- 0
temp_data <- data.frame(rank = CBT_data$rank,
                        X_CBT) %>%
  mutate(rank = (rank - mean(rank))/sd(rank))
form <- formula(paste0("rank~", paste0(colnames(X_CBT), collapse = "+")))

CBT_beta_start <-coef(lm(form, data = temp_data))%>%as.vector()
CBT_beta_start[is.na(CBT_beta_start)] <- 0

mu_beta_start <- apply(cbind(PMT_beta_start, CBT_beta_start), 1, mean) %>%c()

initial_list_noelite <- list(
  beta_rank = c(0,CBT_beta_start[-1]),
  beta_micro = PMT_beta_start,
  mu_beta = mu_beta_start[-1])

initial_list<- list(
  beta_rank = c(0,CBT_beta_start[-c(1, which_noelite)]),
  beta_micro = PMT_beta_start[-which_noelite],
  mu_beta = mu_beta_start[-c(1, which_noelite)])

#create rank matrix: one column per 'ranker' (community)
Tau <- array(NA, dim = c(nrow(CBT_data), R))
j = 0
for ( idx in unique(CBT_data$community_id)){ #loop over columns
  j = j + 1
  Tau[CBT_data$community_id == idx,j] <- CBT_data$rank[CBT_data$community_id == idx]
}

#Run MCMC for Bayesian Consensus Targeting - WITH CORRECTION
Hybridtemp_noelite <- HybridTarget(Tau=Tau, 
                 X_PMT = X_PMT, 
                 X_CBT = X_CBT,
                 X_program = X_program,
                 X_elite = "connected",
                 Y_micro = Y_micro, #needs to be a matrix, not vector
                 iter_keep = iter_keep,
                 iter_burn = iter_burn,
                  print_opt = print_opt,
                 initial.list = initial_list_noelite)

#Run MCMC for Bayesian Consensus Targeting - WITHOUT CORRECTION
Hybridtemp <- HybridTarget(Tau=Tau, 
                           X_PMT = X_PMT[,-which(colnames(X_PMT) == "connected")], 
                           X_CBT = X_CBT[,-which(colnames(X_CBT) == "connected")],
                           X_program = X_program[,-which(colnames(X_program) == "connected")],
                           X_elite = NULL,
                           Y_micro = Y_micro, #needs to be a matrix, not vector
                           iter_keep = iter_keep,
                           iter_burn = iter_burn,
                           print_opt = print_opt,
                           initial.list = initial_list)

#Run MCMC for Bayesian Community Based Targeting -  WITH CORRECTION
CBtemp_noelite <- CBTarget(Tau=Tau, 
                 X_CBT = X_CBT,
                 X_program = X_program,
                 X_elite = "connected",
                 iter_keep =iter_keep,
                 iter_burn = iter_burn,
                 print_opt = print_opt,
                 initial.list = initial_list_noelite)

#Run MCMC for Bayesian Community Based Targeting -  WITHOUT CORRECTION
CBtemp <- CBTarget(Tau=Tau, 
                   X_CBT = X_CBT[,-which(colnames(X_CBT) == "connected")],
                   X_program = X_program[,-which(colnames(X_program) == "connected")],
                   X_elite =NULL,
                   iter_keep =iter_keep,
                   iter_burn = iter_burn,
                   print_opt = print_opt,
                   initial.list = initial_list)


#PMT - no elite bias correction
temp_data <- data.frame(Y_micro = Y_micro,
                        X_PMT)
form <- formula(paste0("Y_micro~", paste0(colnames(X_PMT)[-which_noelite], collapse = "+")))

PMT_beta <-coef(lm(form, data = temp_data))%>%as.vector()


#---Save coefficients from models with/without elite connection accounted for
Hybrid_mu_beta_mean_noelite <- apply(Hybridtemp_noelite$mu_beta, 2, mean)
Hybrid_beta_rank_mean_noelite <- apply(Hybridtemp_noelite$beta_rank, 2, mean)
Hybrid_beta_micro_mean_noelite <- apply(Hybridtemp_noelite$beta_micro, 2, mean)

Hybrid_mu_beta_mean<- append( apply(Hybridtemp$mu_beta, 2, mean), 0, after = which_noelite-2)
Hybrid_beta_rank_mean <- append(apply(Hybridtemp$beta_rank, 2, mean), 0, after = which_noelite-1)
Hybrid_beta_micro_mean <- append(apply(Hybridtemp$beta_micro, 2, mean), 0, after = which_noelite-1)

CB_beta_rank_mean_noelite <- apply(CBtemp_noelite$beta_rank, 2, mean)
CB_beta_rank_mean <- append(apply(CBtemp$beta_rank, 2, mean), 0, after = which_noelite-1)

c[[i]] <- data.frame(parameter = m3,
                    rep = rep,
                    CBT_ncomm=CBT_ncomm,
           Hybrid_mu_beta_mean_noelite = Hybrid_mu_beta_mean_noelite,
           Hybrid_beta_rank_mean_noelite = Hybrid_beta_rank_mean_noelite[-1],
           Hybrid_beta_micro_mean_noelite=Hybrid_beta_micro_mean_noelite[-1],
           
           Hybrid_mu_beta_mean = Hybrid_mu_beta_mean,
           Hybrid_beta_rank_mean = Hybrid_beta_rank_mean[-1],
           Hybrid_beta_micro_mean=Hybrid_beta_micro_mean[-1],
           
           CB_beta_rank_mean_noelite = CB_beta_rank_mean_noelite[-1],
           CB_beta_rank_mean = CB_beta_rank_mean[-1],
           
           PMT_beta = append(PMT_beta[-1], 0, after = which_noelite-2))



#Fit probit regression for Community Based Targeting
lr <- glm(prop_rank<=poverty_rate ~ ., 
          data = CBT_data[,c("prop_rank","poverty_rate", m3[-which(m3 == "connected")])], 
          family =binomial(link = "probit"))

#HYBRID-BASED PREDICTION - WITH CORRECTION
Program_data$hybrid_prediction_noelite <-apply(Hybridtemp_noelite$mu_noelite, 2, mean)
#HYBRID-BASED PREDICTION - WITHOUT CORRECTION
Program_data$hybrid_prediction <-apply(Hybridtemp$mu, 2, mean)

#CBT SCORE-BASED PREDICTION -  WITH CORRECTION
Program_data$cbt_model_prediction_noelite <- apply(CBtemp_noelite$mu_noelite, 2, mean)
#CBT SCORE-BASED PREDICTION -  WITHOUT CORRECTION
Program_data$cbt_model_prediction <- apply(CBtemp$mu, 2, mean)

#PROBIT CBT-BASED PREDICTION -  WITHOUT CORRECTION
Program_data$CBT_LR_prediction<- -predict(lr, as.data.frame(X_program[,-1])) #(LRF NEEDS TO CHANGE TO) logistic regression

#OLS-BASED PMT PREDICTION -  WITHOUT CORRECTION
Program_data$PMT_prediction <- (X_program[,-c(1, which_noelite)]%*%PMT_beta[-1])#beta_start is the OLS estimate of beta

Program_data <- Program_data%>%group_by(village, province, district, subdistrict, poverty_rate) %>%
  mutate(hybrid_noelite_rank =rank(hybrid_prediction_noelite)/length(village),
         hybrid_rank =rank(hybrid_prediction)/length(village),
         pmt_rank =rank(PMT_prediction)/length(village),
         cbt_model_rank = rank(cbt_model_prediction)/length(village),
         cbt_model_rank_noelite = rank(cbt_model_prediction_noelite)/length(village),
         consumption_rank = rank(consumption)/length(village),
         cbt_rank = rank/length(village),
         CBT_LR_rank = rank(CBT_LR_prediction)/length(village)) %>%
  mutate(hybrid_noelite_inclusion = hybrid_noelite_rank <= poverty_rate,
         hybrid_inclusion = hybrid_rank <= poverty_rate,
         pmt_inclusion = pmt_rank <= poverty_rate,
         consumption_inclusion = consumption_rank<=poverty_rate,
         cbt_model_inclusion = cbt_model_rank<=poverty_rate,
         cbt_model_noelite_inclusion = cbt_model_rank_noelite<=poverty_rate,
         CBT_LR_inclusion = CBT_LR_rank<=poverty_rate,
         cbt_inclusion = cbt_rank <= poverty_rate) %>%ungroup() %>%
  mutate_at(vars(matches("inclusion")), as.factor)



r[[i]] <- rbind(
confusionMatrix(Program_data$hybrid_noelite_inclusion,   Program_data$cbt_inclusion,positive = "TRUE")$byClass,
confusionMatrix(Program_data$hybrid_inclusion,           Program_data$cbt_inclusion,positive = "TRUE")$byClass,
confusionMatrix(Program_data$cbt_model_noelite_inclusion,Program_data$cbt_inclusion,positive = "TRUE")$byClass,
confusionMatrix(Program_data$cbt_model_inclusion,        Program_data$cbt_inclusion,positive = "TRUE")$byClass,
confusionMatrix(Program_data$pmt_inclusion,              Program_data$cbt_inclusion,positive = "TRUE")$byClass,
confusionMatrix(Program_data$CBT_LR_inclusion,           Program_data$cbt_inclusion,positive = "TRUE")$byClass) %>%as.data.frame%>%
  mutate(Method = c( "Hybrid Score (corrected)","Hybrid Score","CBT Score (corrected)", "CBT Score", "PMT OLS", "CBT Logit"),
         CBT_ncomm = CBT_ncomm,
         TD = Sensitivity - (1-Specificity),
         rep = rep) 
  }

   return(list(r=r, c=c))
  }, mc.cores = length(CBT_ncomm_list))

 
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
 

write.csv(all_results, "Alatas Analysis/all_results.csv")

write.csv(all_coef, "Alatas Analysis/all_coef.csv")

all_results %>%melt(id.var = c("Method", "CBT_ncomm")) %>%
  mutate(Method = factor(Method, levels = c("Hybrid Score (corrected)","Hybrid Score","CBT Score", "CBT Score (corrected)", "CBT Logit", "PMT OLS"),
                         labels = c("Hybrid Score (corrected)","Hybrid Score","CBT Score","CBT Score (corrected)", "CBT Probit", "PMT OLS"))) %>%
  group_by(Method, CBT_ncomm, variable) %>%
  mutate(mean = median(value ))%>%ungroup%>%
  subset(variable %in% c("TD", "Sensitivity", "Specificity") )%>%
  #subset(Method != "CBT Probit" & Method != "PMT OLS")%>%
  ggplot() +#geom_boxplot(aes(x = Method, y = value,linetype = Method, group = interaction(Method, CBT_prop))) + 
  geom_boxplot(aes(x = Method, y = value, colour = Method, group = interaction(CBT_ncomm, Method))) + 
  stat_summary(aes(x = Method, y = value, colour = Method, group = interaction(CBT_ncomm, Method)),
               fun=mean, geom="point", color="black")+
  #geom_line(aes(x = CBT_prop, y = mean, group = interaction(Method), linetype = Method, colour = Method)) + 
  facet_grid(variable~CBT_ncomm, scales = "free")+ theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = .9))  +
  scale_colour_brewer(type = "qual", palette = "Dark2")# +
  #scale_y_log10()
ggsave("Alatas Analysis/all_results.pdf")

all_results %>%
  mutate(rep = rep(rep(1:10, times = rep(5, 10)),4))%>%
  melt(id.var = c("Method", "CBT_prop", "rep")) %>%
  mutate(Method = factor(Method, levels = c("Hybrid Score (corrected)","Hybrid Score","CBT Score", "CBT Logit", "PMT OLS"),
                         labels = c("Hybrid Score (corrected)","Hybrid Score","CBT Score", "CBT Probit", "PMT OLS"))) %>%
  group_by(Method, CBT_prop, variable) %>%
  mutate(mean = median(value ))%>%ungroup%>%
  subset(Method != "CBT Probit" & Method != "PMT OLS")%>%
  ggplot() +#geom_boxplot(aes(x = Method, y = value,linetype = Method, group = interaction(Method, CBT_prop))) + 
  geom_line(aes(x = Method, y = value,  group = interaction(CBT_prop, rep))) + 
  stat_summary(aes(x = Method, y = value, colour = Method, group = interaction(CBT_prop, Method)),
               fun=mean, geom="point", color="black")+
  #geom_line(aes(x = CBT_prop, y = mean, group = interaction(Method), linetype = Method, colour = Method)) + 
  facet_grid(variable~CBT_prop, scales = "free")+ theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = .9))  +
  scale_colour_brewer(type = "qual", palette = "Dark2")


