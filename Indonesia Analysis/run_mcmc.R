
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


PMT_idx <-which(full_data$pmt == 1)

full_data <- full_data %>%mutate_at(m_num, function(x){(x - mean(x))/(2*sd(x))})

#parallelized across CBT proportions via mcapply
CBT_ncomm_list <- c(10, 50, 100, 200) 
nrep <- 30

 results <-  mclapply(CBT_ncomm_list, function(CBT_ncomm){
   i <- 0
   r <- list()
   c <- list()
  for(rep in c(1:nrep)){
    print(paste("***********Rep ", rep," of CBT proportion ", CBT_ncomm, "**************"))
    i = i + 1

    
#RANDOM SAMPLES OF CBT DATA, PMT DATA, AND PROGRAM (testing) DATA
whats_left <- unique(full_data$community_id[-PMT_idx]) #communities not in PMT
samps <- data.frame(community_id = whats_left,
                    samp =     rep(c("CBT2", "Program", "NA"), 
                                   c(CBT_ncomm, 
                                     200, 
                                     length(whats_left) -CBT_ncomm- 200))[sample.int(length(whats_left))])

print(table(samps$samp))

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
Y_micro <- as.matrix(log(PMT_data$consumption))
Y_micro <- apply(Y_micro, 2, function(x){(x - mean(x))/sd(x)})

R2 = CBT2_data %>% group_by(village, province, district, subdistrict) %>% summarise(n = length(cow))%>%ungroup() %>%nrow

which_noelite <- which(colnames(X_CBT2) == "connected") #NOTE THIS INDEX INCLUDES THE FIRST POSITION OF INTERCEPT

#For initial values
temp_data <- data.frame(Y_micro = Y_micro,
                        X_PMT[,-1])
form <- formula(paste0("Y_micro~", paste0(colnames(temp_data[,-1]), collapse = "+")))

PMT_beta_start <-coef(lm(form, data = temp_data))%>%as.vector()
PMT_beta_start[is.na(PMT_beta_start)] <- 0
temp_data <- data.frame(rank = CBT2_data$rank,
                        X_CBT2[,-1]) %>%
  mutate(rank = (rank - mean(rank))/sd(rank))
form <- formula(paste0("rank~", paste0(colnames(X_CBT2), collapse = "+")))

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

#Create Rank matrix in format required for CBTarget()
Tau2 <- array(NA, dim = c(nrow(CBT2_data), R2))
j = 0
for ( idx in unique(CBT2_data$community_id)){ #loop over columns
  j = j + 1
  Tau2[CBT2_data$community_id == idx,j] <- CBT2_data$rank[CBT2_data$community_id == idx]
}

#Run MCMC for Bayesian Consensus Targeting - WITH CORRECTION
Hybridtemp_noelite <- HybridTarget(Tau=Tau2, 
                 X_PMT = X_PMT, 
                 X_CBT = X_CBT2,
                 X_program = X_program,
                 X_elite = "connected",
                 Y_micro = Y_micro, #needs to be a matrix, not vector
                 iter_keep = iter_keep,
                 iter_burn = iter_burn,
                  print_opt = print_opt,
                 initial.list = initial_list_noelite)

#Run MCMC for Bayesian Consensus Targeting - WITHOUT CORRECTION
Hybridtemp <- HybridTarget(Tau=Tau2, 
                           X_PMT = X_PMT[,-which(colnames(X_PMT) == "connected")], 
                           X_CBT = X_CBT2[,-which(colnames(X_CBT2) == "connected")],
                           X_program = X_program[,-which(colnames(X_program) == "connected")],
                           X_elite = NULL,
                           Y_micro = Y_micro, #needs to be a matrix, not vector
                           iter_keep = iter_keep,
                           iter_burn = iter_burn,
                           print_opt = print_opt,
                           initial.list = initial_list)

#Run MCMC for Bayesian Community Based Targeting -  WITH CORRECTION
CBtemp_noelite <- CBTarget(Tau=Tau2, 
                 X_CBT = X_CBT2,
                 X_program = X_program,
                 X_elite = "connected",
                 iter_keep =iter_keep,
                 iter_burn = iter_burn,
                 print_opt = print_opt,
                 initial.list = initial_list_noelite)



#Run MCMC for Bayesian Community Based Targeting -  WITHOUT CORRECTION
CBtemp <- CBTarget(Tau=Tau2, 
                   X_CBT = X_CBT2[,-which(colnames(X_CBT2) == "connected")],
                   X_program = X_program[,-which(colnames(X_program) == "connected")],
                   X_elite =NULL,
                   iter_keep =iter_keep,
                   iter_burn = iter_burn,
                   print_opt = print_opt,
                   initial.list = initial_list)

# -----PMT-------------------------------------
Y_micro_sub <- as.matrix(log(CBT2_data$consumption))
Y_micro_sub <- apply(Y_micro_sub, 2, function(x){(x - mean(x))/sd(x)})
X_PMT_sub <-     cbind(1,CBT2_data[,m3]) %>%as.matrix()
temp_data <- data.frame(Y_micro = Y_micro_sub,
                        X_PMT_sub)
form <- formula(paste0("Y_micro~", paste0(colnames(X_PMT)[-which_noelite], collapse = "+")))

PMT_model <- lm(form, data = temp_data)
PMT_beta <-coef(PMT_model)%>%as.vector()

# -----CBT probit-------------------------------------

lr <- glm(prop_rank<=poverty_rate ~ ., 
          data = CBT2_data[,c("prop_rank","poverty_rate", m3[-which(m3 == "connected")])], 
          family =binomial(link = "probit"))

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
Program_data$PMT_prediction <- predict(PMT_model, Program_data)


Program_data <- Program_data%>%group_by(village, province, district, subdistrict, poverty_rate) %>%
  mutate(hybrid_noelite_rank =rank(hybrid_prediction_noelite)/length(village),
         hybrid_rank =rank(hybrid_prediction)/length(village),
         pmt_rank =rank(PMT_prediction)/length(village),
         cbt_model_rank = rank(cbt_model_prediction)/length(village),
         cbt_model_rank_noelite = rank(cbt_model_prediction_noelite)/length(village),
         consumption_rank = rank(consumption)/length(village),
         cbt_rank = rank/length(village),
         CBT_LR_rank = rank(CBT_LR_prediction)/length(village)) %>%
  ungroup()

#ith list element contains one row per household, all rank predictions, rep id, and number CBT communities
r[[i]] <- Program_data %>% dplyr::select(c(hhid, village, province, district, subdistrict,poverty_rate, 
                                 hybrid_noelite_rank, hybrid_rank,
                                 cbt_model_rank, cbt_model_rank_noelite,
                                 consumption_rank, cbt_rank,CBT_LR_rank,
                                 pmt_rank)) %>%
                                    mutate(rep = rep,
                                          CBT_ncomm = CBT_ncomm)

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
 
#all_results contains all households from each replication
write.csv(all_results, "Indonesia Analysis/all_results.csv")

write.csv(all_coef, "Indonesia Analysis/all_coef.csv")




#-------COEFFICIENT ESTIMATES USING FULL DATA-------------------

#RANDOM SAMPLES OF CBT DATA, PMT DATA, AND PROGRAM (testing) DATA
whats_left <- unique(full_data$community_id[-PMT_idx]) #communities not in PMT

CBT_data <- full_data %>%subset(community_id %in% whats_left)
  
PMT_data <- CBT_data#full_data[PMT_idx,] #%>% subset(community == 0)

Program_data <- full_data[1:10,]  


X_PMT <-     cbind(1,PMT_data[,m3]) %>%as.matrix()#cbind(1, PMT_data[,m3]%>%apply(2, function(x){(x - mean(x))/sd(x)})) 
X_CBT <-     cbind(1,CBT_data[,m3]) %>%as.matrix()
X_program <- cbind(1,Program_data[,m3]) %>%as.matrix()

R = CBT_data %>% group_by(village, province, district, subdistrict) %>% summarise(n = length(cow))%>%ungroup() %>%nrow

Tau <- array(NA, dim = c(nrow(CBT_data), R))
j = 0
for ( idx in unique(CBT_data$community_id)){ #loop over columns
  j = j + 1
  Tau[CBT_data$community_id == idx,j] <- CBT_data$rank[CBT_data$community_id == idx]
}


#Run MCMC for Bayesian Community Based Targeting -  WITHOUT CORRECTION
CBtemp <- CBTarget(Tau=Tau, 
                   X_CBT = X_CBT[,-which(colnames(X_CBT) == "connected")],
                   X_program = X_program[,-which(colnames(X_program) == "connected")],
                   X_elite =NULL,
                   iter_keep =iter_keep,
                   iter_burn = iter_burn,
                   print_opt = 100)
#Run MCMC for Bayesian Community Based Targeting -  WITH CORRECTION
CBtemp_noelite <- CBTarget(Tau=Tau, 
                           X_CBT = X_CBT,
                           X_program = X_program,
                           X_elite = "connected",
                           iter_keep =iter_keep,
                           iter_burn = iter_burn,
                           print_opt = print_opt)

which_noelite <- which(colnames(X_CBT) == "connected") #NOTE THIS INDEX INCLUDES THE FIRST POSITION OF INTERCEPT

Y_micro <- as.matrix(log(CBT_data$consumption))
Y_micro <- apply(Y_micro, 2, function(x){(x - mean(x))/sd(x)})
X_PMT_sub <-     cbind(1,CBT_data[,m3]) %>%as.matrix()
temp_data <- data.frame(Y_micro = Y_micro,
                        X_CBT)
form <- formula(paste0("Y_micro~", paste0(colnames(X_CBT)[-which_noelite], collapse = "+")))

PMT_beta <-coef(lm(form, data = temp_data))%>%as.vector()

CB_beta_rank_mean_noelite <- apply(CBtemp_noelite$beta_rank, 2, mean)
CB_beta_rank_mean <- append(apply(CBtemp$beta_rank, 2, mean), 0, after = which_noelite-1)

coefs <- data.frame(parameter = m3,
                    #CB_beta_rank_mean_noelite = CB_beta_rank_mean_noelite[-1],
                                CB_beta_rank_mean = CB_beta_rank_mean[-1],
                                PMT_beta = append(PMT_beta[-1], 0, after = which_noelite-2)
                                )
write.csv(coefs, "Indonesia Analysis/coef_total_sample.csv", row.names = FALSE)


CB_beta_rank_pprob_noelite <- apply(CBtemp_noelite$beta_rank[,-1], 2, function(x){mean(x > 0)})
CB_beta_rank_pprob <- apply(CBtemp$beta_rank[,-1], 2, function(x){mean(x > 0)})

data.frame(parameter = m3,
           apply(CBtemp_noelite$beta_rank[,-1], 2, quantile, c(.025, .975)) %>%t(),
           mean = CB_beta_rank_mean_noelite[-1],
           pprob0 = CB_beta_rank_pprob_noelite) %>%
  write.csv( "Indonesia Analysis/CB_beta_rank_CI_noelite.csv")

data.frame(parameter = m3[-(which_noelite-1)],
           apply(CBtemp$beta_rank[,-1], 2, quantile, c(.025, .975)) %>%t(),
           mean = CB_beta_rank_mean[-c(1, which_noelite)],
           pprob0 = CB_beta_rank_pprob) %>%
  write.csv("Indonesia Analysis/CB_beta_rank_CI.csv")

apply(CBtemp_noelite$beta_rank, 2, doESS)
apply(CBtemp$beta_rank, 2, doESS)


#-------COEFFICIENT ESTIMATES USING ELITE SUBTREATMENT == 1/0-------------------
elite_status <- 1
esub_idx <- full_data %>% subset(elite == elite_status) %>% dplyr::select(community_id)
esub_idx <- as.vector(esub_idx$community_id)

#RANDOM SAMPLES OF CBT DATA, PMT DATA, AND PROGRAM (testing) DATA
whats_left <- unique(full_data$community_id[-PMT_idx]) #communities not in PMT

CBT_data <- full_data %>%subset(community_id %in% whats_left & community_id %in% esub_idx)

PMT_data <- CBT_data#full_data[PMT_idx,] #%>% subset(community == 0)

Program_data <- full_data[1:10,]  


X_PMT <-     cbind(1,PMT_data[,m3]) %>%as.matrix()#cbind(1, PMT_data[,m3]%>%apply(2, function(x){(x - mean(x))/sd(x)})) 
X_CBT <-     cbind(1,CBT_data[,m3]) %>%as.matrix()
X_program <- cbind(1,Program_data[,m3]) %>%as.matrix()

R = CBT_data %>% group_by(village, province, district, subdistrict) %>% summarise(n = length(cow))%>%ungroup() %>%nrow

Tau <- array(NA, dim = c(nrow(CBT_data), R))
j = 0
for ( idx in unique(CBT_data$community_id)){ #loop over columns
  j = j + 1
  Tau[CBT_data$community_id == idx,j] <- CBT_data$rank[CBT_data$community_id == idx]
}


#Run MCMC for Bayesian Community Based Targeting -  WITHOUT CORRECTION
CBtemp <- CBTarget(Tau=Tau, 
                   X_CBT = X_CBT[,-which(colnames(X_CBT) == "connected")],
                   X_program = X_program[,-which(colnames(X_program) == "connected")],
                   X_elite =NULL,
                   iter_keep =iter_keep,
                   iter_burn = iter_burn,
                   print_opt = print_opt)
#Run MCMC for Bayesian Community Based Targeting -  WITH CORRECTION
CBtemp_noelite <- CBTarget(Tau=Tau, 
                           X_CBT = X_CBT,
                           X_program = X_program,
                           X_elite = "connected",
                           iter_keep =iter_keep,
                           iter_burn = iter_burn,
                           print_opt = print_opt)

which_noelite <- which(colnames(X_CBT) == "connected") #NOTE THIS INDEX INCLUDES THE FIRST POSITION OF INTERCEPT

Y_micro <- as.matrix(log(CBT_data$consumption))
Y_micro <- apply(Y_micro, 2, function(x){(x - mean(x))/sd(x)})
X_PMT_sub <-     cbind(1,CBT_data[,m3]) %>%as.matrix()
temp_data <- data.frame(Y_micro = Y_micro,
                        X_CBT)
form <- formula(paste0("Y_micro~", paste0(colnames(X_CBT)[-which_noelite], collapse = "+")))

PMT_beta <-coef(lm(form, data = temp_data))%>%as.vector()

CB_beta_rank_mean_noelite <- apply(CBtemp_noelite$beta_rank, 2, mean)
CB_beta_rank_mean <- append(apply(CBtemp$beta_rank, 2, mean), 0, after = which_noelite-1)

coefs <- data.frame(parameter = m3,
                    CB_beta_rank_mean_noelite = CB_beta_rank_mean_noelite[-1],
                    CB_beta_rank_mean = CB_beta_rank_mean[-1],
                    PMT_beta = append(PMT_beta[-1], 0, after = which_noelite-2)
)
write.csv(coefs, paste0("Indonesia Analysis/coef_elite", elite_status, ".csv"), row.names = FALSE)

