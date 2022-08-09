rm(list = ls())
detectCores(logical=FALSE)

iter_keep = 2000   ## Gibbs sampler kept iterations (post burn-in)
iter_burn = 2000   ## Gibbs sampler burn-in iterations 
print_opt = 100  ## print a message every print.opt steps


full_data <- read.csv("Data/Burkina Faso/Cleaning/hillebrecht.csv") %>%
  group_by(community, year)%>%
  mutate(informant1 = ifelse(is.na(informant1), NA, floor(rank(-informant1))),
         informant2 = ifelse(is.na(informant2), NA, floor(rank(-informant2))),
         informant3 = ifelse(is.na(informant3), NA, floor(rank(-informant3))),
         treat_rate = sum(treated)/length(treated)) %>% 
  #aggregation of three rankings for use in poverty rate exercises
 mutate(informant_agg = floor(rank(mean(informant1, informant2, informant3)))) %>%
  ungroup%>%  arrange(community)

#x variables to include in model
m_num <- c("rooms", "hhsize","age1660","age60")
m_bin <- colnames(full_data)[which(colnames(full_data)=="floor"):which(colnames(full_data)=="pig")]
m_bin <- m_bin[!m_bin %in% m_num]
m3 <- c(m_num, m_bin)

#standardize numeric covariates
full_data <- full_data %>%mutate_at(m_num, function(x){(x - mean(x))/(2*sd(x))})

#training data is all of 2008 data
PMT_idx <-which(full_data$year == 2008) 
full_data_left <- full_data[-PMT_idx,]
PMT_data <- full_data[PMT_idx,] 

#-----PRELIMINARY: RUN PHASE 1 FOR DYNAMIC UPDATING
CBT1_data <- PMT_data
R1 = 3*(CBT1_data %>% group_by(community) %>% summarise(n = length(floor))%>%ungroup() %>%nrow)
X_CBT1 <-     cbind(1,CBT1_data[,m3]) %>%as.matrix()

#Create Rank matrix in format required for CBTarget()
Tau1 <- array(NA, dim = c(nrow(CBT1_data), R1))
j = 0
for ( idx in unique(CBT1_data$community)){ #loop over columns
  for (infmt in c("informant1", "informant2", "informant3")){
    j = j + 1
    Tau1[CBT1_data$community == idx,j] <- pull(CBT1_data[CBT1_data$community == idx,], infmt)
  }
}

CBtemp_P1 <- CBTarget(Tau=Tau1, 
                      X_CBT = X_CBT1[,-which(colnames(X_CBT1) == "minority")],
                      X_program = X_CBT1[,-which(colnames(X_CBT1) == "minority")],
                      X_elite =NULL,
                      iter_keep =iter_keep,
                      iter_burn = iter_burn,
                      print_opt = print_opt)

DU_beta_rank_mean <- apply(CBtemp_P1$beta_rank[,-1], 2, mean)
DU_beta_rank_vcov <- cov(CBtemp_P1$beta_rank[,-1])
DU_omega_rank_probs <- c(mean(CBtemp_P1$omega_rank[,1]==0.5),mean(CBtemp_P1$omega_rank[,1]==1),mean(CBtemp_P1$omega_rank[,1]==2))
DU_omega_rank_probs <- apply(cbind(c(1/3,1/3,1/3), DU_omega_rank_probs), 1, mean)




CBT_ncomm_list <- c(5,10,15,25)
nrep <- 30
results <-  mclapply(CBT_ncomm_list, function(CBT_ncomm){
  i <- 0
  r <- list()
  c <- list()
  for(rep in c(1:nrep)){
    print(paste("***********Rep ", rep," of CBT proportion ", CBT_ncomm, "**************"))
    i = i + 1

    
    #RANDOM SAMPLES OF CBT DATA, PMT DATA, AND PROGRAM (testing) DATA
    whats_left <- unique(full_data_left$community) #communities not in PMT
    samps <- data.frame(community = whats_left,
                        samp =     rep(c("CBT2", "Program", "NA"), 
                                       c(CBT_ncomm, 
                                         25, max(51-CBT_ncomm - 25, 0)))[sample.int(length(whats_left))])

    
    CBT2_data <- full_data_left %>%subset(community %in% samps$community[samps$samp == "CBT2"])
    Program_data <- full_data_left %>%subset(community %in% samps$community[samps$samp == "Program"])    

    #print(c(dim(CBT1_data)[1],dim(CBT2_data)[1],dim(Program_data)))
    #print(table(samps$samp))
    
    while(any(apply(CBT2_data[,m3], 2, var) == 0)|any(apply(CBT2_data[,m3], 2, var) == 0)|any(apply(PMT_data[,m3], 2, var) == 0)){ #have to do to deal with complete separation and ML estimation of logistic regression (#shouldadonebayes)
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
    
    initial_list<- list(
      beta_rank = c(0,CBT_beta_start[-c(1, which_noelite)]),
      beta_micro = PMT_beta_start[-which_noelite],
      mu_beta = mu_beta_start[-c(1, which_noelite)])
    
    #create rank matrix: one column per 'ranker' (community)
    
    Tau2 <- array(NA, dim = c(nrow(CBT2_data), R2))
    j = 0
    for ( idx in unique(CBT2_data$community)){ #loop over columns
      for (infmt in c("informant1", "informant2", "informant3")){
        j = j + 1
        Tau2[CBT2_data$community == idx,j] <- pull(CBT2_data[CBT2_data$community == idx,], infmt)
      }
    }
    
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
                                       initial.list = initial_list_noelite)
    
    #Run MCMC for Bayesian Consensus Targeting - WITHOUT CORRECTION
    Hybridtemp <- HybridTarget(Tau=Tau2, 
                               X_PMT = X_PMT[,-which(colnames(X_PMT) == "minority")], 
                               X_CBT = X_CBT2[,-which(colnames(X_CBT2) == "minority")],
                               X_program = X_program[,-which(colnames(X_program) == "minority")],
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
                               X_elite = "minority",
                               iter_keep =iter_keep,
                               iter_burn = iter_burn,
                               print_opt = print_opt,
                               initial.list = initial_list_noelite)
    
    #Run MCMC for Bayesian Community Based Targeting -  WITHOUT CORRECTION
    CBtemp <- CBTarget(Tau=Tau2, 
                       X_CBT = X_CBT2[,-which(colnames(X_CBT2) == "minority")],
                       X_program = X_program[,-which(colnames(X_program) == "minority")],
                       X_elite =NULL,
                       iter_keep =iter_keep,
                       iter_burn = iter_burn,
                       print_opt = print_opt,
                       initial.list = initial_list)
    
    #DYNAMIC UPDATING------------------------------------------
    
    CBtemp_DU <- CBTarget(Tau=Tau2, 
                          X_CBT = X_CBT2[,-which(colnames(X_CBT2) == "minority")],
                          X_program = X_program[,-which(colnames(X_program) == "minority")],
                          X_elite = NULL,
                          iter_keep =iter_keep,
                          iter_burn = iter_burn,
                          print_opt = print_opt,
                          initial.list = initial_list,
                          delta_prior_mean = DU_beta_rank_mean,
                          delta_prior_var = diag(DU_beta_rank_vcov), #feed it the whole vector
                          prior_prob_rank = DU_omega_rank_probs)
    
    #PMT - no elite bias correction

    Y_micro_sub <- as.matrix(ihs_trans(CBT2_data$consumption))
    Y_micro_sub <- apply(Y_micro_sub, 2, function(x){(x - mean(x))/sd(x)})
    X_PMT_sub <-     cbind(1,CBT2_data[,m3]) %>%as.matrix()
    temp_data <- data.frame(Y_micro = Y_micro_sub,
                            X_PMT_sub)
    form <- formula(paste0("Y_micro~", paste0(colnames(X_PMT_sub)[-which_noelite], collapse = "+")))
    
    PMT_model <- lm(form, data = temp_data)
    PMT_beta <-coef(PMT_model)%>%as.vector()
    
    # -----CBT probit-------------------------------------

    lr <- glm(treated ~ ., 
              data = CBT2_data[,c("treated", m3[-which(m3 == "minority")])], 
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
    CB_DU_beta_rank_mean <- append(apply(CBtemp_DU$beta_rank, 2, mean), 0, after = which_noelite-1)
    
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
                         CB_DU_beta_rank_mean = CB_DU_beta_rank_mean[-1],
                         
                         PMT_beta = append(PMT_beta[-1], 0, after = which_noelite-2))
    

    #HYBRID-BASED PREDICTION - WITH CORRECTION
    Program_data$hybrid_prediction_noelite <-apply(Hybridtemp_noelite$mu_noelite, 2, mean)
    #HYBRID-BASED PREDICTION - WITHOUT CORRECTION
    Program_data$hybrid_prediction <-apply(Hybridtemp$mu, 2, mean)
    
    #CBT SCORE-BASED PREDICTION -  WITH CORRECTION
    Program_data$cbt_model_prediction_noelite <- apply(CBtemp_noelite$mu_noelite, 2, mean)
    #CBT SCORE-BASED PREDICTION -  WITHOUT CORRECTION
    Program_data$cbt_model_prediction <- apply(CBtemp$mu, 2, mean)
    #CBT SCORE-BASED PREDICTION -  WITHOUT CORRECTION, WITH DYNAMIC UPDATING
    Program_data$cbt_DU_model_prediction <- apply(CBtemp_DU$mu, 2, mean)
    
    #PROBIT CBT-BASED PREDICTION -  WITHOUT CORRECTION
    Program_data$CBT_LR_prediction<- -predict(lr, as.data.frame(X_program[,-1])) #(LRF NEEDS TO CHANGE TO) logistic regression
    
    #OLS-BASED PMT PREDICTION -  WITHOUT CORRECTION
    Program_data$PMT_prediction <- predict(PMT_model, Program_data)#(X_program[,-c(1, which_noelite)]%*%PMT_beta[-1])#beta_start is the OLS estimate of beta
    
    Program_data <- Program_data%>%group_by(community) %>%
      mutate(hybrid_noelite_rank =rank(hybrid_prediction_noelite)/length(consumption),
             hybrid_rank =rank(hybrid_prediction)/length(consumption),
             pmt_rank =rank(PMT_prediction)/length(consumption),
             cbt_model_rank = rank(cbt_model_prediction)/length(consumption),
             cbt_DU_model_rank = rank(cbt_DU_model_prediction)/length(consumption),
             cbt_model_rank_noelite = rank(cbt_model_prediction_noelite)/length(consumption),
             consumption_rank = rank(consumption)/length(consumption),
             CBT_LR_rank = rank(CBT_LR_prediction)/length(consumption)) 
    
    
 
    
    
    
    r[[i]] <- Program_data %>% dplyr::select(c(hhid, year,community,treat_rate, 
                                                         hybrid_noelite_rank, hybrid_rank,
                                                         cbt_model_rank, cbt_model_rank_noelite,
                                                        cbt_DU_model_rank,
                                                         consumption_rank, cbt_rank,CBT_LR_rank,
                                                         pmt_rank)) %>%
      mutate(rep = rep,
             CBT_ncomm = CBT_ncomm)
    
    
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


write.csv(all_results, "Burkina Faso Analysis/all_results.csv")

write.csv(all_coef, "Burkina Faso Analysis/all_coef.csv")

#---------------------------------------------
#RANDOM SAMPLES OF CBT DATA, PMT DATA, AND PROGRAM (testing) DATA
whats_left <- unique(full_data$community[-PMT_idx]) #communities not in PMT

CBT_data <- full_data %>%subset(community %in% whats_left)

PMT_data <- CBT_data#full_data[PMT_idx,] #%>% subset(community == 0)

Program_data <- full_data[1:10,]  


X_PMT <-     cbind(1,PMT_data[,m3]) %>%as.matrix()#cbind(1, PMT_data[,m3]%>%apply(2, function(x){(x - mean(x))/sd(x)})) 
X_CBT <-     cbind(1,CBT_data[,m3]) %>%as.matrix()
X_program <- cbind(1,Program_data[,m3]) %>%as.matrix()

R = 3*(CBT_data %>% group_by(community) %>% summarise(n = length(community))%>%ungroup() %>%nrow)

Tau <- array(NA, dim = c(nrow(CBT_data), R))
j = 0
for ( idx in unique(CBT_data$community)){ #loop over columns
  for (infmt in c("informant1", "informant2", "informant3")){
    j = j + 1
    Tau[CBT_data$community == idx,j] <- pull(CBT_data[CBT_data$community == idx,], infmt)
  }
}

which_noelite <- which(colnames(X_CBT) == "minority") #NOTE THIS INDEX INCLUDES THE FIRST POSITION OF INTERCEPT


#Run MCMC for Bayesian Community Based Targeting -  WITHOUT CORRECTION
CBtemp <- CBTarget(Tau=Tau, 
                   X_CBT = X_CBT[,-which(colnames(X_CBT) == "minority")],
                   X_program = X_program[,-which(colnames(X_program) == "minority")],
                   X_elite =NULL,
                   iter_keep =iter_keep,
                   iter_burn = iter_burn,
                   print_opt = print_opt)
#Run MCMC for Bayesian Community Based Targeting -  WITH CORRECTION
CBtemp_noelite <- CBTarget(Tau=Tau, 
                           X_CBT = X_CBT,
                           X_program = X_program,
                           X_elite = "minority",
                           iter_keep =iter_keep,
                           iter_burn = iter_burn,
                           print_opt = print_opt)

which_noelite <- which(colnames(X_CBT) == "minority") #NOTE THIS INDEX INCLUDES THE FIRST POSITION OF INTERCEPT

Y_micro <- as.matrix(ihs_trans(CBT_data$consumption))
Y_micro <- apply(Y_micro, 2, function(x){(x - mean(x))/sd(x)})
X_PMT_sub <-     cbind(1,X_CBT[,m3]) %>%as.matrix()
temp_data <- data.frame(Y_micro = Y_micro,
                        X_PMT_sub)
form <- formula(paste0("Y_micro~", paste0(colnames(X_PMT_sub)[-which_noelite], collapse = "+")))
PMT_model <- lm(form, data = temp_data)
PMT_beta <-coef(PMT_model)%>%as.vector()


CB_beta_rank_mean_noelite <- apply(CBtemp_noelite$beta_rank, 2, mean)
CB_beta_rank_mean <- append(apply(CBtemp$beta_rank, 2, mean), 0, after = which_noelite-1)


coefs <- data.frame(parameter = m3,
                    CB_beta_rank_mean_noelite = CB_beta_rank_mean_noelite[-1],
                    CB_beta_rank_mean = CB_beta_rank_mean[-1],
                    PMT_beta = append(PMT_beta[-1], 0, after = which_noelite-2)
)
write.csv(coefs, "Burkina Faso Analysis/coef_total_sample.csv", row.names=FALSE)


data.frame(parameter = m3,
           apply(CBtemp_noelite$beta_rank[,-1], 2, quantile, c(.025, .975)) %>%t(),
           mean = CB_beta_rank_mean_noelite[-1]) %>%
  write.csv( "Burkina Faso Analysis/CB_beta_rank_CI_noelite.csv")

data.frame(parameter = m3[-(which_noelite-1)],
           apply(CBtemp$beta_rank[,-1], 2, quantile, c(.025, .975)) %>%t(),
           mean = CB_beta_rank_mean[-c(1, which_noelite)]) %>%
  write.csv("Burkina Faso Analysis/CB_beta_rank_CI.csv")

apply(CBtemp_noelite$beta_rank, 2, doESS)
apply(CBtemp$beta_rank, 2, doESS)

