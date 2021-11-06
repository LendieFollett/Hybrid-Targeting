


#SIMULATION CODE ASSUMES EACH COMMUNITY HAS R/K RANKING SYSTEMS FOR THEIR MEMBERS 

community_PMT <- rep(1:K, N_PMT/K) #each training + testing household assigned to a community
community_CBT <- rep(1:K, N_CBT/K) #each training + testing household assigned to a community
community_Program <- rep(1:K, N_Program/K) #each training + testing household assigned to a community

# X_MICRO - both training and testing has same micro covariate information
X_PMT = cbind(rmvnorm(N_PMT, mean = rep(0, P-1), sigma = diag(P-1)), #regular demographic
                rbinom(n = N_PMT, size = 1, p = .1)) #indicator for elite connection
X_PMT = cbind(1, X_PMT) #add intercept

X_CBT = cbind(rmvnorm(N_CBT, mean = rep(0, P-1), sigma = diag(P-1)), #regular demographic
              rbinom(n = N_CBT, size = 1, p = .1)) #indicator for elite connection
X_CBT = cbind(1, X_CBT) #add intercept


X_Program = cbind(rmvnorm(N_Program, mean = rep(0, P-1), sigma = diag(P-1)), #regular demographic
              rbinom(n = N_Program, size = 1, p = .1)) #indicator for elite connection
X_Program <- cbind(1, X_Program) #add intercept

Z <- array(NA, dim = c(N_CBT, R)) #only testing has latent ranks (e.g., consumption)

#parameter values
omega_micro_true <- 0.5
set.seed(4632340)
omega_rank_true <- rep(sample(x=c(.5, 1, 2), replace=TRUE, size = R/K),K)#, prob = c(0,0,1)
set.seed(3527357)
beta_rank_true = 2*c(0,sample(c(0,1), size = P, replace=TRUE))  #first column is intercept
set.seed(68568)
beta_micro_true = beta_rank_true + c(0,rnorm(P, mean = 0, sd = .25))
mu_beta <- apply(cbind(beta_rank_true,beta_micro_true), 1, mean)
sigma2_alpha <- 1
#Fill "responses"

#gamma_micro_true <-  rnorm(N0, 0, sigma2_micro^.5)
Y_micro <-  rnorm(N_PMT, X_PMT %*% beta_micro_true, sqrt(1/omega_micro_true)) #+ gamma_micro_true

set.seed(4296724)
alpha_true <-  rnorm(N_CBT, 0, sigma2_alpha^.5)
for (r in 1:R){ #fill latent Z scores
  Z[,r] <-  rnorm(N_CBT, X_CBT %*% beta_rank_true, sqrt(1/omega_rank_true[r])) +alpha_true
}



#Tau <- apply(Z, 2, rank)

#incomplete rankings
Tau <- data.frame(Z, community=community_CBT) %>% 
  group_by(community) %>%
  mutate_all(rank)  #rank within community
#If there are R rankers and K communities,
#then each community had R/K unique rankers
Ztemp <- data.frame(Z, community=community_CBT)

j=0
for ( idx in 1:K){ #loop over columns
  for (idx2 in 1:(R/K)){
    j = j + 1
  Tau[Tau$community != idx,j] <- NA
  Ztemp[Ztemp$community!=idx,j] <- NA
  }
}

Tau <- subset(Tau, select = -c(community)) %>%
  as.matrix() #R rankings (what we actually observe)
Ztemp <- subset(Ztemp, select = -c(community)) %>%
  as.matrix() #R rankings (what we actually observe)

#real life observe: Tau, Y_micro ('test' only), Y_comm 
#do NOT observe: Z



