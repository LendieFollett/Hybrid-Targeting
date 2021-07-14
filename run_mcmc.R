#Run MCMC, convergence diagnostics
#see results.R for analysis
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

source("Bayes Consensus Ranking/functions.R")
#parameters for simulation

full_data <- read.csv("Empirical Study/alatas.csv")
#add a community id
full_data <- full_data%>% mutate(community_id = as.numeric(factor(interaction(village, province, district, subdistrict))))
#groups of x variables
m1 <- c("elite","hhsize","hhsize_ae","hhage","hhmale","hhmarried","hhage2", "hhsize2", "hhmalemarr",
        "hheduc2","hheduc3","hheduc4",
        "age04","higheduc2","higheduc3","higheduc4","depratio")
m2.1 <- c("pcfloor", "tfloor","twall", "toilet","water","lighting", "troof",
          "fcook","house", "ac","computer","radio","tv", "dvd","satellite", 
          "gas", "refrigerator", "bicycle", "motorcycle", "auto", "hp", 
          "jewelry","chicken","cow")
m2 <- c(m1, m2.1)
m3.1 <- c("credit","hhsector1", "hhsector2","hhsector3",
          "formal","informal", "eschild","jschild","sschild")
m3 <- c(m2,m3.1) #full collection
test_data <- full_data %>% subset(community == 1)
train_data <- full_data %>% subset(community == 0)

X_micro0 <- cbind(1,train_data[,m3]) %>%as.matrix
X_micro1 <- cbind(1, test_data[,m3]) %>%as.matrix

Y_micro <- train_data$consumption

R = test_data %>% group_by(village, province, district, subdistrict) %>% summarise(n = length(cow))%>%ungroup() %>%nrow
M = 1  ## just consumption
N0 = train_data %>%nrow
N1 = test_data %>%nrow
P = ncol(X_micro0)-1 #(-1 since i've included the intercept)


#starting values for random effects
temp_data <- data.frame(y = (Y_micro - mean(Y_micro,na.rm=TRUE))/sd(Y_micro),
                         X_micro0)
form <- formula(paste0("y~-1+", paste0(colnames(X_micro0), collapse = "+")))
#gamma_start <- ranef(lmer(form, data = temp_data))[[1]]$`(Intercept)` 
beta_start <-coef(lm(form, data = temp_data))%>%as.vector()
initial_list <- list(#gamma_rank = gamma_start,
  beta_rank = beta_start,
  beta_micro = beta_start,
  mu_beta = beta_start)

Tau <- array(NA, dim = c(nrow(test_data), R))
j = 0
for ( idx in unique(test_data$community_id)){ #loop over columns
  j = j + 1
  Tau[test_data$community_id == idx,j] <- test_data$rank[test_data$community_id == idx]
}


iter_keep = 1000   ## Gibbs sampler kept iterations (post burn-in)
iter_burn =500   ## Gibbs sampler burn-in iterations 
print_opt = 100  ## print a message every print.opt steps


#Run MCMC for Bayesian Consensus Targeting
temp <- BCTarget(Tau=Tau, 
                 X_micro0 = X_micro0, 
                 X_micro1 = X_micro1,
                 X_elite = "elite",
                 Y_micro = Y_micro,
                 iter_keep = iter_keep,
                 iter_burn = iter_burn,
                  print_opt = print_opt)






