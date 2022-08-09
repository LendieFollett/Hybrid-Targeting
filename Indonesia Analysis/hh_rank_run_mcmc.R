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


iter_keep = 20   ## Gibbs sampler kept iterations (post burn-in)
iter_burn = 20   ## Gibbs sampler burn-in iterations 
print_opt = 100  ## print a message every print.opt steps


all_ranks <- read.csv("Data/Indonesia/Cleaning/hh_rankings.csv") 

#Create Rank matrix in format required for CBTarget()
Tau2 <- array(NA, dim = c(length(unique(all_ranks$hhid_ranked)), length(unique(all_ranks$hhid_ranker))))
j = 0
for ( idx in unique(all_ranks$hhid_ranked)){ #loop over rows
  j = j + 1
  Tau2[all_ranks$hhid_ranker == idx,j] <- all_ranks$rank[all_ranks$hhid_ranker == idx]
}





