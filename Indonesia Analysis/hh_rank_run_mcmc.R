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


all_ranks <- read.csv("Data/Indonesia/Cleaning/hh_rankings.csv") %>%
  arrange(hhid_ranked)

ranked <- unique(all_ranks$hhid_ranked)
rankers <- unique(all_ranks$hhid_ranker)
#Create Rank matrix in format required for CBTarget()
Tau2 <- array(NA, dim = c(length(unique(all_ranks$hhid_ranked)), length(unique(all_ranks$hhid_ranker))))
j = 0
for ( idx in rankers){ #loop over columns
  all_ranks_sub <- subset(all_ranks, hhid_ranker == idx)
  j = j + 1
  for (idx2 in all_ranks_sub$hhid_ranked){
    Tau2[ranked == idx2, j] <- all_ranks_sub$rank[all_ranks_sub$hhid_ranked == idx2][1] 
    #the [1] is because some households ranked another household twice... i just took the first for now
  }
}





