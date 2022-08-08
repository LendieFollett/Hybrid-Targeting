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

#Run Indonesia simulation study
source("Indonesia Analysis/run_mcmc.R")

#Output: 
## Indonesia Analysis/all_results.csv
## Indonesia Analysis/all_coef.csv
## Indonesia Analysis/coef_total_sample.csv
## Indonesia Analysis/CB_beta_rank_CI_noelite.csv
## Indonesia Analysis/CB_beta_rank_CI.csv

#Run Burkina Faso simulation study
source("Burkina Faso/run_mcmc.R")

#Output: 
## Burkina Faso Analysis/all_results.csv
## Burkina Faso Analysis/all_coef.csv
## Burkina Faso Analysis/coef_total_sample.csv
## Burkina Faso Analysis/CB_beta_rank_CI_noelite.csv
## Burkina Faso Analysis/CB_beta_rank_CI.csv

