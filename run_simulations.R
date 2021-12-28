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
library(DirichletReg)
detectCores(logical=FALSE)


sourceCpp("functions.cpp")
source("Bayes Consensus Ranking/functions.R")
source("Bayes Consensus Ranking/HybridTarget.R")
source("Bayes Consensus Ranking/CBTarget.R")

#Run Indonesia simulation study
source("Alatas Analysis/run_mcmc.R")
#Run Burkina Faso simulation study
source("Hillebrecht Analysis/run_mcmc.R")



