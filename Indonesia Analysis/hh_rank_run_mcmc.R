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


iter_keep = 20   ## Gibbs sampler kept iterations (post burn-in)
iter_burn = 20   ## Gibbs sampler burn-in iterations 
print_opt = 100  ## print a message every print.opt steps


all_ranks <- read.csv("Data/Indonesia/Cleaning/hh_rankings.csv") %>%
  arrange(hhid_ranked)


ranked1 <- unique(all_ranks$hhid_ranked)
rankers1 <- unique(all_ranks$hhid_ranker)

#only contains ranked hh that are present in alatas.csv and hh_rankings.csc
full_data <- read.csv("Data/Indonesia/Cleaning/alatas.csv") %>%
  dplyr::select(-c("hhsize_ae")) %>% arrange(village, province, district, subdistrict)%>% 
  mutate(community_id = as.numeric(factor(interaction(village, province, district, subdistrict))))%>%
  group_by(village, province, district, subdistrict)%>%
  mutate(prop_rank = rank,
         poverty_rate = mean(treated)) %>%
  mutate(rank = ifelse(is.na(rank), NA, floor(rank(rank)))) %>%ungroup %>%
  merge(data.frame(ranked1), by.x = "hhid", by.y = "ranked1", all.x = FALSE, all.y = FALSE)%>%
  arrange(hhid) 


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

#RANDOM SAMPLES OF CBT DATA, PMT DATA, AND PROGRAM (testing) DATA
whats_left <- unique(full_data$community_id[-PMT_idx]) #communities not in PMT

CBT_data <- full_data %>%subset(community_id %in% whats_left)

PMT_data <- CBT_data#full_data[PMT_idx,] #%>% subset(community == 0)


all_ranks <- subset(all_ranks, hhid_ranked %in% unique(CBT_data$hhid))
ranked <- unique(all_ranks$hhid_ranked)
rankers <- unique(all_ranks$hhid_ranker)


all_ranks_sum <- all_ranks %>% group_by(hhid_ranked) %>%
  summarise(rank = mean(rank, na.rm=TRUE))


#Create Rank matrix in format required for CBTarget()
Tau2 <- array(NA, dim = c(length(unique(all_ranks$hhid_ranked)), length(unique(all_ranks$hhid_ranker)))) %>%as.data.frame()
Tau2$they_ranked <- NA
j = 0
for ( idx in rankers){ #loop over columns
  all_ranks_sub <- subset(all_ranks, hhid_ranker == idx)
  who_they_ranked <- all_ranks_sub$hhid_ranked
  Tau2$they_ranked[ranked == idx] <- paste0(paste0("h",who_they_ranked), collapse = "")
  j = j + 1
  colnames(Tau2)[j] <- paste0("h", idx)
  for (idx2 in who_they_ranked){
    Tau2[ranked == idx2, j] <- all_ranks_sub$rank[all_ranks_sub$hhid_ranked == idx2][1] 
    #the [1] is because some households ranked another household twice... i just took the first for now
  }
}

Tau2$ranked_them <- apply(Tau2, 1, function(x){paste0(colnames(Tau2)[!is.na(x)], collapse = "")})
Tau2$group <- rnorm(nrow(Tau2))
Tau2$group[1] <- 1
for (i in 2:nrow(Tau2)){
  if (any(c(
    grepl(as.character(ranked[i]), Tau2$ranked_them[Tau2$group == Tau2$group[i-1] ] ),
    grepl(as.character(ranked[i]), Tau2$they_ranked[Tau2$group == Tau2$group[i-1] ] )))){
    Tau2$group[i] <- Tau2$group[i-1] 
  }else{
    Tau2$group[i] <- Tau2$group[i-1] +1
  }
}

head(Tau2$group, 8)
head(Tau2$ranked_them, 8)
head(ranked, 8)


Tau3 <- data.frame(mean_rank = apply(Tau2[,1:(ncol(Tau2)-3)], 1, mean, na.rm=TRUE), group = Tau2$group) %>%
  group_by(group) %>%
  mutate(rank = floor(rank(mean_rank))) %>% ungroup()

Tau <- array(NA, dim = c(length(ranked),max(Tau3$group)))
for (i in unique(Tau3$group)){
  Tau[Tau3$group == i,i] <- Tau3$rank[Tau3$group == i]
}


rm(Tau2, Tau3)

Program_data <- full_data[1:10,]  

X_PMT <-     cbind(1,PMT_data[,m3]) %>%as.matrix()#cbind(1, PMT_data[,m3]%>%apply(2, function(x){(x - mean(x))/sd(x)})) 
X_CBT <-     cbind(1,CBT_data[,m3]) %>%as.matrix()
X_program <- cbind(1,Program_data[,m3]) %>%as.matrix()

R = ncol(Tau)

#run with elite capture variable 

CBtemp <- CBTarget(Tau=Tau, 
                      X_CBT = X_CBT,
                      X_program = X_program,
                      X_elite ="connected",
                      iter_keep =iter_keep,
                      iter_burn = iter_burn,
                      print_opt = print_opt)





