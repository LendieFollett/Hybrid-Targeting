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


all_ranks_sum <- all_ranks %>% group_by(hhid_ranked, hhea) %>%
  summarise(rank = mean(rank, na.rm=TRUE)) %>%
  group_by(hhea) %>%
  mutate(rank = floor(rank(rank))) %>%
  arrange(hhid_ranked) 


#Create Rank matrix in format required for CBTarget()
Tau <- array(NA, dim = c(length(unique(all_ranks$hhid_ranked)), length(unique(all_ranks$hhid_ranker)))) %>%as.data.frame()

for (i in 1:length(unique(all_ranks_sum$hhea))){
  hhea <- unique(all_ranks_sum$hhea)[i]
  Tau[CBT_data$hhid %in% all_ranks_sum$hhid_ranked[all_ranks_sum$hhea == hhea],i]<-all_ranks_sum$rank[all_ranks_sum$hhea == hhea]
}



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


#run without it
CBtemp <- CBTarget(Tau=Tau, 
                   X_CBT = X_CBT[,-which(colnames(X_CBT) == "connected")],
                   X_program = X_program[,-which(colnames(X_program) == "connected")],
                   X_elite =NULL,
                   iter_keep =iter_keep,
                   iter_burn = iter_burn,
                   print_opt = print_opt)


which_noelite <- which(colnames(X_CBT) == "connected") #NOTE THIS INDEX INCLUDES THE FIRST POSITION OF INTERCEPT

CB_beta_rank_mean_noelite <- apply(CBtemp_noelite$beta_rank, 2, mean)
CB_beta_rank_mean <- append(apply(CBtemp$beta_rank, 2, mean), 0, after = which_noelite-1)

all_coef_hh <- data.frame(parameter = m3,
                    CB_beta_rank_mean_noelite = CB_beta_rank_mean_noelite[-1],
                    CB_beta_rank_mean = CB_beta_rank_mean[-1])

all_coef_comm <- read.csv("Indonesia Analysis/coef_total_sample.csv")

variable_labels <- read.csv("Data/Indonesia/Cleaning/variables.csv")

#variable_labels_add <- data.frame(Name = "connected", Definition = "Elite connection")

#variable_labels <- rbind(variable_labels, variable_labels_add)
library(lmomco)
score_order <- all_coef %>% merge(variable_labels, by.x = "parameter", by.y = "Name") %>% 
  dplyr::select(Definition,Order, CB_beta_rank_mean) %>%
  #subset(CB_beta_rank_mean != 0)%>% #remove elite connection 0
  melt(id.vars = c("Definition", "Order")) %>%
  mutate(par_est = ifelse(abs(value) < 0.01, 0, value)) %>%
  group_by(Order,variable) %>%
  mutate(std_mean =par_est/harmonic.mean(par_est)$harmean) %>%
  arrange(-Order) 


all_coef %>%merge(variable_labels, by.x = "parameter", by.y = "Name") %>%
  dplyr::select(Definition, Order, CB_beta_rank_mean, CB_beta_rank_mean_noelite) %>%
  #subset(CB_beta_rank_mean != 0)%>% #remove elite connection 0
  melt(id.vars = c("Definition", "Order")) %>%
  group_by(variable) %>%
  mutate(std_mean = value/mean(abs(value))) %>%
  mutate(Definition = factor(Definition, levels = score_order$Definition),
         variable = factor(variable, levels = c("CB_beta_rank_mean", "CB_beta_rank_mean_noelite"),
                           labels = c("Hybrid", "Hybrid-EC")))%>%
  ggplot() + 
  geom_col(aes(x = Definition, y = std_mean, fill = variable ), position = position_dodge(width = 0.5)) +
  #facet_grid(Category~., scales = "free_y")+
  coord_flip() + 
  theme_bw() +
  labs(x = "", y = "Standardized Coefficient Estimate") +
  scale_fill_grey("Method")+
  theme(legend.position = c(.9,.9), 
        legend.box.background = element_rect(colour = "black"))


ggsave("Indonesia Analysis/coef_score_HH_rankings.pdf", width = 12, height = 12)


