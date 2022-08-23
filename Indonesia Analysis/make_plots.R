rm(list = ls())
library(dplyr)
library(ggplot2)
library(reshape2)
library(lmomco)


all_results_hh <- read.csv("Indonesia Analysis/all_results.csv")

#pmt_corrected <- read.csv("Indonesia Analysis/PMT_nonconverge_corrected_cbt.csv") %>%
#  mutate(variable = as.factor("EER"), Method = as.factor("PMT")) %>%
#  rename(mean=mean_EER ) %>%
#  relocate(Method, CBT_ncomm, variable, mean)

all_coef <- read.csv("Indonesia Analysis/coef_total_sample.csv")

elite1_coef <- read.csv("Indonesia Analysis/coef_elite1.csv")

elite0_coef <- read.csv("Indonesia Analysis/coef_elite0.csv")

#Sensitivity= P(beneficiary | true poor)
#Specificity= P(non-beneficiary | true non-poor)
#Precision = P(true poor|beneficiary)
#1-Precision = P(true non-poor|beneficiary)
#EE = E2/P = C/(A+C) = 1-sensitivity
#IE = E1/B = B/(A+B) = 1-precision
#TD = C1/P - E1/NP = A/(A+C)-B/(B+D)

#Hybrid (original proposed hybrid)
#"Hybrid-EC" (i.e., hybrid + elite capture adjustment)
#"Hybrid-AI" (i.e., hybrid + auxiliary information) 
#"Hybrid-DU" (i.e., hybrid + dynamic updating)

#vary poverty rate .2, .3, .4
PR <- 0.2
#multiplicative constant shifts community-level poverty rate up or down
multiplicative_constant <- PR/0.3


#calculate inclusions based on chosen poverty rate
all_results_hh <- all_results_hh %>%
  mutate(hybrid_noelite_inclusion = hybrid_noelite_rank <= poverty_rate*multiplicative_constant,
         hybrid_inclusion = hybrid_rank <= poverty_rate*multiplicative_constant,
         pmt_inclusion = pmt_rank <= poverty_rate*multiplicative_constant,
         consumption_inclusion = consumption_rank<=poverty_rate*multiplicative_constant,
         cbt_model_inclusion = cbt_model_rank<=poverty_rate*multiplicative_constant,
         cbt_model_noelite_inclusion = cbt_model_rank_noelite<=poverty_rate*multiplicative_constant,
         CBT_LR_inclusion = CBT_LR_rank<=poverty_rate*multiplicative_constant,
         cbt_inclusion = cbt_rank <= poverty_rate*multiplicative_constant) %>%ungroup() %>%
  mutate_at(vars(matches("inclusion")), as.factor)

#### --- PREP DATA ----------------------------------


#across community analysis
r <- list()
i = 0
for (reps in unique(all_results_hh$rep)){
  for (n in unique(all_results_hh$CBT_ncomm)){
    print(paste0("sample size = ",n, ", rep = ", reps))
    i = i + 1
    all_results_sub <- subset(all_results_hh, rep == reps & CBT_ncomm == n)
    r[[i]] <- rbind(
      confusionMatrix(all_results_sub$hybrid_noelite_inclusion,   all_results_sub$cbt_inclusion,positive = "TRUE")$byClass,
      confusionMatrix(all_results_sub$hybrid_inclusion,           all_results_sub$cbt_inclusion,positive = "TRUE")$byClass,
      confusionMatrix(all_results_sub$cbt_model_noelite_inclusion,all_results_sub$cbt_inclusion,positive = "TRUE")$byClass,
      confusionMatrix(all_results_sub$cbt_model_inclusion,        all_results_sub$cbt_inclusion,positive = "TRUE")$byClass,
      confusionMatrix(all_results_sub$pmt_inclusion,              all_results_sub$cbt_inclusion,positive = "TRUE")$byClass,
      confusionMatrix(all_results_sub$CBT_LR_inclusion,           all_results_sub$cbt_inclusion,positive = "TRUE")$byClass) %>%as.data.frame%>%
      mutate(Method = c( "Hybrid Score (corrected)","Hybrid Score","CBT Score (corrected)", "CBT Score", "PMT OLS", "CBT Logit"),
             rep = reps, 
             CBT_ncomm = n)
    
    r[[i]]$spearman <- c(cor.test(x=all_results_sub$cbt_rank, y=all_results_sub$hybrid_noelite_rank, method = 'spearman')$estimate,
                         cor.test(x=all_results_sub$cbt_rank, y=all_results_sub$hybrid_rank, method = 'spearman')$estimate,
                         cor.test(x=all_results_sub$cbt_rank, y=all_results_sub$cbt_model_rank_noelite, method = 'spearman')$estimate,
                         cor.test(x=all_results_sub$cbt_rank, y=all_results_sub$cbt_model_rank, method = 'spearman')$estimate,
                         cor.test(x=all_results_sub$cbt_rank, y=all_results_sub$pmt_rank, method = 'spearman')$estimate,
                         cor.test(x=all_results_sub$cbt_rank, y=all_results_sub$CBT_LR_rank, method = 'spearman')$estimate)
    
  }
}
#across communities - get one row per rep
all_results <- do.call(rbind, r)



#the loop below takes a while (maybe an hour or so?)
#you can read the csv instead:
all_results_comm <- read.csv(paste0("Indonesia Analysis/all_results_community_level",PR*100,".csv"))
all_results_comm <- read.csv(paste0("Indonesia Analysis/all_results_community_level.csv"))

#within community analysis - get one row per rep per community sampled in test
r <- list()
i = 0
for (reps in unique(all_results_hh$rep)){
  for (n in unique(all_results_hh$CBT_ncomm)){
    print(paste0("sample size = ",n, ", rep = ", reps))
    temp <- subset(all_results_hh, rep == reps & CBT_ncomm == n)
    for (comm in unique(interaction(temp$village, temp$province, temp$district, temp$subdistrict))){
    i = i + 1
    all_results_sub <- subset(all_results_hh, rep == reps & CBT_ncomm == n & 
                                interaction(all_results_hh$village, all_results_hh$province, all_results_hh$district, all_results_hh$subdistrict) == comm)
    r[[i]] <- rbind(
      confusionMatrix(all_results_sub$hybrid_noelite_inclusion,   all_results_sub$cbt_inclusion,positive = "TRUE")$byClass,
      confusionMatrix(all_results_sub$hybrid_inclusion,           all_results_sub$cbt_inclusion,positive = "TRUE")$byClass,
      confusionMatrix(all_results_sub$cbt_model_noelite_inclusion,all_results_sub$cbt_inclusion,positive = "TRUE")$byClass,
      confusionMatrix(all_results_sub$cbt_model_inclusion,        all_results_sub$cbt_inclusion,positive = "TRUE")$byClass,
      #NOTE: PMT being compared to consumption-based truth, not CBT-based truth
      confusionMatrix(all_results_sub$pmt_inclusion,              all_results_sub$consumption_inclusion,positive = "TRUE")$byClass,
      confusionMatrix(all_results_sub$CBT_LR_inclusion,           all_results_sub$cbt_inclusion,positive = "TRUE")$byClass) %>%as.data.frame%>%
      mutate(Method = c( "Hybrid Score (corrected)","Hybrid Score","CBT Score (corrected)", "CBT Score", "PMT OLS", "CBT Logit"),
             rep = reps, 
             CBT_ncomm = n,
             village=all_results_sub$village[1], province=all_results_sub$province[1], district = all_results_sub$district[1], subdistrict = all_results_sub$subdistrict[1])
    
    r[[i]]$spearman <- c(cor.test(x=all_results_sub$cbt_rank, y=all_results_sub$hybrid_noelite_rank, method = 'spearman')$estimate,
                         cor.test(x=all_results_sub$cbt_rank, y=all_results_sub$hybrid_rank, method = 'spearman')$estimate,
                         cor.test(x=all_results_sub$cbt_rank, y=all_results_sub$cbt_model_rank_noelite, method = 'spearman')$estimate,
                         cor.test(x=all_results_sub$cbt_rank, y=all_results_sub$cbt_model_rank, method = 'spearman')$estimate,
                         cor.test(x=all_results_sub$cbt_rank, y=all_results_sub$pmt_rank, method = 'spearman')$estimate,
                         cor.test(x=all_results_sub$cbt_rank, y=all_results_sub$CBT_LR_rank, method = 'spearman')$estimate)
    
      }
    }
}

#across communities - get one row per rep per community sampled in test
all_results_comm <- do.call(rbind, r)

write.csv(all_results_comm, paste0("Indonesia Analysis/all_results_community_level",PR*100,".csv"))
#the following plots show overall rep-level variability

#### --- ERROR RATE PLOTS ----------------------------------

plot_data <- all_results %>%  mutate(IER = 1-Precision,
                        EER = 1-Sensitivity) %>%
  melt(id.var = c("Method", "CBT_ncomm")) %>%
  mutate(Method = factor(Method, levels = c("Hybrid Score (corrected)","Hybrid Score","CBT Score", "CBT Score (corrected)","CBT DU", "CBT Logit", "PMT OLS"),
                         labels = c("Hybrid-AI-EC","Hybrid-AI","Hybrid","Hybrid-EC","Hybrid-DU", "Probit", "PMT"))) %>%
  group_by(Method, CBT_ncomm, variable) %>%
  summarise(mean = mean(value ))%>%ungroup %>%
  subset(variable %in% c( "EER"))

#plot_data[plot_data$Method == "PMT",] <- pmt_corrected


plot_data %>%
  subset( Method %in% c("Hybrid", "PMT", "Probit")  )%>%
  ggplot() + geom_line(aes(x = CBT_ncomm, y = mean, linetype = Method)) +
  geom_point(aes(x = CBT_ncomm, y = mean)) +
  #geom_linerange(aes(x = CBT_ncomm, ymin = min,ymax=max, linetype = Method))+
  theme_bw() +
  labs(x = "Number of Ranking Communities", y = "Average Error Rate")+ 
  theme(legend.position = c(0.9, 0.8)) +
  theme(legend.box.background = element_rect(colour = "black"))

ggsave(paste0("Indonesia Analysis/ER_hybrid",PR*100,".pdf"), width = 8, height = 5)


plot_data %>%
  subset( Method %in% c("Hybrid", "Hybrid-AI")  )%>%
  ggplot() + geom_line(aes(x = CBT_ncomm, y = mean, linetype = Method)) +
  geom_point(aes(x = CBT_ncomm, y = mean)) +
  #geom_linerange(aes(x = CBT_ncomm, ymin = min,ymax=max, linetype = Method))+
  theme_bw() +
  labs(x = "Number of Ranking Communities", y = "Average Error Rate")+ 
  theme(legend.position = c(0.9, 0.8))+
  theme(legend.box.background = element_rect(colour = "black"))

ggsave(paste0("Indonesia Analysis/ER_hybrid_AI",PR*100,".pdf"), width = 8, height = 5)

plot_data %>%
  subset( Method %in% c("Hybrid", "Hybrid-EC")  )%>%
  ggplot() + geom_line(aes(x = CBT_ncomm, y = mean, linetype = Method)) +
  geom_point(aes(x = CBT_ncomm, y = mean)) +
  #geom_linerange(aes(x = CBT_ncomm, ymin = min,ymax=max, linetype = Method))+
  theme_bw() +
  labs(x = "Number of Ranking Communities", y = "Average Error Rate")+ 
  theme(legend.position = c(0.9, 0.8))+
  theme(legend.box.background = element_rect(colour = "black"))

ggsave(paste0("Indonesia Analysis/ER_hybrid_EC",PR*100,".pdf"), width = 8, height = 5)

#(no DU for Indonesia)

#### --- CORRELATION PLOTS - COMMUNITY LEVEL ----------------------------------
#RANK CORRELATION ANALYSIS

plot_data_corr <- all_results_comm %>%  mutate(IER = 1-Precision,
                        EER = 1-Sensitivity) %>%
  melt(id.var = c("Method", "CBT_ncomm")) %>%
  mutate(Method = factor(Method, levels = c("Hybrid Score (corrected)","Hybrid Score","CBT Score", "CBT Score (corrected)","CBT DU", "CBT Logit", "PMT OLS"),
                         labels = c("Hybrid-AI-EC","Hybrid-AI","Hybrid","Hybrid-EC","Hybrid-DU", "Probit", "PMT"))) %>%
  group_by(Method, CBT_ncomm, variable) %>%
  summarise(mean = mean(value ))%>%ungroup %>%
  subset(variable %in% c( "spearman"))


plot_data_corr %>%
  subset( Method %in% c("Hybrid", "PMT", "Probit")  )%>%
  ggplot() + geom_line(aes(x = CBT_ncomm, y = mean, linetype = Method)) +
  geom_point(aes(x = CBT_ncomm, y = mean)) +
  #geom_linerange(aes(x = CBT_ncomm, ymin = min,ymax=max, linetype = Method))+
  theme_bw() +
  labs(x = "Number of Ranking Communities", y = "Average Rank Correlation")+ 
  theme(legend.position = c(0.9, 0.8)) +
  theme(legend.box.background = element_rect(colour = "black"))

ggsave("Indonesia Analysis/CORR_hybrid.pdf", width = 8, height = 5)


#### --- ERROR RATE PLOTS - COMMUNITY VARIABILITY ----------------------------------

#PMT comparison needs to be predicting

plot_data <- all_results_comm %>%  mutate(IER = 1-Precision,
                                     EER = 1-Sensitivity) %>%
  melt(id.var = c("Method", "CBT_ncomm", "village", "province", "district", "subdistrict")) %>%
  mutate(Method = factor(Method, levels = c("Hybrid Score (corrected)","Hybrid Score","CBT Score", "CBT Score (corrected)","CBT DU", "CBT Logit", "PMT OLS"),
                         labels = c("Hybrid-AI-EC","Hybrid-AI","Hybrid","Hybrid-EC","Hybrid-DU", "Probit", "PMT"))) 



plot_data %>% subset(variable == "EER" & CBT_ncomm %in% c(10, 200) & Method %in% c("Hybrid", "Probit", "PMT" )) %>%
  ggplot() +
  geom_boxplot(aes(x = Method,y = value,  colour = as.factor(CBT_ncomm)),draw_quantiles = c(.25, .5, .75)) +
  scale_colour_grey("Number of\nRanking\nCommunities") +
  theme_bw() +
  labs(x = "Method", y = "Community-level Error Rate")

ggsave(paste0("Indonesia Analysis/ER_community_level",PR*100,".pdf"), width = 8, height = 5)

plot_data %>% subset(variable == "EER" & CBT_ncomm %in% c(10, 200)) %>%
  group_by(CBT_ncomm,Method) %>% 
  summarise(mean = mean(value, na.rm=TRUE),
            var = var(value, na.rm=TRUE)) %>%
  write.csv(paste0("Indonesia Analysis/ER_community_level",PR*100,"csv"))

plot_data %>% subset(variable == "EER" & Method %in% c("Hybrid", "PMT")) %>%
  group_by(CBT_ncomm,Method) %>% 
  summarise(mean = mean(value, na.rm=TRUE),
            var = var(value, na.rm=TRUE)) %>%
  ggplot() +
  geom_point(aes(x = CBT_ncomm, y = var))+
  geom_line(aes(x = CBT_ncomm, y = var, linetype = Method)) +
  theme_bw()  +
  labs(x = "Number of Ranking Communities", y = "Variance of Community-Level Error Rates")
#ggsave(paste0("Indonesia Analysis/ER_commlevel",PR*100,".pdf"), width = 8, height = 5)


dodge <- position_dodge(width=.9)
plot_data %>% subset(variable == "EER" & Method %in% c("Hybrid", "PMT")) %>%
  group_by(CBT_ncomm,Method) %>% 
  summarise(mean = mean(value, na.rm=TRUE),
            sd = sd(value, na.rm=TRUE)) %>%
  mutate(CBT_ncomm  = as.factor(CBT_ncomm)) %>%
  ggplot(aes(x = CBT_ncomm)) +
  geom_col(aes( y = mean, fill = Method), position = dodge) +
  geom_errorbar(aes( ymin = mean - sd, ymax = mean + sd, colour = Method),position = dodge, width = 0.25) +
  scale_colour_grey(start = .1, end = .12) +
  scale_fill_grey(start = .3, end = .8) +
  theme_bw()+
  labs(x = "Number of Ranking Communities", y = "Community-Level Error Rates")
ggsave(paste0("Indonesia Analysis/ER_commlevel",PR*100,".pdf"), width = 8, height = 5)



#### --- COEFFICIENT PLOTS ----------------------------------

variable_labels <- read.csv("Data/Indonesia/Cleaning/variables.csv")

#variable_labels_add <- data.frame(Name = "connected", Definition = "Elite connection")

#variable_labels <- rbind(variable_labels, variable_labels_add)

score_order <- all_coef %>% merge(variable_labels, by.x = "parameter", by.y = "Name") %>% 
  dplyr::select(Definition,Order, CB_beta_rank_mean) %>%
  subset(CB_beta_rank_mean != 0)%>% #remove elite connection 0
  melt(id.vars = c("Definition", "Order")) %>%
  mutate(par_est = ifelse(abs(value) < 0.01, 0, value)) %>%
  group_by(Order,variable) %>%
  mutate(std_mean =par_est/harmonic.mean(par_est)$harmean) %>%
  arrange(-Order) 


all_coef %>%merge(variable_labels, by.x = "parameter", by.y = "Name") %>% 
  dplyr::select(Definition, Order, CB_beta_rank_mean, PMT_beta) %>%
  subset(CB_beta_rank_mean != 0)%>% #remove elite connection 0
  melt(id.vars = c("Definition", "Order")) %>%
  group_by(variable) %>%
  mutate(std_mean = value/mean(abs(value)),
         mean_abs = mean(abs(value))) %>%
  mutate(Definition = factor(Definition, levels = score_order$Definition),
         variable = factor(variable, levels = c("CB_beta_rank_mean", "PMT_beta"),
                           labels = c("Hybrid", "PMT")))%>%
  ggplot() + 
  geom_col(aes(x = Definition, y = std_mean, fill = variable ), position = position_dodge(width = 0.5)) +
  #facet_grid(Category~., scales = "free_y")+
  coord_flip() + 
  theme_bw() +
  labs(x = "", y = "Standardized Coefficient Estimate") +
  scale_fill_grey("Method")+
  theme(legend.position = c(.9,.9), 
        legend.box.background = element_rect(colour = "black"))

ggsave("Indonesia Analysis/coef_score_alatas.pdf", width = 12, height = 12)

all_coef %>%merge(variable_labels, by.x = "parameter", by.y = "Name") %>%
  dplyr::select(Definition, Order, CB_beta_rank_mean, CB_beta_rank_mean_noelite) %>%
  subset(CB_beta_rank_mean != 0)%>% #remove elite connection 0
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

ggsave("Indonesia Analysis/coef_score_EC_alatas.pdf", width = 12, height = 12)




#### COEFFICIENT PLOTS WITH ELITE = 1 -------------

elite1_coef %>%merge(variable_labels, by.x = "parameter", by.y = "Name") %>% 
  dplyr::select(Definition, Order, CB_beta_rank_mean, PMT_beta) %>%
  subset(CB_beta_rank_mean != 0)%>% #remove elite connection 0
  melt(id.vars = c("Definition", "Order")) %>%
  group_by(variable) %>%
  mutate(std_mean = value/mean(abs(value)),
         mean_abs = mean(abs(value))) %>%
  mutate(Definition = factor(Definition, levels = score_order$Definition),
         variable = factor(variable, levels = c("CB_beta_rank_mean", "PMT_beta"),
                           labels = c("Hybrid", "PMT")))%>%
  ggplot() + 
  geom_col(aes(x = Definition, y = std_mean, fill = variable ), position = position_dodge(width = 0.5)) +
  #facet_grid(Category~., scales = "free_y")+
  coord_flip() + 
  theme_bw() +
  labs(x = "", y = "Standardized Coefficient Estimate") +
  scale_fill_grey("Method")+
  theme(legend.position = c(.9,.9), 
        legend.box.background = element_rect(colour = "black"))

ggsave("Indonesia Analysis/coef_score_elite1_alatas.pdf", width = 12, height = 12)

elite1_coef %>%merge(variable_labels, by.x = "parameter", by.y = "Name") %>%
dplyr::select(Definition, Order, CB_beta_rank_mean, CB_beta_rank_mean_noelite) %>%
  subset(CB_beta_rank_mean != 0)%>% #remove elite connection 0
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

ggsave("Indonesia Analysis/coef_score_EC_elite1_alatas.pdf", width = 12, height = 12)


#### COEFFICIENT PLOTS WITH ELITE = 0 -------------

elite0_coef %>%merge(variable_labels, by.x = "parameter", by.y = "Name") %>% 
  dplyr::select(Definition, Order, CB_beta_rank_mean, PMT_beta) %>%
  subset(CB_beta_rank_mean != 0)%>% #remove elite connection 0
  melt(id.vars = c("Definition", "Order")) %>%
  group_by(variable) %>%
  mutate(std_mean = value/mean(abs(value)),
         mean_abs = mean(abs(value))) %>%
  mutate(Definition = factor(Definition, levels = score_order$Definition),
         variable = factor(variable, levels = c("CB_beta_rank_mean", "PMT_beta"),
                           labels = c("Hybrid", "PMT")))%>%
  ggplot() + 
  geom_col(aes(x = Definition, y = std_mean, fill = variable ), position = position_dodge(width = 0.5)) +
  #facet_grid(Category~., scales = "free_y")+
  coord_flip() + 
  theme_bw() +
  labs(x = "", y = "Standardized Coefficient Estimate") +
  scale_fill_grey("Method")+
  theme(legend.position = c(.9,.9), 
        legend.box.background = element_rect(colour = "black"))

ggsave("Indonesia Analysis/coef_score_elite0_alatas.pdf", width = 12, height = 12)

elite0_coef %>%merge(variable_labels, by.x = "parameter", by.y = "Name") %>%
  dplyr::select(Definition, Order, CB_beta_rank_mean, CB_beta_rank_mean_noelite) %>%
  subset(CB_beta_rank_mean != 0)%>% #remove elite connection 0
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

ggsave("Indonesia Analysis/coef_score_EC_elite0_alatas.pdf", width = 12, height = 12)


#include elite1 and elite0 in same plot instead of hybrid and hybrid ec
#use hybrid ec for both

score_order <- all_coef %>% merge(variable_labels, by.x = "parameter", by.y = "Name") %>% 
  dplyr::select(Definition,Order, CB_beta_rank_mean) %>%
  #subset(CB_beta_rank_mean != 0)%>% #remove elite connection 0
  melt(id.vars = c("Definition", "Order")) %>%
  mutate(par_est = ifelse(abs(value) < 0.01, 0, value)) %>%
  group_by(Order,variable) %>%
  mutate(std_mean =par_est/harmonic.mean(par_est)$harmean) %>%
  arrange(-Order) 


merge(elite1_coef[,c(1,2)], elite0_coef[,c(1,2)], by = "parameter") %>%
  merge(variable_labels, by.x = "parameter", by.y = "Name") %>%
  dplyr::select(Definition, Order, CB_beta_rank_mean_noelite.x, CB_beta_rank_mean_noelite.y) %>%
  #subset(CB_beta_rank_mean != 0)%>% #remove elite connection 0
  melt(id.vars = c("Definition", "Order")) %>%
  group_by(variable) %>%
  mutate(std_mean = value/mean(abs(value))) %>%
  mutate(Definition = factor(Definition, levels = score_order$Definition),
         variable = factor(variable, levels = c("CB_beta_rank_mean_noelite.x", "CB_beta_rank_mean_noelite.y"),
                           labels = c("Hybrid-EC with Elite", "Hybrid-EC without Elite")))%>%
  ggplot() + 
  geom_col(aes(x = Definition, y = std_mean, fill = variable ), position = position_dodge(width = 0.5)) +
  #facet_grid(Category~., scales = "free_y")+
  coord_flip() + 
  theme_bw() +
  labs(x = "", y = "Standardized Coefficient Estimate") +
  scale_fill_grey("Method")+
  theme(legend.position = c(.9,.9), 
        legend.box.background = element_rect(colour = "black"))
  
ggsave("Indonesia Analysis/coef_score_EC_elitecomp_alatas.pdf", width = 12, height = 12)



