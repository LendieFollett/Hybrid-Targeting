library(dplyr)
library(ggplot2)
library(reshape2)

all_results_hh <- read.csv("Burkina Faso Analysis/all_results.csv")

#pmt_corrected <- read.csv("Burkina Faso Analysis/PMT_nonconverge_corrected_cbt.csv") %>%
#  mutate(variable = as.factor("EER"), Method = as.factor("PMT")) %>%
#  rename(mean=mean_EER ) %>%
#  relocate(Method, CBT_ncomm, variable, mean)

all_coef <- read.csv("Burkina Faso Analysis/coef_total_sample.csv")


#vary poverty rate .2, .3, .4
PR <- 0.3
#multiplicative constant shifts community-level poverty rate up or down
multiplicative_constant <- PR/0.3


all_results_hh <- all_results_hh %>%
  mutate(hybrid_noelite_inclusion = hybrid_noelite_rank <= treat_rate*multiplicative_constant,
         hybrid_inclusion = hybrid_rank <= treat_rate*multiplicative_constant,
         pmt_inclusion = pmt_rank <= treat_rate*multiplicative_constant,
         consumption_inclusion = consumption_rank<=treat_rate*multiplicative_constant,
         cbt_model_inclusion = cbt_model_rank<=treat_rate*multiplicative_constant,
         cbt_DU_model_inclusion = cbt_DU_model_rank <= treat_rate*multiplicative_constant,
         cbt_model_noelite_inclusion = cbt_model_rank_noelite<=treat_rate*multiplicative_constant,
         CBT_LR_inclusion = CBT_LR_rank<=treat_rate*multiplicative_constant,
         cbt_inclusion = cbt_rank <= treat_rate*multiplicative_constant)%>%#,
         #cbt_inclusion = ifelse(treated == 1, TRUE, FALSE)) 
  ungroup() %>%
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
      confusionMatrix(all_results_sub$cbt_DU_model_inclusion,     all_results_sub$cbt_inclusion,positive = "TRUE")$byClass,
      confusionMatrix(all_results_sub$pmt_inclusion,              all_results_sub$cbt_inclusion,positive = "TRUE")$byClass,
      confusionMatrix(all_results_sub$CBT_LR_inclusion,           all_results_sub$cbt_inclusion,positive = "TRUE")$byClass)%>%as.data.frame%>%
      mutate(Method = c( "Hybrid Score (corrected)",
                         "Hybrid Score",
                         "CBT Score (corrected)", 
                         "CBT Score","CBT DU", "PMT OLS", "CBT Logit"),
             rep = reps, 
             CBT_ncomm = n)
    
    
    r[[i]]$spearman <- c(cor.test(x=all_results_sub$cbt_rank, y=all_results_sub$hybrid_noelite_rank, method = 'spearman')$estimate,
                         cor.test(x=all_results_sub$cbt_rank, y=all_results_sub$hybrid_rank, method = 'spearman')$estimate,
                         cor.test(x=all_results_sub$cbt_rank, y=all_results_sub$cbt_model_rank_noelite, method = 'spearman')$estimate,
                         cor.test(x=all_results_sub$cbt_rank, y=all_results_sub$cbt_model_rank, method = 'spearman')$estimate,
                         cor.test(x=all_results_sub$cbt_rank, y=all_results_sub$cbt_DU_model_rank, method = 'spearman')$estimate,
                         cor.test(x=all_results_sub$cbt_rank, y=all_results_sub$pmt_rank, method = 'spearman')$estimate,
                         cor.test(x=all_results_sub$cbt_rank, y=all_results_sub$CBT_LR_rank, method = 'spearman')$estimate)
    
  }
}
#across communities - get one row per rep
all_results <- do.call(rbind, r)


#within community analysis - get one row per rep per community sampled in test

#some communities are too small to compute error rates: 122 and 4, in particular
all_results_hh <- subset(all_results_hh, !community %in% c(4, 122))

r <- list()
i = 0
for (reps in unique(all_results_hh$rep)){
  for (n in unique(all_results_hh$CBT_ncomm)){
    print(paste0("sample size = ",n, ", rep = ", reps))
    temp <- subset(all_results_hh, rep == reps & CBT_ncomm == n)
    for (comm in unique(temp$community)){
      i = i + 1
      all_results_sub <- subset(all_results_hh, rep == reps & CBT_ncomm == n & 
                                  all_results_hh$community == comm)
      r[[i]] <- rbind(
        confusionMatrix(all_results_sub$hybrid_noelite_inclusion,   all_results_sub$cbt_inclusion,positive = "TRUE")$byClass,
        confusionMatrix(all_results_sub$hybrid_inclusion,           all_results_sub$cbt_inclusion,positive = "TRUE")$byClass,
        confusionMatrix(all_results_sub$cbt_model_noelite_inclusion,all_results_sub$cbt_inclusion,positive = "TRUE")$byClass,
        confusionMatrix(all_results_sub$cbt_model_inclusion,        all_results_sub$cbt_inclusion,positive = "TRUE")$byClass,
        #NOTE: PMT being compared to consumption-based truth, not CBT-based truth
        confusionMatrix(all_results_sub$cbt_DU_model_inclusion,     all_results_sub$cbt_inclusion,positive = "TRUE")$byClass,
        confusionMatrix(all_results_sub$pmt_inclusion,              all_results_sub$consumption_inclusion,positive = "TRUE")$byClass,
        confusionMatrix(all_results_sub$CBT_LR_inclusion,           all_results_sub$cbt_inclusion,positive = "TRUE")$byClass) %>%as.data.frame%>%
        mutate(Method = c("Hybrid Score (corrected)","Hybrid Score","CBT Score", "CBT Score (corrected)","CBT DU","PMT OLS", "CBT Logit"),
               rep = reps, 
               CBT_ncomm = n,
               community = all_results_sub$community[1])
      
      r[[i]]$spearman <- c(cor.test(x=all_results_sub$cbt_rank, y=all_results_sub$hybrid_noelite_rank, method = 'spearman')$estimate,
                           cor.test(x=all_results_sub$cbt_rank, y=all_results_sub$hybrid_rank, method = 'spearman')$estimate,
                           cor.test(x=all_results_sub$cbt_rank, y=all_results_sub$cbt_model_rank_noelite, method = 'spearman')$estimate,
                           cor.test(x=all_results_sub$cbt_rank, y=all_results_sub$cbt_model_rank, method = 'spearman')$estimate,
                           cor.test(x=all_results_sub$cbt_rank, y=all_results_sub$cbt_DU_model_rank, method = 'spearman')$estimate,
                           cor.test(x=all_results_sub$cbt_rank, y=all_results_sub$pmt_rank, method = 'spearman')$estimate,
                           cor.test(x=all_results_sub$cbt_rank, y=all_results_sub$CBT_LR_rank, method = 'spearman')$estimate)
      
    }
  }
}
#across communities - get one row per rep per community sampled in test
all_results_comm <- do.call(rbind, r)




#### --- ERROR RATE PLOTS ----------------------------------

plot_data <- all_results %>%  mutate(IER = 1-Precision,
                                     EER = 1-Sensitivity) %>%
  melt(id.var = c("Method", "CBT_ncomm")) %>%
  mutate(Method = factor(Method, levels = c("Hybrid Score (corrected)","Hybrid Score","CBT Score", "CBT Score (corrected)","CBT DU", "CBT Logit", "PMT OLS"),
                         labels = c(        "Hybrid-AI-EC",            "Hybrid-AI",    "Hybrid",   "Hybrid-EC",           "Hybrid-DU", "Probit", "PMT"))) %>%
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
  theme(legend.position = c(0.9, 0.65))+
  theme(legend.box.background = element_rect(colour = "black"))

ggsave(paste0("Burkina Faso Analysis/ER_hybrid",PR*100,".pdf"), width = 8, height = 5)


plot_data %>%
  subset( Method %in% c("Hybrid", "Hybrid-AI")  )%>%
  ggplot() + geom_line(aes(x = CBT_ncomm, y = mean, linetype = Method)) +
  geom_point(aes(x = CBT_ncomm, y = mean)) +
  #geom_linerange(aes(x = CBT_ncomm, ymin = min,ymax=max, linetype = Method))+
  theme_bw() +
  labs(x = "Number of Ranking Communities", y = "Average Error Rate")+ 
  theme(legend.position = c(0.9, 0.8))+
  theme(legend.box.background = element_rect(colour = "black"))

ggsave(paste0("Burkina Faso Analysis/ER_hybrid_AI",PR*100,".pdf"), width = 8, height = 5)

plot_data %>%
  subset( Method %in% c("Hybrid", "Hybrid-DU")  )%>%
  ggplot() + geom_line(aes(x = CBT_ncomm, y = mean, linetype = Method)) +
  geom_point(aes(x = CBT_ncomm, y = mean)) +
  #geom_linerange(aes(x = CBT_ncomm, ymin = min,ymax=max, linetype = Method))+
  theme_bw() +
  labs(x = "Number of Ranking Communities", y = "Average Error Rate")+ 
  theme(legend.position = c(0.9, 0.8))+
  theme(legend.box.background = element_rect(colour = "black"))


ggsave(paste0("Burkina Faso Analysis/ER_hybrid_DU",PR*100,".pdf"), width = 8, height = 5)


#### --- CORRELATION PLOTS - COMMUNITY LEVEL ----------------------------------
#RANK CORRELATION ANALYSIS

plot_data_corr <- all_results_comm %>%  mutate(IER = 1-Precision,
                                          EER = 1-Sensitivity) %>%
  melt(id.var = c("Method", "CBT_ncomm")) %>%
  mutate(Method = factor(Method, levels = c("Hybrid Score (corrected)","Hybrid Score","CBT Score", "CBT Score (corrected)","CBT DU", "CBT Logit", "PMT OLS"),
                         labels = c(        "Hybrid-AI-EC",            "Hybrid-AI",    "Hybrid",   "Hybrid-EC",           "Hybrid-DU", "Probit", "PMT"))) %>%
  group_by(Method, CBT_ncomm, variable) %>%
  summarise(mean = mean(value))%>%ungroup %>%
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

ggsave("Burkina Faso Analysis/CORR_hybrid.pdf", width = 8, height = 5)




#### --- ERROR RATE PLOTS - COMMUNITY VARIABILITY ----------------------------------

#PMT comparison needs to be predicting

plot_data <- all_results_comm %>%  mutate(IER = 1-Precision,
                                          EER = 1-Sensitivity) %>%
  melt(id.var = c("Method", "CBT_ncomm", "community")) %>%
  mutate(Method = factor(Method, levels = c("Hybrid Score (corrected)","Hybrid Score","CBT Score", "CBT Score (corrected)","CBT DU", "CBT Logit", "PMT OLS"),
                         labels = c(        "Hybrid-AI-EC",            "Hybrid-AI",    "Hybrid",   "Hybrid-EC",           "Hybrid-DU", "Probit", "PMT"))) 


plot_data %>% subset(variable == "EER" & CBT_ncomm %in% c(5, 25) & Method %in% c("Hybrid", "Probit", "PMT" )) %>%
  ggplot() +
  geom_boxplot(aes(x = Method,y = value,  colour = as.factor(CBT_ncomm)),draw_quantiles = c(.25, .5, .75)) +
  scale_colour_grey("Number of\nRanking\nCommunities") +
  theme_bw() +
  labs(x = "Method", y = "Community-level Error Rate")

ggsave(paste0("Burkina Faso Analysis/ER_community_level",PR*100,".pdf"), width = 8, height = 5)



plot_data %>% subset(variable == "EER" & CBT_ncomm %in% c(5, 25) ) %>%
  group_by(CBT_ncomm,Method) %>% 
  summarise(mean = mean(value, na.rm=TRUE),
            var = var(value, na.rm=TRUE)) %>%
  write.csv(paste0("Burkina Faso Analysis/ER_community_level",PR*100,".pdf"))


plot_data %>% subset(variable == "EER" & Method %in% c("Hybrid", "PMT")) %>%
  group_by(CBT_ncomm,Method) %>% 
  summarise(mean = mean(value, na.rm=TRUE),
            var = var(value, na.rm=TRUE)) %>%
  ggplot() +
  geom_point(aes(x = CBT_ncomm, y = var))+
  geom_line(aes(x = CBT_ncomm, y = var, linetype = Method)) +
  theme_bw()  +
  labs(x = "Number of Ranking Communities", y = "Variance of Community-Level Error Rates")
ggsave(paste0("Indonesia Analysis/ER_commlevel",PR*100,".pdf"), width = 8, height = 5)


#### --- COEFFICIENT PLOTS ----------------------------------

variable_labels <- read.csv("Data/Burkina Faso/Cleaning/variables.csv")
#variable_labels_add <- data.frame(Name = "connected", Definition = "Elite connection")

#variable_labels <- rbind(variable_labels, variable_labels_add)

score_order <- all_coef %>% merge(variable_labels, by.x = "parameter", by.y = "Name") %>% 
  dplyr::select(Definition, Order, CB_beta_rank_mean) %>%
  subset(CB_beta_rank_mean != 0)%>% #remove elite connection 0
  melt(id.vars = c("Definition", "Order")) %>%
  group_by(Definition,variable, Order) %>%
  summarise(mean = mean(value)) %>%
  group_by(variable, Order) %>%
  mutate(std_mean =mean/mean(abs(mean)),
         mean_abs = mean(abs(mean))) %>%
  arrange(-Order) 

all_coef %>%merge(variable_labels, by.x = "parameter", by.y = "Name") %>% 
  dplyr::select(Definition, Order, CB_beta_rank_mean, PMT_beta) %>%
  subset(CB_beta_rank_mean != 0)%>% #remove elite connection 0
  melt(id.vars = c("Definition", "Order")) %>%
  group_by(Definition, Order, variable) %>%
  summarise(mean = mean(value)) %>% ungroup() %>%
  group_by(variable) %>%
  mutate(std_mean = mean/mean(abs(mean))) %>%
  ungroup() %>%
  mutate(Definition = factor(Definition, levels = score_order$Definition),
         variable = factor(variable, levels = c("CB_beta_rank_mean", "PMT_beta"),
                           labels = c("Hybrid", "PMT")))%>%
  ggplot() + 
  geom_col(aes(x = Definition, y = std_mean, fill = variable ),position =  position_dodge(width = 0.5)) +
  #facet_grid(Category~., scales = "free_y")+
  coord_flip() + 
  theme_bw() +
  labs(x = "", y = "Standardized Coefficient Estimate") +
  scale_fill_grey("Method") +
  theme(legend.position = c(.9,.9), 
        legend.box.background = element_rect(colour = "black"))


ggsave("Burkina Faso Analysis/coef_score_hillebrecht.pdf", width = 12, height = 12)

all_coef %>%merge(variable_labels, by.x = "parameter", by.y = "Name") %>%
  dplyr::select(Definition, Order, CB_beta_rank_mean, CB_beta_rank_mean_noelite) %>%
  subset(CB_beta_rank_mean != 0)%>% #remove elite connection 0
  melt(id.vars = c("Definition", "Order")) %>%
  group_by(Definition, Order, variable) %>%
  summarise(mean = mean(value)) %>% ungroup() %>%
  group_by(variable) %>%
  mutate(std_mean = mean/mean(abs(mean))) %>%
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

ggsave("Burkina Faso Analysis/coef_score_EC_hillebrecht.pdf", width = 12, height = 12)


