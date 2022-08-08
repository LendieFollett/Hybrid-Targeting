rm(list = ls())
library(dplyr)
library(ggplot2)
library(reshape2)
library(lmomco)


all_results <- read.csv("Indonesia Analysis/all_results.csv")

#pmt_corrected <- read.csv("Indonesia Analysis/PMT_nonconverge_corrected_cbt.csv") %>%
#  mutate(variable = as.factor("EER"), Method = as.factor("PMT")) %>%
#  rename(mean=mean_EER ) %>%
#  relocate(Method, CBT_ncomm, variable, mean)

all_coef <- read.csv("Indonesia Analysis/coef_total_sample.csv")

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

%>%
  mutate(hybrid_noelite_inclusion = hybrid_noelite_rank <= poverty_rate,
         hybrid_inclusion = hybrid_rank <= poverty_rate,
         pmt_inclusion = pmt_rank <= poverty_rate,
         consumption_inclusion = consumption_rank<=poverty_rate,
         cbt_model_inclusion = cbt_model_rank<=poverty_rate,
         cbt_model_noelite_inclusion = cbt_model_rank_noelite<=poverty_rate,
         CBT_LR_inclusion = CBT_LR_rank<=poverty_rate,
         cbt_inclusion = cbt_rank <= poverty_rate) %>%ungroup() %>%
  mutate_at(vars(matches("inclusion")), as.factor)

#Within community
temp = Program_data%>%group_by(village, province, district, subdistrict, poverty_rate) %>%
  summarise(hybrid_noelite_sens = confusionMatrix(hybrid_noelite_inclusion,   cbt_inclusion,positive = "TRUE")$byClass[1],
            hybrid_noelite_spec = confusionMatrix(hybrid_noelite_inclusion,   cbt_inclusion,positive = "TRUE")$byClass[2],
            hybrid_sens = confusionMatrix(hybrid_inclusion,   cbt_inclusion,positive = "TRUE")$byClass[1],
            hybrid_spec = confusionMatrix(hybrid_inclusion,   cbt_inclusion,positive = "TRUE")$byClass[2],
            cbt_model_noelite_sens = confusionMatrix(cbt_model_noelite_inclusion,   cbt_inclusion,positive = "TRUE")$byClass[1],
            cbt_model_noelite_spec = confusionMatrix(cbt_model_noelite_inclusion,   cbt_inclusion,positive = "TRUE")$byClass[2],
            cbt_model_sens = confusionMatrix(cbt_model_inclusion,   cbt_inclusion,positive = "TRUE")$byClass[1],
            cbt_model_spec = confusionMatrix(cbt_model_inclusion,   cbt_inclusion,positive = "TRUE")$byClass[2],
            pmt_sens = confusionMatrix(pmt_inclusion,   cbt_inclusion,positive = "TRUE")$byClass[1],
            pmt_spec = confusionMatrix(pmt_inclusion,   cbt_inclusion,positive = "TRUE")$byClass[2],
            CBT_LR_sens = confusionMatrix(CBT_LR_inclusion,   cbt_inclusion,positive = "TRUE")$byClass[1],
            CBT_LR_spec = confusionMatrix(CBT_LR_inclusion,   cbt_inclusion,positive = "TRUE")$byClass[2],
  )

#across communities
r[[i]] <- rbind(
  confusionMatrix(Program_data$hybrid_noelite_inclusion,   Program_data$cbt_inclusion,positive = "TRUE")$byClass,
  confusionMatrix(Program_data$hybrid_inclusion,           Program_data$cbt_inclusion,positive = "TRUE")$byClass,
  confusionMatrix(Program_data$cbt_model_noelite_inclusion,Program_data$cbt_inclusion,positive = "TRUE")$byClass,
  confusionMatrix(Program_data$cbt_model_inclusion,        Program_data$cbt_inclusion,positive = "TRUE")$byClass,
  confusionMatrix(Program_data$pmt_inclusion,              Program_data$cbt_inclusion,positive = "TRUE")$byClass,
  confusionMatrix(Program_data$CBT_LR_inclusion,           Program_data$cbt_inclusion,positive = "TRUE")$byClass) %>%as.data.frame%>%
  mutate(Method = c( "Hybrid Score (corrected)","Hybrid Score","CBT Score (corrected)", "CBT Score", "PMT OLS", "CBT Logit"),
         CBT_ncomm = CBT_ncomm,
         TD = Sensitivity - (1-Specificity),
         rep = rep) 


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

ggsave("Indonesia Analysis/ER_hybrid.pdf", width = 8, height = 5)


plot_data %>%
  subset( Method %in% c("Hybrid", "Hybrid-AI")  )%>%
  ggplot() + geom_line(aes(x = CBT_ncomm, y = mean, linetype = Method)) +
  geom_point(aes(x = CBT_ncomm, y = mean)) +
  #geom_linerange(aes(x = CBT_ncomm, ymin = min,ymax=max, linetype = Method))+
  theme_bw() +
  labs(x = "Number of Ranking Communities", y = "Average Error Rate")+ 
  theme(legend.position = c(0.9, 0.8))+
  theme(legend.box.background = element_rect(colour = "black"))

ggsave("Indonesia Analysis/ER_hybrid_AI.pdf", width = 8, height = 5)

plot_data %>%
  subset( Method %in% c("Hybrid", "Hybrid-EC")  )%>%
  ggplot() + geom_line(aes(x = CBT_ncomm, y = mean, linetype = Method)) +
  geom_point(aes(x = CBT_ncomm, y = mean)) +
  #geom_linerange(aes(x = CBT_ncomm, ymin = min,ymax=max, linetype = Method))+
  theme_bw() +
  labs(x = "Number of Ranking Communities", y = "Average Error Rate")+ 
  theme(legend.position = c(0.9, 0.8))+
  theme(legend.box.background = element_rect(colour = "black"))

ggsave("Indonesia Analysis/ER_hybrid_EC.pdf", width = 8, height = 5)

#(no DU for alatas)

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



