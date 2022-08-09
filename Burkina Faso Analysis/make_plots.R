library(dplyr)
library(ggplot2)
library(reshape2)

all_results <- read.csv("Burkina Faso Analysis/all_results.csv")

#pmt_corrected <- read.csv("Burkina Faso Analysis/PMT_nonconverge_corrected_cbt.csv") %>%
#  mutate(variable = as.factor("EER"), Method = as.factor("PMT")) %>%
#  rename(mean=mean_EER ) %>%
#  relocate(Method, CBT_ncomm, variable, mean)

all_coef <- read.csv("Burkina Faso Analysis/coef_total_sample.csv")



%>%
  mutate(hybrid_noelite_inclusion = hybrid_noelite_rank <= treat_rate,
         hybrid_inclusion = hybrid_rank <= treat_rate,
         pmt_inclusion = pmt_rank <= treat_rate,
         consumption_inclusion = consumption_rank<=treat_rate,
         cbt_model_inclusion = cbt_model_rank<=treat_rate,
         cbt_DU_model_inclusion = cbt_DU_model_rank<=treat_rate,
         cbt_model_noelite_inclusion = cbt_model_rank_noelite<=treat_rate,
         CBT_LR_inclusion = CBT_LR_rank<=treat_rate,
         cbt_inclusion = ifelse(treated == 1, TRUE, FALSE)) %>%ungroup() %>%
  mutate_at(vars(matches("inclusion")), as.factor)


rbind(
  confusionMatrix(Program_data$hybrid_noelite_inclusion,   Program_data$cbt_inclusion,positive = "TRUE")$byClass,
  confusionMatrix(Program_data$hybrid_inclusion,           Program_data$cbt_inclusion,positive = "TRUE")$byClass,
  confusionMatrix(Program_data$cbt_model_noelite_inclusion,Program_data$cbt_inclusion,positive = "TRUE")$byClass,
  confusionMatrix(Program_data$cbt_model_inclusion,        Program_data$cbt_inclusion,positive = "TRUE")$byClass,
  confusionMatrix(Program_data$cbt_DU_model_inclusion,    Program_data$cbt_inclusion,positive = "TRUE")$byClass,
  confusionMatrix(Program_data$pmt_inclusion,              Program_data$cbt_inclusion,positive = "TRUE")$byClass,
  confusionMatrix(Program_data$CBT_LR_inclusion,           Program_data$cbt_inclusion,positive = "TRUE")$byClass) %>%as.data.frame%>%
  mutate(Method = c( "Hybrid Score (corrected)","Hybrid Score","CBT Score (corrected)", "CBT Score","CBT DU", "PMT OLS", "CBT Logit"),
         CBT_ncomm = CBT_ncomm,
         TD = Sensitivity - (1-Specificity),
         rep = rep)


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

ggsave("Burkina Faso Analysis/ER_hybrid.pdf", width = 8, height = 5)


plot_data %>%
  subset( Method %in% c("Hybrid", "Hybrid-AI")  )%>%
  ggplot() + geom_line(aes(x = CBT_ncomm, y = mean, linetype = Method)) +
  geom_point(aes(x = CBT_ncomm, y = mean)) +
  #geom_linerange(aes(x = CBT_ncomm, ymin = min,ymax=max, linetype = Method))+
  theme_bw() +
  labs(x = "Number of Ranking Communities", y = "Average Error Rate")+ 
  theme(legend.position = c(0.9, 0.8))+
  theme(legend.box.background = element_rect(colour = "black"))

ggsave("Burkina Faso Analysis/ER_hybrid_AI.pdf", width = 8, height = 5)

plot_data %>%
  subset( Method %in% c("Hybrid", "Hybrid-DU")  )%>%
  ggplot() + geom_line(aes(x = CBT_ncomm, y = mean, linetype = Method)) +
  geom_point(aes(x = CBT_ncomm, y = mean)) +
  #geom_linerange(aes(x = CBT_ncomm, ymin = min,ymax=max, linetype = Method))+
  theme_bw() +
  labs(x = "Number of Ranking Communities", y = "Average Error Rate")+ 
  theme(legend.position = c(0.9, 0.8))+
  theme(legend.box.background = element_rect(colour = "black"))

ggsave("Burkina Faso Analysis/ER_hybrid_DU.pdf", width = 8, height = 5)



all_results %>%  mutate(IER = 1-Precision,
                        EER = 1-Sensitivity) %>%
  melt(id.var = c("Method", "CBT_ncomm")) %>%
  mutate(Method = factor(Method, levels = c("Hybrid Score (corrected)","Hybrid Score","CBT Score", "CBT Score (corrected)","CBT DU", "CBT Logit", "PMT OLS"),
                         labels = c(        "Hybrid-AI-EC",            "Hybrid-AI",    "Hybrid",   "Hybrid-EC",           "Hybrid-DU", "Probit", "PMT"))) %>%
  subset(variable %in% c( "IER") & Method %in% c("Hybrid", "Hybrid-DU")  )%>%
  ggplot() + geom_boxplot(aes(x = as.factor(CBT_ncomm), y = value,colour = Method)) +
  #geom_point(aes(x = CBT_ncomm, y = mean)) +
  #geom_linerange(aes(x = CBT_ncomm, ymin = min,ymax=max, linetype = Method))+
  theme_bw() +
  labs(x = "Number of Ranking Communities", y = "Average Error Rate")+ 
  theme(legend.position = c(0.9, 0.9))

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


