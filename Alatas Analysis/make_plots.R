rm(list = ls())
library(dplyr)
library(ggplot2)
library(reshape2)
library(lmomco)
all_results <- read.csv("Alatas Analysis/all_results.csv")

pmt_corrected <- read.csv("Alatas Analysis/PMT_nonconverge_corrected_cbt.csv") %>%
  mutate(variable = as.factor("EER"), Method = as.factor("PMT")) %>%
  rename(mean=mean_EER ) %>%
  relocate(Method, CBT_ncomm, variable, mean)

all_coef <- read.csv("Alatas Analysis/coef_total_sample.csv")

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


#### --- ERROR RATE PLOTS ----------------------------------

plot_data <- all_results %>%  mutate(IER = 1-Precision,
                        EER = 1-Sensitivity) %>%
  melt(id.var = c("Method", "CBT_ncomm")) %>%
  mutate(Method = factor(Method, levels = c("Hybrid Score (corrected)","Hybrid Score","CBT Score", "CBT Score (corrected)","CBT DU", "CBT Logit", "PMT OLS"),
                         labels = c("Hybrid-AI-EC","Hybrid-AI","Hybrid","Hybrid-EC","Hybrid-DU", "Probit", "PMT"))) %>%
  group_by(Method, CBT_ncomm, variable) %>%
  summarise(mean = mean(value ))%>%ungroup %>%
  subset(variable %in% c( "EER"))

plot_data[plot_data$Method == "PMT",] <- pmt_corrected


plot_data %>%
  subset( Method %in% c("Hybrid", "PMT", "Probit")  )%>%
  ggplot() + geom_line(aes(x = CBT_ncomm, y = mean, linetype = Method)) +
  geom_point(aes(x = CBT_ncomm, y = mean)) +
  #geom_linerange(aes(x = CBT_ncomm, ymin = min,ymax=max, linetype = Method))+
  theme_bw() +
  labs(x = "Number of Ranking Communities", y = "Average Error Rate")+ 
  theme(legend.position = c(0.9, 0.8)) +
  theme(legend.box.background = element_rect(colour = "black"))

ggsave("Alatas Analysis/ER_hybrid.pdf", width = 8, height = 5)


plot_data %>%
  subset( Method %in% c("Hybrid", "Hybrid-AI")  )%>%
  ggplot() + geom_line(aes(x = CBT_ncomm, y = mean, linetype = Method)) +
  geom_point(aes(x = CBT_ncomm, y = mean)) +
  #geom_linerange(aes(x = CBT_ncomm, ymin = min,ymax=max, linetype = Method))+
  theme_bw() +
  labs(x = "Number of Ranking Communities", y = "Average Error Rate")+ 
  theme(legend.position = c(0.9, 0.8))+
  theme(legend.box.background = element_rect(colour = "black"))

ggsave("Alatas Analysis/ER_hybrid_AI.pdf", width = 8, height = 5)

plot_data %>%
  subset( Method %in% c("Hybrid", "Hybrid-EC")  )%>%
  ggplot() + geom_line(aes(x = CBT_ncomm, y = mean, linetype = Method)) +
  geom_point(aes(x = CBT_ncomm, y = mean)) +
  #geom_linerange(aes(x = CBT_ncomm, ymin = min,ymax=max, linetype = Method))+
  theme_bw() +
  labs(x = "Number of Ranking Communities", y = "Average Error Rate")+ 
  theme(legend.position = c(0.9, 0.8))+
  theme(legend.box.background = element_rect(colour = "black"))

ggsave("Alatas Analysis/ER_hybrid_EC.pdf", width = 8, height = 5)

#(no DU for alatas)

#### --- COEFFICIENT PLOTS ----------------------------------

variable_labels <- read.csv("Data/Indonesia/Cleaning/variables.csv")

#variable_labels_add <- data.frame(Name = "connected", Definition = "Elite connection")

#variable_labels <- rbind(variable_labels, variable_labels_add)

score_order <- all_coef %>% merge(variable_labels, by.x = "parameter", by.y = "Name") %>% 
  dplyr::select(Definition,Category, CB_beta_rank_mean) %>%
  subset(CB_beta_rank_mean != 0)%>% #remove elite connection 0
  melt(id.vars = c("Definition", "Category")) %>%
  mutate(par_est = ifelse(abs(value) < 0.01, 0, value)) %>%
  group_by(Category,variable) %>%
  mutate(std_mean =par_est/harmonic.mean(par_est)$harmean) %>%
  arrange(Category,std_mean) 


all_coef %>%merge(variable_labels, by.x = "parameter", by.y = "Name") %>% 
  dplyr::select(Definition, Category, CB_beta_rank_mean, PMT_beta) %>%
  subset(CB_beta_rank_mean != 0)%>% #remove elite connection 0
  melt(id.vars = c("Definition", "Category")) %>%
  mutate(par_est = ifelse(abs(value) < 0.01, 0, value)) %>%
  #group_by(Definition, Category, variable) %>%
  #summarise(par_est = mean(value)) %>% ungroup() %>%
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

ggsave("Alatas Analysis/coef_score_alatas.pdf", width = 12, height = 12)

all_coef %>%merge(variable_labels, by.x = "parameter", by.y = "Name") %>%
  dplyr::select(Definition, Category, CB_beta_rank_mean, CB_beta_rank_mean_noelite) %>%
  subset(CB_beta_rank_mean != 0)%>% #remove elite connection 0
  melt(id.vars = c("Definition", "Category")) %>%
  #group_by(Definition, Category, variable) %>%
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

ggsave("Alatas Analysis/coef_score_EC_alatas.pdf", width = 12, height = 12)



