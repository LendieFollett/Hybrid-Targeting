rm(list = ls())
library(dplyr)
library(ggplot2)
library(reshape2)
all_results <- read.csv("Alatas Analysis/all_results.csv")

all_coef <- read.csv("Alatas Analysis/all_coef.csv")

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
  mutate(Method = factor(Method, levels = c("Hybrid Score (corrected)","Hybrid Score","CBT Score", "CBT Score (corrected)", "CBT Logit", "PMT OLS"),
                         labels = c("Hybrid-AI-EC","Hybrid-AI","Hybrid","Hybrid-EC", "Probit", "PMT"))) %>%
  group_by(Method, CBT_ncomm, variable) %>%
  mutate(mean = median(value ),
         min = min(value),
         max = max(value))%>%ungroup


plot_data %>%
  subset(variable %in% c( "IER") & Method %in% c("Hybrid", "PMT", "Probit")  )%>%
  ggplot() + geom_line(aes(x = CBT_ncomm, y = mean, linetype = Method)) +
  geom_point(aes(x = CBT_ncomm, y = mean)) +
  #geom_linerange(aes(x = CBT_ncomm, ymin = min,ymax=max, linetype = Method))+
  theme_bw() +
  labs(x = "Number of Ranking Communities", y = "Average Error Rate")+ 
  theme(legend.position = c(0.9, 0.8))

ggsave("Alatas Analysis/ER_hybrid.pdf", width = 8, height = 4)


plot_data %>%
  subset(variable %in% c( "IER") & Method %in% c("Hybrid", "Hybrid-AI")  )%>%
  ggplot() + geom_line(aes(x = CBT_ncomm, y = mean, linetype = Method)) +
  geom_point(aes(x = CBT_ncomm, y = mean)) +
  #geom_linerange(aes(x = CBT_ncomm, ymin = min,ymax=max, linetype = Method))+
  theme_bw() +
  labs(x = "Number of Ranking Communities", y = "Average Error Rate")+ 
  theme(legend.position = c(0.9, 0.8))

ggsave("Alatas Analysis/ER_hybrid_AI.pdf", width = 8, height = 4)

plot_data %>%
  subset(variable %in% c( "IER") & Method %in% c("Hybrid", "Hybrid-EC")  )%>%
  ggplot() + geom_line(aes(x = CBT_ncomm, y = mean, linetype = Method)) +
  geom_point(aes(x = CBT_ncomm, y = mean)) +
  #geom_linerange(aes(x = CBT_ncomm, ymin = min,ymax=max, linetype = Method))+
  theme_bw() +
  labs(x = "Number of Ranking Communities", y = "Average Error Rate")+ 
  theme(legend.position = c(0.9, 0.8))

ggsave("Alatas Analysis/ER_hybrid_EC.pdf", width = 8, height = 4)


#### --- COEFFICIENT PLOTS ----------------------------------

variable_labels <- read.csv("Data/Indonesia/Cleaning/variables.csv")

variable_labels_add <- data.frame(Name = "connected", Definition = "Elite connection")

variable_labels <- rbind(variable_labels, variable_labels_add)

score_order <- all_coef %>% merge(variable_labels, by.x = "parameter", by.y = "Name") %>% 
  subset(CBT_ncomm == 200) %>%
  dplyr::select(Definition, CB_beta_rank_mean) %>%
  melt(id.vars = c("Definition")) %>%
  group_by(Definition, variable) %>%
  summarise(mean = mean(value)) %>%
  group_by(variable) %>%
  mutate(std_mean = mean/mean[Definition == "Household floor area per capita"]) %>%
  arrange(mean) 

all_coef %>%merge(variable_labels, by.x = "parameter", by.y = "Name") %>% subset(CBT_ncomm == 200) %>%
  dplyr::select(Definition, CB_beta_rank_mean, PMT_beta) %>%
  melt(id.vars = c("Definition")) %>%
  group_by(Definition, variable) %>%
  summarise(mean = mean(value)) %>% ungroup() %>%
  group_by(variable) %>%
  mutate(std_mean = mean/mean[Definition == "Household floor area per capita"]) %>%
  mutate(Definition = factor(Definition, levels = score_order$Definition),
         variable = factor(variable, levels = c("CB_beta_rank_mean", "PMT_beta"),
                           labels = c("Hybrid", "PMT")))%>%
  ggplot() + 
  geom_line(aes(x = Definition, y = std_mean, linetype = variable, group = variable )) +
  coord_flip() + 
  theme_bw() +
  labs(x = "", y = "Standardized Coefficient Estimate \n (Relative to household floor area per capita)") +
  scale_linetype("Method")

ggsave("Hillebrecht Analysis/coef_score.pdf", width = 8, height = 8)





all_coef %>%merge(variable_labels, by.x = "parameter", by.y = "Name") %>% subset(CBT_ncomm == 200) %>%
  dplyr::select(Definition, Hybrid_beta_rank_mean, Hybrid_beta_rank_mean_noelite) %>%
  melt(id.vars = c("Definition")) %>%
  group_by(Definition, variable) %>%
  summarise(mean = mean(value)) %>% ungroup() %>%
  mutate(Definition = factor(Definition, levels = score_order$Definition),
         variable = factor(variable, levels = c("Hybrid_beta_rank_mean", "Hybrid_beta_rank_mean_noelite"),
                           labels = c("Hybrid", "Hybrid_EC")))%>%
  ggplot() + 
  geom_line(aes(x = Definition, y = mean, linetype = variable, group = variable )) +
  coord_flip() + 
  theme_bw() +
  labs(x = "", y = "Average Coefficient Estimate") +
  scale_linetype("Method")




