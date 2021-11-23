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

#### --- ERROR RATE PLOTS ----------------------------------

plot_data <- all_results %>%  mutate(IER = 1-Precision,
                        EER = 1-Sensitivity) %>%
  melt(id.var = c("Method", "CBT_ncomm")) %>%
  mutate(Method = factor(Method, levels = c("Hybrid Score (corrected)","Hybrid Score","CBT Score", "CBT Score (corrected)", "CBT Logit", "PMT OLS"),
                         labels = c("Hybrid Score (corrected)","Hybrid Score","CBT Score","CBT Score (corrected)", "CBT Probit", "PMT OLS"))) %>%
  group_by(Method, CBT_ncomm, variable) %>%
  mutate(mean = median(value ),
         min = min(value),
         max = max(value))%>%ungroup

plot_data%>%
  subset(variable %in% c( "IER") & Method %in% c("CBT Score", "PMT OLS", "CBT Probit")  )%>%#& !Method %in% c("CBT Probit", "PMT OLS")
  ggplot() +
  geom_boxplot(aes(x = as.factor(CBT_ncomm), y = value, group = interaction(CBT_ncomm, Method), colour = Method)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = .9))  +
  scale_colour_grey()+theme_bw()+
  labs(x = "Number of Ranking Communities", y = "Average Error Rate")+ 
  theme(legend.position = c(0.9, 0.8))

ggsave("Alatas Analysis/ER_score_v1.pdf")

plot_data %>%
  subset(variable %in% c( "IER") & Method %in% c("CBT Score", "PMT OLS", "CBT Probit")  )%>%
  ggplot() + geom_line(aes(x = CBT_ncomm, y = mean, linetype = Method)) +
  geom_point(aes(x = CBT_ncomm, y = mean)) +
  #geom_linerange(aes(x = CBT_ncomm, ymin = min,ymax=max, linetype = Method))+
  theme_bw() +
  labs(x = "Number of Ranking Communities", y = "Average Error Rate")+ 
  theme(legend.position = c(0.9, 0.8))

ggsave("Alatas Analysis/ER_score_v2.pdf")


#### --- COEFFICIENT PLOTS ----------------------------------

score_order <- all_coef %>% subset(CBT_ncomm == 200) %>%
  dplyr::select(parameter, CB_beta_rank_mean) %>%
  melt(id.vars = c("parameter")) %>%
  group_by(parameter, variable) %>%
  summarise(mean = mean(value)) %>%
  arrange(mean) 

all_coef %>% subset(CBT_ncomm == 200) %>%
  dplyr::select(parameter, CB_beta_rank_mean, PMT_beta) %>%
  melt(id.vars = c("parameter")) %>%
  group_by(parameter, variable) %>%
  summarise(mean = mean(value)) %>% ungroup() %>%
  mutate(parameter = factor(parameter, levels = score_order$parameter),
         variable = factor(variable, levels = c("CB_beta_rank_mean", "PMT_beta"),
                           labels = c("CBT Score Model", "PMT")))%>%
  ggplot() + 
  geom_line(aes(x = parameter, y = mean, linetype = variable, group = variable )) +
  coord_flip() + 
  theme_bw() +
  labs(x = "", y = "Average Coefficient Estimate") +
  scale_linetype("Method")

