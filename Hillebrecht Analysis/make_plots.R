library(dplyr)
library(ggplot2)
library(reshape2)
all_results <- read.csv("Hillebrecht Analysis/all_results.csv")

all_coef <- read.csv("Hillebrecht Analysis/all_coef.csv")

all_results %>%  mutate(IER = 1-Precision,
                        EER = 1-Sensitivity) %>%
  melt(id.var = c("Method", "CBT_ncomm")) %>%
  mutate(Method = factor(Method, levels = c("Hybrid Score","CBT Score", "CBT Probit", "PMT OLS"),
                         labels = c("Hybrid Score","CBT Score", "CBT Probit", "PMT OLS"))) %>%
  group_by(Method, CBT_ncomm, variable) %>%
  mutate(mean = median(value ))%>%ungroup%>%
  subset(variable %in% c("TD", "IER", "EER")  )%>%#& !Method %in% c("CBT Probit", "PMT OLS")
  ggplot() + 
  geom_boxplot(aes(x = Method, y = value, colour = Method, group = interaction(CBT_ncomm, Method))) + 
  stat_summary(aes(x = Method, y = value, colour = Method, group = interaction(CBT_ncomm, Method)),
               fun=mean, geom="point", color="black")+
  facet_grid(variable~CBT_ncomm, scales = "free")+ theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = .9))  +
  scale_colour_brewer(type = "qual", palette = "Dark2")
ggsave("Hillebrecht Analysis/all_results.pdf")
