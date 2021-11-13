
all_results <- read.csv("Hillebrecht Analysis/all_results.csv")

all_coef <- read.csv("Hillebrecht Analysis/all_coef.csv")

all_results %>%melt(id.var = c("Method", "CBT_ncomm")) %>%
  mutate(Method = factor(Method, levels = c("Hybrid Score (corrected)","Hybrid Score","CBT Score", "CBT Score (corrected)", "CBT Logit", "PMT OLS"),
                         labels = c("Hybrid Score (corrected)","Hybrid Score","CBT Score","CBT Score (corrected)", "CBT Probit", "PMT OLS"))) %>%
  group_by(Method, CBT_ncomm, variable) %>%
  mutate(mean = median(value ))%>%ungroup%>%
  subset(variable %in% c("TD", "Sensitivity", "Specificity") )%>%
  #subset(Method != "CBT Probit" & Method != "PMT OLS")%>%
  ggplot() +#geom_boxplot(aes(x = Method, y = value,linetype = Method, group = interaction(Method, CBT_prop))) + 
  geom_boxplot(aes(x = Method, y = value, colour = Method, group = interaction(CBT_ncomm, Method))) + 
  stat_summary(aes(x = Method, y = value, colour = Method, group = interaction(CBT_ncomm, Method)),
               fun=mean, geom="point", color="black")+
  #geom_line(aes(x = CBT_prop, y = mean, group = interaction(Method), linetype = Method, colour = Method)) + 
  facet_grid(variable~CBT_ncomm, scales = "free")+ theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = .9))  +
  scale_colour_brewer(type = "qual", palette = "Dark2")# +
#scale_y_log10()
ggsave("Alatas Analysis/all_results.pdf")
