######------FIGURE 1C)------######
##load packages
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(tidyr)) install.packages("tidyr")
if(!require(dplyr)) install.packages("dplyr")
if(!require(ggpubr)) install.packages("ggpubr")
if(!require(ggpattern)) install.packages("ggpattern")
if(!require(ggh4x)) install.packages("ggh4x")
if(!require(ggrepel)) install.packages("ggrepel")
library(scales)

##set working directory
setwd("~/OneDrive/Documents/Postdoc/2022/10. reporter score/manuscript/mSystems/")
##get simulated dataset
data = read.csv("Figure/data1.csv")

##compute zscore
data2 = data %>%
  group_by(pathway) %>%
  mutate(zscore_KO = qnorm(1-pvalues),
         KO_id = KO_id+1) %>%
  ungroup()
##compute reporter score for each pathway
zscore_Pathway = data2 %>%
  group_by(pathway) %>%
  summarise(n = n(),
            zscore_Pathway = sum(zscore_KO)/sqrt(n)) %>%
  ungroup()
##get the KO list
zscore_KO_all = data2 %>%
  select(KO_id, zscore_KO) %>%
  distinct() %>%
  arrange(KO_id)

##get random mean and sd; and compute adjusted reporter score
adjusted_rs = c()
for (i in 1:10) {
  reporter_score_random = c()
  set.seed(123)
  for (j in 1:10000) {
    zscore_KO_random = sample(zscore_KO_all$zscore_KO, zscore_Pathway$n[i], replace = FALSE)
    zscore_Pathway_random = sum(zscore_KO_random)/sqrt(zscore_Pathway$n[i])
    reporter_score_random = append(reporter_score_random, zscore_Pathway_random)
  }
  print(i)
  print(mean(reporter_score_random))
  print(sd(reporter_score_random))
  adjusted_rs = append(adjusted_rs, (zscore_Pathway$zscore_Pathway[i]-mean(reporter_score_random))/sd(reporter_score_random))
}
zscore_Pathway$adjusted_rs = adjusted_rs

##hypergeometric test
nsig = length(unique(data2[data2$pvalues<0.05,]$KO_id))
nonsig = 200-nsig
p_values = c()
for (i in c(1:10)) {
  a = nrow(filter(data2, pathway==i & pvalues<0.05))
  b = nrow(filter(data2, pathway==i & pvalues>=0.05))
  c = 50-a
  d = nonsig - b
  p_value = fisher.test(matrix(c(a,b,c,d), 2, 2, byrow = T), alternative = "greater") ##corresponds to a one-sided version of Fisher’s exact test
  p_values = append(p_values, p_value$p.value)
}

####plot figure 1C####
data4fig1C = zscore_Pathway %>%
  mutate(pathway = paste("pathway", pathway),
         hyperG = p_values) 
fig1C = ggplot(data4fig1C, aes(adjusted_rs, -log10(hyperG))) +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "#999999", size = 0.4, alpha = 1)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#cc340c", size = 0.5, alpha = 0.5)+
  geom_vline(xintercept = 1.64, linetype = "dashed", color = "#cc340c", size = 0.4, alpha = 1)+
  geom_point() +
  geom_text_repel(aes(label = pathway), size = 3.5, direction = "both") +
  scale_x_continuous(limits = c(min(data4fig1C$adjusted_rs), 4), breaks = c(-4, -1.64, 0, 1.64, 4))+
  xlab("reporter score") +
  ylab(~atop(paste(bold("-log10 "), bolditalic("P"), bold(" value")))) +
  
  theme_classic()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 12, color = "black", face = "bold"),
        axis.title.y = element_text(size = 12, color = "black",
                                    margin = margin(0,-0.5,0,0, "cm"), face = "bold"))
fig1C
ggsave("Figure/figure1C.pdf", width = 5, height=5)
