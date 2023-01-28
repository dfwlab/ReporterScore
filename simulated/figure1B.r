######------FIGURE 1B)------######
##load packages
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(tidyr)) install.packages("tidyr")
if(!require(dplyr)) install.packages("dplyr")
if(!require(ggpubr)) install.packages("ggpubr")
if(!require(ggpattern)) install.packages("ggpattern")
if(!require(ggh4x)) install.packages("ggh4x")
if(!require(ggdist)) install.packages("ggdist")

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
pathway_level = (zscore_Pathway %>%
  mutate(pathway = paste0("pathway ", 1:10)) %>%
  arrange(adjusted_rs))$pathway ##pathway level for figure

####plot figure 1B####
##figure 1B-1: reporter scores
fig1B_1 = zscore_Pathway %>%
  mutate(pathway = paste0("pathway ", 1:10)) %>%
  arrange(adjusted_rs) %>%
  mutate(pathway = factor(pathway, levels = pathway)) %>%
  # mutate(pathway = factor(pathway, levels = paste0("pathway ", 10:1))) %>%
  ggplot() +
  geom_tile(aes(x = "", y = pathway, fill = adjusted_rs), color = "white", width = 0.3, height = 1) +
  # scale_fill_viridis_c(option = "A", name="Reporter\nScore") +
  scale_fill_gradient2(n.breaks = 6,
                       low = "#003399", high = "#cc0000", mid = "white",
                       midpoint = 0, limit = c(-4.2,4.2), space = "Lab",
                       name="Reporter\nScore") +
  scale_y_discrete(position = "left") +
  xlab("reporter\nscore") +
  ylab(NULL) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(0,0,0,0, "cm"),
        panel.border = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10, 
                                   face = c("bold", rep("plain", 9)),
                                   margin = margin(0,-0.5,0,0, "cm")),
        legend.title = element_text(size = 10),
        legend.position = "left"
        )

##figure 1B-2: pvalue distribution
fig1B_2_data = data2 %>% select(pathway, pvalues, bin = p_cat) %>%
  mutate(pathway = paste0("pathway ", pathway)) %>%
  mutate(pathway = factor(pathway, levels = pathway_level)) %>%
  mutate(bin = factor(bin, levels = c("<0.05", "[0.05, 0.1)", "[0.1, 0.2)", "[0.2, 0.3)",
                                      "[0.3, 0.4)", "[0.4, 0.5]", ">0.5"))) %>%
  mutate(direct = ifelse(pathway %in% paste0("pathway ", c(2, 9, 3, 10, 4)), "up", "down")) ##all KOs in pathway2, 9, 3, 10, 4 are up-regulated

# fig1B_2 = ggplot(fig1B_2_data) +
#   geom_bar_pattern(aes(bin, pattern = direct, fill = direct), width = 0.5,
#                    color = "black",
#                    pattern_density = 0.08,
#                    pattern_spacing = 0.08,
#                    pattern_size = 0.5,
#                    pattern_fill = "black",
#                    alpha = 0.7)+
#   scale_y_continuous(breaks = c(0, 10, 20))+
#   scale_pattern_discrete(choices = c('stripe', 'circle'))+
#   scale_fill_manual(values = c("up" = "#f46f20", "down" = "#156077"))+
#   facet_grid(rows = "pathway")+
#   theme_bw()+
#   theme(axis.ticks.x = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.title = element_blank(),
#         legend.position = "none",
#         strip.text = element_blank(),
#         axis.text = element_text(family = "Arial", size = 10, color = "black")
#   )

# fig1B_2 = ggplot(fig1B_2_data, aes(x = -log10(pvalues), y = pathway))+
#   ggdist::stat_halfeye(aes(fill = direct), alpha = 0.7, adjust = .5, width = .6, .width = 0, justification = 0, point_colour = NA)+
#   geom_point(## draw horizontal lines instead of points
#              shape = "|", size = 2, alpha = .4) +
#   scale_fill_manual(values = c("up" = "#f46f20", "down" = "#156077"))+
#   theme_bw()+
#   theme(axis.ticks = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.grid.major.y = element_blank(),
#         panel.border = element_blank(),
#         axis.title = element_blank(),
#         legend.position = "none",
#         strip.text = element_blank(),
#         axis.text = element_text(size = 10, color = "black")
#   )
# ggsave("Figure/temp.pdf", width = 4, height = 6)

fig1B_2 = ggplot(fig1B_2_data, aes(x = -log10(pvalues), y = pathway))+
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "#cc340c", size = 0.4)+
  ggdist::stat_halfeye(aes(fill = direct), alpha = 0.8, adjust = .5, width = .6, .width = 0, justification = 0, point_colour = NA)+
  ggbeeswarm::geom_beeswarm(alpha=1, cex = 1, priority = "descending", groupOnX = FALSE, size = 0.5, color = "black")+
  scale_fill_manual(values = c("up" = "#f46f20", "down" = "#156077"))+
  theme_bw()+
  theme(axis.ticks = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        strip.text = element_blank(),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_blank()
  )

fig1B = ggarrange(fig1B_1, fig1B_2, ncol = 2, widths = c(2,4))

ggsave("Figure/figure1B.pdf", bg = "white", width = 7, height = 6)

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
print(p_values) 

#fig1C
# fig1C = zscore_Pathway %>%
#   mutate(pathway = paste0("pathway ", 1:10)) %>%
#   arrange(adjusted_rs) %>%
#   mutate(pathway = factor(pathway, levels = pathway)) %>%
#   mutate(colors = ifelse(adjusted_rs>=0, "P", "N")) %>%
#   mutate(y = ifelse(adjusted_rs>=0, -1.6, 1.6)) %>%
#   mutate(y = ifelse(pathway=="pathway 10", 1.8, y)) %>%
#   ggplot() +
#   geom_hline(yintercept = -1.6, size=0.3, linetype = "dashed", color = "#999999")+
#   geom_hline(yintercept = 1.6, size=0.3, linetype = "dashed", color = "#999999")+
#   geom_bar(aes(x = pathway, y = adjusted_rs, fill = colors), width = 0.5, 
#            stat = "identity", show.legend = FALSE) +
#   scale_fill_manual(values = c("P" = "#cc340c", "N" = "#3f60aa")) +
#   scale_y_continuous(limits = c(-4.2, 4), breaks = seq(-4,4), 
#                      minor_breaks = seq(-4.6,4.4,0.2),
#                      guide = "axis_minor")+
#   geom_hline(yintercept = 0, size = 0.3)+
#   geom_text(aes(x = pathway, y = y, label = pathway), 
#             family = "Arial", size = 3.5) +
#   coord_flip() +
#   xlab(NULL) +
#   ylab("reporter score") +
#   theme_bw() +
#   theme(axis.ticks.y = element_blank(),
#         ggh4x.axis.ticks.length.minor = rel(0.5),
#         panel.grid = element_blank(),
#         panel.background = element_rect(size = 0.5, color = "black"),
#         plot.margin = margin(0,0,0,0, "cm"),
#         panel.border = element_blank(),
#         axis.text.x = element_text(family = "Arial", size = 10),
#         axis.text.y = element_blank(),
#         axis.title.x = element_text(family = "Arial", size = 10, face = "bold")
#   )
# fig1C
# ggsave("Figure/figure1C.png", bg = "white", width = 2, height = 6)

