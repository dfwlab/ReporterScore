######------FIGURE 2------######
##load packages
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(tidyr)) install.packages("tidyr")
if(!require(dplyr)) install.packages("dplyr")
if(!require(ggbeeswarm)) install.packages("ggbeeswarm")
if(!require(ggpubr)) install.packages("ggpubr")
if(!require(ggrepel)) install.packages("ggrepel")
##set working directory
setwd("~/OneDrive/Documents/Postdoc/2022/10. reporter score/manuscript/iMeta/20230122_revised version-DF-LL/")

##read data
fix_names <- function(x) gsub("\\s+", ".", x)
data = readxl::read_xlsx("supplementary material-DF.xlsx", sheet = 2, .name_repair = fix_names)
##data processing
colnames(data)[2] = "Impact.factors"
colnames(data)[14] = "2023"
data4fig2A = data[,c(1:5)] %>%
  mutate(Year = as.character(Year)) %>%
  # add_row(Year = "Jan 2023") %>%
  mutate(Year = factor(Year),
         Color = if_else(Year == 2005, "O", if_else(Impact.factors>20, "Y", "N"))) ##set circle filling color
label = filter(data4fig2A, Color %in% c("O", "Y"))

CitationByYear = colSums(data[-1,c(6:14)])
data4fig2B = data.frame(Year = c(2005, 2015:2022, "Jan 2023"), nCitations = c(0, CitationByYear))

##fig2A
seg <- tibble(x = c(0),
              xend = c(10.6),
              y = c(0),
              yend = c(0))
fig2A = ggplot() +
  geom_beeswarm(data = data4fig2A,
                aes(Year, Impact.factors, size = Times.cited, fill = Color, color = Color),
                alpha=0.7, cex = 2, priority = "descending",
                shape = 21, stroke = 0.3) +
  geom_segment(data = seg[1,c(1:4)], 
               aes(x = x, xend = xend,
                   y = y, yend = yend), 
               color = "black", linetype = "solid", lineend = "round", size=0.2)+
  geom_text_repel(data = label, aes(factor(Year), Impact.factors, label = Paper),
                  size = 4.2, direction = "y", nudge_y = c(-2.3, 2.4,-2.3))+
  ylab("5-Year Journal Impact Factor")+
  scale_y_continuous(limits = c(0, 70))+
  scale_fill_manual(values = c("Y" = "#cc340c", "N" = "#999999", "O" = "#ffbc14"))+
  scale_color_manual(values = c("Y" = "white", "N" = "white", "O" = "white"))+
  scale_size("Times cited",limits=c(0,1600),breaks=c(0, 50, 100, 200, 500, 1000, 1500), range=c(2,6))+
  theme_void()+
  guides(size = guide_legend(nrow = 1),
         fill = "none", color = "none") +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", size = 0.4),
        plot.margin = margin(0,0.1,0.1,0.1, "cm"),
        legend.position = "top",
        axis.text.y = element_text(size = 13, color = "black", hjust = 1),
        axis.text.x = element_text(size = 13, color = "black", face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, color = "black", angle = 90,
                                    margin = margin(0,0.33,0,0, "cm"), face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12)
  )
fig2A

####fig2B
seg2 <- tibble(x = c(0),
              xend = c(10.6),
              y = c(0),
              yend = c(0))

fig2B = ggplot(data4fig2B) +
  geom_bar(aes(factor(Year), cumsum(nCitations)), stat = "identity", width = 0.4, fill="#003366", alpha=0.75)+
  geom_segment(data = seg[1,c(1:4)],
               aes(x = x, xend = xend,
                   y = y, yend = yend),
               color = "white", linetype = "solid", lineend = "round", size=0.2)+
  geom_text(data = data4fig2B[c(2:10),],
            aes(factor(Year), cumsum(nCitations), label = cumsum(nCitations)), 
            position = position_dodge2(width = 0.8, preserve = 'single'), 
            vjust = 1.5, hjust = 0.5, size = 4.3, color = "black")+
  scale_y_reverse(limits=c(5500, 0), breaks=c(5500, 4125, 2750, 1375, 0), 
                  labels=c("", "", "", 100, 0))+
  ylab("Cumulative Number of Citations")+
  xlab("")+
  theme_void()+
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        # panel.grid.major.y = element_line(color = "grey90", size = 0.4),
        plot.margin = margin(0,0.1,0.1,0.1, "cm"),
        legend.position = "top",
        # axis.line.y = element_line(color = "black", size = 0.4),
        axis.text.y = element_text(size = 12, color = "white", hjust = 1),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 12, color = "black", 
                                    margin = margin(0.15,0,0,0, "cm"), face = "bold"),
        axis.title.y = element_text(size = 12, color = "black", angle = 90,
                                    margin = margin(0,0.15,0,0, "cm"), face = "bold"),
        legend.title = element_text(size = 8, face = "bold"),
        legend.text = element_text(size = 8)
  )

##arrange two subfigures
fig2 = ggarrange(fig2A, fig2B, nrow = 2, heights = c(2,1))
##save figure file
ggsave("figure2.pdf", bg = "white", width = 7, height=7)
