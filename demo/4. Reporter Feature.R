######------Demo------######
##load packages
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(tidyr)) install.packages("tidyr")
if(!require(dplyr)) install.packages("dplyr")
if(!require(ggpubr)) install.packages("ggpubr")
if(!require(ggpattern)) install.packages("ggpattern")
if(!require(ggh4x)) install.packages("ggh4x")
if(!require(ggdist)) install.packages("ggdist")
if(!require(ggrepel)) install.packages("ggrepel")
library("piano")

##set working directory
setwd("~/OneDrive/Documents/Postdoc/2022/10. reporter score/manuscript/iMeta/20230122_revised version-DF-LL/demo")

##load data
kos_res = read.csv("KO差异分析.tsv", sep = "\t", header = 1, row.names = 1)
path2ko = read.csv("KEGG_path2ko.csv")
path2name = read.csv("KEGG_pathways.csv", row.names = 1)

##compute reporter feature
gsc = cbind(path2ko$KO, path2ko$Pathway)
gsc = loadGSC(gsc) ##compute gene-set collection 
pvalue = kos_res$P.Value
names(pvalue) = kos_res$ID
direction = ifelse(kos_res$FC<1, -1, 1)
names(direction) = kos_res$ID

set.seed(1)
res = runGSA(pvalue, 
             direction,
             geneSetStat = "reporter",
             signifMethod = "nullDist",
             gsc = gsc,
             nPerm = 500)
# exploreGSAres(res)
GSAsummaryTable(res, save = TRUE, file = "reporterfeature.xls") ##save to a file


##figureS1B: heatmap
reslist <- list(res, res)
nod <- consensusScores(resList = reslist, class = "non", 
                       adjusted = TRUE, method = "median", n = length(reslist[[1]]$gsc), 
                       plot = FALSE)
ddu <- consensusScores(resList = reslist, class = "distinct", 
                       adjusted = TRUE, direction = "up", method = "median", 
                       n = length(reslist[[1]]$gsc), plot = FALSE)
ddd <- consensusScores(resList = reslist, class = "distinct", 
                       adjusted = TRUE, direction = "down", method = "median", 
                       n = length(reslist[[1]]$gsc), plot = FALSE)
mdu <- consensusScores(resList = reslist, class = "mixed", 
                       adjusted = TRUE, direction = "up", method = "median", 
                       n = length(reslist[[1]]$gsc), plot = FALSE)
mdd <- consensusScores(resList = reslist, class = "mixed", 
                       adjusted = TRUE, direction = "down", method = "median", 
                       n = length(reslist[[1]]$gsc), plot = FALSE)
nodPval <- nod$pMat
dduPval <- ddu$pMat
dddPval <- ddd$pMat
mduPval <- mdu$pMat
mddPval <- mdd$pMat
nod <- nod$rankMat
ddu <- ddu$rankMat
ddd <- ddd$rankMat
mdu <- mdu$rankMat
mdd <- mdd$rankMat
nod <- cbind(rownames(nod), nod[, 2])
colnames(nod) <- c("", "nod")
ddu <- cbind(rownames(ddu), ddu[, 2])
colnames(ddu) <- c("", "ddu")
ddd <- cbind(rownames(ddd), ddd[, 2])
colnames(ddd) <- c("", "ddd")
mdu <- cbind(rownames(mdu), mdu[, 2])
colnames(mdu) <- c("", "mdu")
mdd <- cbind(rownames(mdd), mdd[, 2])
colnames(mdd) <- c("", "mdd")
tmp <- merge(nod, ddu, by = 1)
tmp <- merge(tmp, ddd, by = 1)
tmp <- merge(tmp, mdu, by = 1)
tmp <- merge(tmp, mdd, by = 1)
consRankMat <- tmp[, 2:6]
rownames(consRankMat) <- tmp[, 1]
for (i in 1:ncol(consRankMat)) {
  consRankMat[, i] <- as.numeric(as.character(consRankMat[,i]))
}
consRankMat <- as.matrix(consRankMat)
plotmat <- consRankMat[which(apply(consRankMat, 1, min) <= 5), ]
colnames(plotmat) <- c("Non-directional", "Distinct-directional (up)", 
                       "Distinct-directional (dn)", "Mixed-directional (up)", 
                       "Mixed-directional (dn)")
myorder <- c(3, 5, 1, 4, 2)
gs_names <- rownames(plotmat[, myorder])
pMat <- matrix(ncol = 5, nrow = length(gs_names))
colnames(pMat) <- c("Distinct-directional (dn)", "Mixed-directional (dn)", 
                    "Non-directional", "Mixed-directional (up)", "Distinct-directional (up)")
rownames(pMat) <- gs_names
iGs <- match(gs_names, rownames(dddPval))
pMat[, 1] <- apply(dddPval[iGs, ], 1, median, na.rm = TRUE)
iGs <- match(gs_names, rownames(mddPval))
pMat[, 2] <- apply(mddPval[iGs, ], 1, median, na.rm = TRUE)
iGs <- match(gs_names, rownames(nodPval))
pMat[, 3] <- apply(nodPval[iGs, ], 1, median, na.rm = TRUE)
iGs <- match(gs_names, rownames(mduPval))
pMat[, 4] <- apply(mduPval[iGs, ], 1, median, na.rm = TRUE)
iGs <- match(gs_names, rownames(dduPval))
pMat[, 5] <- apply(dduPval[iGs, ], 1, median, na.rm = TRUE)
#cell color
set.seed(1)
tmp <- round(rnorm(10000, 0, 1000))
tmp <- rev(unique(sort(tmp[tmp > 0])))
clrs <- c("red3", "red", "orange", "yellow", "lightyellow", 
            "white")
mycol <- colorRampPalette(rev(clrs), interpolate = "linear")(max(tmp))[tmp]
#data matrix for heatmap
tmpMat <- plotmat
yclust <- hclust(dist(tmpMat))
tmpMat = data.frame(tmpMat[, myorder])
tmpMat$ID = rownames(tmpMat)
plotData = gather(tmpMat, group, Rank, -ID) %>%
  mutate(group = factor(group, levels = colnames(tmpMat)))
pMat = data.frame(pMat) #pvalue
pMat$ID = rownames(pMat)
pData = gather(pMat, group, pvalue, -ID) %>%
  mutate(group = factor(group, levels = colnames(pMat)))
#generate heatmap
ggplot(plotData, aes(group, ID)) +
  geom_tile(aes(fill = Rank)) +
  scale_fill_gradientn(colors = mycol) +
  geom_text(data = pData, aes(group, ID, label = scales::scientific(pvalue, digits = 3)), size = 3) +
  # geom_text(aes(label = Rank)) +
  scale_y_dendrogram(hclust = yclust) +
  theme_bw()+
  theme(
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, size = 10,hjust = 1),
    axis.text.y = element_text(size = 10)
  ) 
ggsave(file = "figureS1B.pdf", width = 6, height = 6)

