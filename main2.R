library("dplyr")
library("openxlsx")
library("tidyverse")
library("ggplot2")
library("FactoMineR")
library("factoextra")
library("corrplot")
library(Hmisc)
library(missMDA)
library(vegan)
library(pheatmap)
library(data.table)
library(cluster)
library(dendextend)
library(ggplot2)
library(gridExtra)
library(viridis)
source("functions.R")

df0 = read.xlsx("PROPPR_longitudinal_data_dictionary_edm_5.13.20.xlsx", sheet = "timepoint_0")
df0 <- df0 %>%
  filter(INJ_MECH != "Both Types of Injury")


df0.p <- pre_process2(df0)
std_cols <- c(234 : 276)
std_biomarkers <- colnames(df0.p[std_cols])
cols <- c(51:93)
biomarkers <- colnames(df0.p[std_cols])

df0.p <- remove_na_patients2(df0.p)

df0.p[std_biomarkers] <- lapply(df0.p[std_biomarkers], as.numeric)

df0.p.t <- transpose_u(df0.p %>%
  select(all_of(std_biomarkers)))

# Remove the 4 outliers
df0.p.t$`576` <- NULL
df0.p.t$`587` <- NULL
df0.p.t$`355` <- NULL
df0.p.t$`579` <- NULL

color <- viridis(50)
breaks.low <- c(seq(-8, -5, by = 3))
rampcolors.low <- colorRampPalette(color[1:2])(length(breaks.low))
breaks.mid <- c(seq(-4.9, 4.9, by = 0.1))
rampcolors.mid <- colorRampPalette(color[2:49])(length(breaks.mid))
breaks.high <- c(seq(5, 8, by = 3))
rampcolors.high <- colorRampPalette(color[49:50])(length(breaks.high))

breaks <- c(breaks.low, breaks.mid, breaks.high)
rampcolors <- c(rampcolors.low, rampcolors.mid, rampcolors.high)

map_all <- pheatmap(
  df0.p.t,
  cluster_rows = F, cluster_cols = TRUE,
  cellwidth = 1,
  cellheight = 5,
  fontsize = 5,
  color = rampcolors,
  breaks = breaks,
  filename = "R/heatmap_allv2.png" 
)

fviz_nbclust(transpose_u(df0.p.t), hcut, method = "silhouette")

res.clust <- rbind(transpose_u(df0.p), cluster = cutree(map_all$tree_col, k = 2))
res.clust <- transpose_u(res.clust)

res.clust[std_biomarkers] <- lapply(res.clust[std_biomarkers], as.numeric)

df0.c1 <- res.clust %>%
  filter(cluster == 1)
df0.c2 <- res.clust %>%
  filter(cluster == 2)
df0.c3 <- res.clust %>%
  filter(cluster == 3)
df0.c4 <- res.clust %>%
  filter(cluster == 4)
df0.c5 <- res.clust %>%
  filter(cluster == 5)
df0.c6 <- res.clust %>%
  filter(cluster == 6)

df0.c1.biomarkers <- df0.c1 %>%
  select(all_of(biomarkers))
df0.c2.biomarkers <- df0.c2 %>%
  select(all_of(biomarkers))

df0.c1.biomarkers <- as.data.frame(lapply(df0.c1.biomarkers, as.numeric))
df0.c2.biomarkers <- as.data.frame(lapply(df0.c2.biomarkers, as.numeric))

t.test(df0.c1.biomarkers[0], df0.c2.biomarkers[0])

cluster.list <- list(df0.c1, df0.c2) #, df0.c3, df0.c4, df0.c5, df0.c5)

outcomes <- colnames(df0[173:199])
outcomes <- append(outcomes, c("RBC_10", "RAN_3HRST", "RAN_24HRST"))

res.list <- compare_perc_occur(cluster.list, outcomes)

plots <- make_graphs(res.list$df, res.list$categories)

n <- length(plots)
nCol <- floor(sqrt(n))
g <- arrangeGrob(grobs = plots, ncol = nCol)
ml <- marrangeGrob(grobs = plots, nrow = 1, ncol = 2)
ggsave(file="occurrence.pdf", ml)
