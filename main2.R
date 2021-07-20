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
source("functions2.R")

df0 = read.xlsx("PROPPR_longitudinal_data_dictionary_edm_5.13.20.xlsx", sheet = "timepoint_0")
df0 <- df0 %>%
  filter(INJ_MECH != "Both Types of Injury")


res_pre <- pre_process(df0, T)
df0.p <- res_pre$df0
std_biomarkers <- res_pre$std_cols
biomarkers <- res_pre$cols

df0.p[std_biomarkers] <- lapply(df0.p[std_biomarkers], as.numeric)

df0.p.t <- transpose_u(df0.p)

# Remove the 4 outliers
df0.p.t$`576` <- NULL
df0.p.t$`587` <- NULL
df0.p.t$`355` <- NULL
df0.p.t$`579` <- NULL

df0.p.map <- transpose_u(df0.p.t) %>% select(all_of(std_biomarkers))
df0.p.map <- as.data.frame(lapply(df0.p.map, as.numeric))

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
  transpose_u(df0.p.map),
  cluster_rows = F, cluster_cols = TRUE,
  cellwidth = 5,
  cellheight = 5,
  fontsize = 5,
  color = rampcolors,
  breaks = breaks,
  cutree_cols = 2,
  filename = "R/heatmap_allv4.png" 
)

fviz_nbclust(df0.p.map, hcut, method = "silhouette")

res.clust <- rbind(transpose_u(df0.p), cluster = cutree(map_all$tree_col, k = 2))
res.clust <- transpose_u(res.clust)

res.clust[std_biomarkers] <- lapply(res.clust[std_biomarkers], as.numeric)

df0.c1 <- res.clust %>%
  filter(cluster == 1)
df0.c2 <- res.clust %>%
  filter(cluster == 2)

df0.c1.export <- df0.c1 %>%
  select(biomarker_key_SUBJECTID, cluster)
df0.c2.export <- df0.c2 %>%
  select(biomarker_key_SUBJECTID, cluster)

df0.c1.biomarkers <- df0.c1 %>%
  select(all_of(biomarkers))
df0.c2.biomarkers <- df0.c2 %>%
  select(all_of(biomarkers))

df0.c1.biomarkers <- as.data.frame(lapply(df0.c1.biomarkers, as.numeric))
df0.c2.biomarkers <- as.data.frame(lapply(df0.c2.biomarkers, as.numeric))

df0.c1.stdb <- df0.c1 %>%
  select(all_of(std_biomarkers))
df0.c2.stdb <- df0.c2 %>%
  select(all_of(std_biomarkers))

df0.c1.stdb <- as.data.frame(lapply(df0.c1.stdb, as.numeric))
df0.c2.stdb <- as.data.frame(lapply(df0.c2.stdb, as.numeric))

df0.c1.export <- as.data.frame(lapply(df0.c1.export, as.numeric))
df0.c2.export <- as.data.frame(lapply(df0.c2.export, as.numeric))

df.export.list <- list("Cluster1" = df0.c1.export, "Cluster2" = df0.c2.export)
write.xlsx(df.export.list, file = "hclust_clusters.xlsx")

df0.blunt <- res.clust %>%
  filter(INJ_MECH == "Blunt Injury Only")

blunt_clusters <- df0.blunt %>%
  group_by(cluster) %>%
  count(cluster)

pen_clusters <- df0.pen %>%
  group_by(cluster) %>%
  count(cluster)

df0.pen <- res.clust %>%
  filter(INJ_MECH == "Penetrating Injury Only")

df0.blunt.biomarkers <- df0.blunt %>%
  select(all_of(biomarkers))
df0.pen.biomarkers <- df0.pen %>%
  select(all_of(biomarkers))

df0.blunt.biomarkers <- as.data.frame(lapply(df0.blunt.biomarkers, as.numeric))
df0.pen.biomarkers <- as.data.frame(lapply(df0.pen.biomarkers, as.numeric))
 
clust1 <- pheatmap(
  transpose_u(df0.c1.stdb),
  cluster_rows = F, cluster_cols = F,
  cellwidth = 1,
  cellheight = 5,
  fontsize = 5,
  color = rampcolors,
  breaks = breaks,
  filename = "R/heatmap_clust1.png" 
)

clust2 <- pheatmap(
  transpose_u(df0.c2.stdb),
  cluster_rows = F, cluster_cols = F,
  cellwidth = 1,
  cellheight = 5,
  fontsize = 5,
  color = rampcolors,
  breaks = breaks,
  filename = "R/heatmap_clust2.png" 
)

cluster.list <- list(df0.c1, df0.c2) #, df0.c3, df0.c4, df0.c5, df0.c5)

outcomes <- colnames(df0[173:199])
outcomes <- append(outcomes, c("RBC_10", "RAN_3HRST", "RAN_24HRST"))

res.list <- compare_perc_occur(cluster.list, outcomes)

plots <- make_graphs(res.list$df, res.list$categories)

n <- length(plots)
nCol <- floor(sqrt(n))
g <- arrangeGrob(grobs = plots, ncol = nCol)
ml <- marrangeGrob(grobs = plots, nrow = 1, ncol = 2)
ggsave(file="occurrence_hclust.pdf", ml)


test.res <- compare_t_test(df0.c1.biomarkers, df0.c2.biomarkers)
test.res.all <- data.frame(biomarker = rep(names(test.res)), bind_rows(test.res))

test.inj.res <- compare_t_test(df0.blunt.biomarkers, df0.pen.biomarkers)
test.inj.res.all <- data.frame(biomarker = rep(names(test.inj.res)), bind_rows(test.inj.res))

test.res.all <- test.res.all %>%
  filter(p.value < (0.05 / 36))

test.inj.res.all <- test.inj.res.all %>%
  filter(p.value < (0.05 / 36))

print(test.res.all)
write.xlsx(test.inj.res.all, file = "t_test_hclust_inj.xlsx")
write.xlsx(test.res.all, file = "t_test_hclust.xlsx")
