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
library(broom)
source("functions2.R")

df.hclust.1 <- read.xlsx("Compare_K_Hclust/hclust_clusters.xlsx", sheet = "Cluster1")
df.hclust.2 <- read.xlsx("Compare_K_Hclust/hclust_clusters.xlsx", sheet = "Cluster2")
df.kmeans.1 <- read.xlsx("Compare_K_Hclust/kmeans_clusters.xlsx", sheet = "Cluster1")
df.kmeans.2 <- read.xlsx("Compare_K_Hclust/kmeans_clusters.xlsx", sheet = "Cluster2")

df.inner.join.1 <- inner_join(df.hclust.1, df.kmeans.1, by = "biomarker_key_SUBJECTID")

df.inner.join.2 <- inner_join(df.hclust.2, df.kmeans.2, by = "biomarker_key_SUBJECTID")


# PCA on hclust clusters

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

map_all <- pheatmap(
  transpose_u(df0.p.map),
  cluster_rows = F, cluster_cols = TRUE,
  cellwidth = 5,
  cellheight = 5,
  fontsize = 5,
  cutree_cols = 2
)

# PCA Start

res.clust <- rbind(transpose_u(df0.p), cluster = cutree(map_all$tree_col, k = 2))
res.clust <- transpose_u(res.clust)

res.clust[std_biomarkers] <- lapply(res.clust[std_biomarkers], as.numeric)

res.pca <- PCA(res.clust[std_biomarkers], graph = F)
fviz_screeplot(res.pca, addlabels = T, ylim = c(0, 50))
var_pca <- get_pca_var(res.pca)
fviz_contrib(res.pca, choice = "var", axes = 1)
fviz_contrib(res.pca, choice = "var", axes = 2)
fviz_pca_ind(res.pca,
             label = "none",
             col.ind = as.character(res.clust$cluster),
             addEllipses =T)
# PCA End


# PCA on Kmeans clusters

df0.p2 <- df0.p %>%
  select(all_of(std_biomarkers))
df0.p2[std_biomarkers] <- lapply(df0.p2[std_biomarkers], as.numeric)
k5 <- kmeans(df0.p2, centers = 2, nstart = 25)
res.clust <- cbind(df0.p, cluster = k5$cluster)

# PCA Start
res.clust.pca <- res.clust %>%
  select(all_of(std_biomarkers))
res.clust.pca[std_biomarkers] <- lapply(res.clust.pca[std_biomarkers], as.numeric)

res.pca <- PCA(res.clust.pca, graph = F)
fviz_screeplot(res.pca, addlabels = T, ylim = c(0, 50))
var_pca <- get_pca_var(res.pca)
fviz_contrib(res.pca, choice = "var", axes = 1)
fviz_contrib(res.pca, choice = "var", axes = 2)
fviz_pca_ind(res.pca,
             label = "none",
             col.ind = as.character(res.clust$cluster),
             addEllipses =T)
# PCA End