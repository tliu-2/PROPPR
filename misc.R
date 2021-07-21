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

df.hclust.1 <- read.xlsx("Compare_K_Hclust/compare_clusters.xlsx", sheet = "Cluster1.H")
df.hclust.2 <- read.xlsx("Compare_K_Hclust/compare_clusters.xlsx", sheet = "Cluster2.H")
df.kmeans.1 <- read.xlsx("Compare_K_Hclust/compare_clusters.xlsx", sheet = "Cluster1.K")
df.kmeans.2 <- read.xlsx("Compare_K_Hclust/compare_clusters.xlsx", sheet = "Cluster2.K")

df.inner.join.1 <- inner_join(df.hclust.1, df.kmeans.1, by = "biomarker_key_SUBJECTID")

df.inner.join.2 <- inner_join(df.hclust.2, df.kmeans.2, by = "biomarker_key_SUBJECTID")


df0 = read.xlsx("PROPPR_longitudinal_data_dictionary_edm_5.13.20.xlsx", sheet = "timepoint_0")
df0 <- df0 %>%
  filter(INJ_MECH != "Both Types of Injury")

# Preprocess
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

# End preprocess

# PCA on hclust clusters

df0.p.map <- transpose_u(df0.p.t) %>% select(all_of(std_biomarkers))
df0.p.map <- as.data.frame(lapply(df0.p.map, as.numeric))

map.all.h <- pheatmap(
  transpose_u(df0.p.map),
  cluster_rows = F, cluster_cols = TRUE,
  cellwidth = 5,
  cellheight = 5,
  fontsize = 5,
  cutree_cols = 2
)

# PCA Start

clusters <- cutree(map.all.h$tree_col, k = 2)
res.clust.h <- cbind(transpose_u(df0.p.t), cluster = cutree(map.all.h$tree_col, k = 2))


res.clust.h[std_biomarkers] <- lapply(res.clust.h[std_biomarkers], as.numeric)

res.pca.h <- PCA(res.clust.h[std_biomarkers], graph = F)
fviz_screeplot(res.pca.h, addlabels = T, ylim = c(0, 50))
var_pca <- get_pca_var(res.pca.h)
fviz_contrib(res.pca.h, choice = "var", axes = 1)
fviz_contrib(res.pca.h, choice = "var", axes = 2)
fviz_pca_ind(res.pca.h,
             label = "none",
             col.ind = as.character(res.clust.h$cluster),
             addEllipses =T)
# PCA End

df0.c1.h <- res.clust.h %>%
  filter(cluster == 1)
df0.c2.h <- res.clust.h %>%
  filter(cluster == 2)



# PCA on Kmeans clusters

df0.p <- transpose_u(df0.p.t)

df0.p2 <- df0.p %>%
  select(all_of(std_biomarkers))
df0.p2[std_biomarkers] <- lapply(df0.p2[std_biomarkers], as.numeric)

map.all.k <- pheatmap(
  df0.p.map,
  kmeans_k=2,
  cluster_rows = F, cluster_cols = F,
  cellwidth = 5,
  cellheight = 5,
  fontsize = 5,
  cutree_cols = 2
)


clusters.k <- map.all.k$kmeans$cluster
res.clust.k <- cbind(df0.p, cluster = map.all.k$kmeans$cluster)

df0.c1.k <- res.clust.k %>% 
  filter(cluster == 1)
df0.c2.k <- res.clust.k %>%
  filter(cluster == 2)

# PCA Start
res.clust.pca.k <- res.clust.k %>%
  select(all_of(std_biomarkers))
res.clust.pca.k[std_biomarkers] <- lapply(res.clust.pca.k[std_biomarkers], as.numeric)

res.pca.k <- PCA(res.clust.pca.k, graph = F)
fviz_screeplot(res.pca.k, addlabels = T, ylim = c(0, 50))
var_pca.k <- get_pca_var(res.pca)
fviz_contrib(res.pca.k, choice = "var", axes = 1)
fviz_contrib(res.pca.k, choice = "var", axes = 2)
fviz_pca_ind(res.pca.k,
             label = "none",
             col.ind = as.character(res.clust.k$cluster),
             addEllipses =T)
# PCA End


# Compare clusters in hclust and kmeans
df0.c1.export.hclust <- df0.c1.h %>%
  select(biomarker_key_SUBJECTID, cluster)
df0.c2.export.hclust <- df0.c2.h %>%
  select(biomarker_key_SUBJECTID, cluster)

df0.c1.export.k <- df0.c1.k %>%
  select(biomarker_key_SUBJECTID, cluster)
df0.c2.export.k <- df0.c2.k %>%
  select(biomarker_key_SUBJECTID, cluster)
  
df.export.list <- list("Cluster1.H" = df0.c1.export.hclust, "Cluster2.H" = df0.c2.export.hclust,
                       "Cluster1.K" = df0.c1.export.k, "Cluster2.K" = df0.c2.export.k)
write.xlsx(df.export.list, file = "compare_clusters.xlsx")

# Occurence Output

outcomes <- colnames(df0[173:199])
outcomes <- append(outcomes, c("RBC_10", "RAN_3HRST", "RAN_24HRST"))

df0.c1 <- res.clust %>%
  filter(cluster == 1)
df0.c2 <- res.clust %>%
  filter(cluster == 2)

cluster.list <- list(df0.c1, df0.c2) #, df0.c3, df0.c4, df0.c5, df0.c5)

res.list <- compare_perc_occur(cluster.list, outcomes)

plots <- make_graphs(res.list$df, res.list$categories)

n <- length(plots)
nCol <- floor(sqrt(n))
g <- arrangeGrob(grobs = plots, ncol = nCol)
ml <- marrangeGrob(grobs = plots, nrow = 1, ncol = 2)
ggsave(file="occurence_kmeans.pdf", ml)
