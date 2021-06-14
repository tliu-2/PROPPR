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
source("functions.R")

df0 = read.xlsx("PROPPR_longitudinal_data_dictionary_edm_5.13.20.xlsx", sheet = "timepoint_0")

dataset_info <- count_data(df0)
# View / print out dataframes.
View(dataset_info$mean_age_gender)
View(dataset_info$num_gender)
View(dataset_info$num_sepsis)
View(dataset_info$num_ards)
View(dataset_info$num_death)
View(dataset_info$num_gender_inj)
View(dataset_info$num_inj_type)


table1_data <- get_table1_data(df0)
print(table1_data$all_age_mean)
print(table1_data$all_age_sd)
print(table1_data$blunt_age_mean)
print(table1_data$blunt_age_sd)
print(table1_data$pen_age_mean)
print(table1_data$pen_age_sd)
print(table1_data$num_sex)

res_list <- pca_plot(df0)

# View plots from function
print(res_list$scree_plot)
print(res_list$var_contrib)
print(res_list$plot)
print(res_list$plot_ellipse)

# View k-means plots
kmeans_res <- k_means(df0)
print(kmeans_res$silhouette_plot)
print(kmeans_res$kmeans_plot)

# Conduct k-means clustering on only blunt and only penetrating injury.
kmeans_cluster_res <- k_means_sep(df0)
print(kmeans_cluster_res$Blunt_Heat)
print(kmeans_cluster_res$Pen_Heat)
print(kmeans_cluster_res$Blunt_Cluster)
print(kmeans_cluster_res$Pen_Cluster)

hierarch_res <- hierarch_cluster(df0)
print(hierarch_res$cluster_plot)

hierarch_sep_res <- hierarch_cluster_sep(df0)
print(hierarch_sep_res$blunt_plot)
print(hierarch_sep_res$pen_plot)
hierarch_sep_res$dendro_b_agg
hierarch_sep_res$dendro_p_agg

df_b <- hierarch_sep_res$df_b
df_p <- hierarch_sep_res$df_p

df_b %>%
  select(cluster) %>%
  count(cluster)

gen_heatmaps(df0)

cluster_res <- cluster_heatmap(df_b, df_p)
print(cluster_res$mapb)
print(cluster_res$mapp)


# Try clustering by rows in heatmap.
df0.p <- pre_process(df0)
rownames.df0 <- rownames(df0.p)
colnames.df0 <- colnames(df0.p)
df0.t <- transpose(df0.p)
rownames(df0.t) <- colnames.df0
colnames(df0.t) <- rownames.df0

map_all <- pheatmap(
  df0.t,
  cluster_rows = TRUE, cluster_cols = TRUE,
  cellwidth = 10,
  cellheight = 10,
  fontsize = 10,
  filename = "R/heatmap_all.png" 
)

fviz_nbclust(df0.t, kmeans, method = "silhouette")

row.clust <- cutree(map_all$tree_row, k = 8)
res.clust <- rbind(df0.t, cluster = cutree(map_all$tree_row, k = 8))
res.clust.t <- transpose(res.clust)
rownames(res.clust.t) <- colnames(res.clust)
colnames(res.clust.t) <- rownames(res.clust)

biomarkers <- c(0:36)
res.1.t <- res.clust.t %>%
  filter(cluster == 1) %>%
  select(all_of(biomarkers))

res.1 <- transpose(res.1.t)
rownames(res.1) <- colnames(res.1.t)
colnames(res.1) <- rownames(res.1.t)

res.2.t <- res.clust.t %>%
  filter(cluster == 2) %>%
  select(all_of(biomarkers))

res.2 <- transpose(res.2.t)
rownames(res.2) <- colnames(res.2.t)
colnames(res.2) <- rownames(res.2.t)


res.3.t <- res.clust.t %>%
  filter(cluster == 3) %>%
  select(all_of(biomarkers))

res.3 <- transpose(res.3.t)
rownames(res.3) <- colnames(res.3.t)
colnames(res.3) <- rownames(res.3.t)


res.4.t <- res.clust.t %>%
  filter(cluster == 4) %>%
  select(all_of(biomarkers))

res.4 <- transpose(res.4.t)
rownames(res.4) <- colnames(res.4.t)
colnames(res.4) <- rownames(res.4.t)


res.5.t <- res.clust.t %>%
  filter(cluster == 5) %>%
  select(all_of(biomarkers))

res.5 <- transpose(res.5.t)
rownames(res.5) <- colnames(res.5.t)
colnames(res.5) <- rownames(res.5.t)


res.6.t <- res.clust.t %>%
  filter(cluster == 6) %>%
  select(all_of(biomarkers))

res.6 <- transpose(res.6.t)
rownames(res.6) <- colnames(res.6.t)
colnames(res.6) <- rownames(res.6.t)


res.7.t <- res.clust.t %>%
  filter(cluster == 7) %>%
  select(all_of(biomarkers))

res.7 <- transpose(res.7.t)
rownames(res.7) <- colnames(res.7.t)
colnames(res.7) <- rownames(res.7.t)


res.8.t <- res.clust.t %>%
  filter(cluster == 8) %>%
  select(all_of(biomarkers))

res.8 <- transpose(res.8.t)
rownames(res.8) <- colnames(res.8.t)
colnames(res.8) <- rownames(res.8.t)


pheatmap(
  res.1,
  cluster_rows = TRUE, cluster_cols = TRUE,
  cellwidth = 10,
  cellheight = 10,
  fontsize = 10,
  filename = "R/test1.png" 
)

