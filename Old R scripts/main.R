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
df0 <- df0 %>%
  filter(INJ_MECH != "Both Types of Injury")
rownames.df0 <- rownames(df0.p)
colnames.df0 <- colnames(df0.p)
df0.t <- transpose(df0.p)
rownames(df0.t) <- colnames.df0
colnames(df0.t) <- rownames.df0

map_all <- pheatmap(
  df0.t,
  cluster_rows = TRUE, cluster_cols = TRUE,
  cellwidth = 1,
  cellheight = 10,
  fontsize = 10,
  color = viridis(50),
  filename = "R/heatmap_all.png" 
)

fviz_nbclust(df0.t, hcut, method = "silhouette")
fviz_nbclust(df0.p, hcut, method = "silhouette")

col.clust <- cutree(map_all$tree_col, k = 2)
col.clust <- cbind(df0.p, cluster = cutree(map_all$tree_col, k = 2))
col.clust <- cbind(col.clust, Survival = df0$`_30DAYST`)

# Examine any clusters inside the original clusters.
df0.t2 <- transpose_u(df0)
res.clust.all <- rbind(df0.t2, cluster = cutree(map_all$tree_row, k = 2))
res.clust.all <- transpose_u(res.clust.all)

biomarker_cols <- df0[51:93]
biomarkers <- colnames(biomarker_cols)
res.clust.all[biomarkers] <- lapply(res.clust.all[biomarkers], as.numeric)
res.clust.all$AGE <- as.numeric(res.clust.all$AGE)
res.clust.all <- pre_process2(res.clust.all)

res.clust.all.mean <- res.clust.all %>%
  group_by(cluster) %>%
  select(all_of(biomarkers)) %>%
  summarise_all(mean, na.rm=T) %>%
  ungroup()
res.clust.all.mean$cluster <- NULL

# Clusters on columns of heatmap.
pheatmap(
  transpose_u(res.clust.all.mean),
  cluster_rows = FALSE, cluster_cols = FALSE,
  cellwidth = 10,
  cellheight = 10,
  fontsize = 10,
  filename = "R/test_sub.png" 
)

df0.c1 <- res.clust.all %>%
  filter(cluster == 1)
df0.c2 <- res.clust.all %>%
  filter(cluster == 2)

outcomes <- colnames(res.clust.all[173:199])
outcomes <- append(outcomes, c("RBC_10", "RAN_3HRST", "RAN_24HRST"))

res.list <- compare_perc_occur(df0.c1, df0.c2, outcomes, df0)
res1 <- data.frame(res.list$c1)
res2 <- data.frame(res.list$c2)

n <- length(res.list$plots)
nCol <- floor(sqrt(n))
g <- arrangeGrob(grobs = res.list$plots, ncol = nCol)
ml <- marrangeGrob(grobs = res.list$plots, nrow = 1, ncol = 2)
ggsave(file="occurrence.pdf", ml)


df.blunt <- res.clust.all %>%
  filter(INJ_MECH == "Blunt Injury Only")
df.pen <- res.clust.all %>%
  filter(INJ_MECH == "Penetrating Injury Only")
print(df.blunt %>%
        count(cluster))
print(df.pen %>%
        count(cluster))

c <- list()
for (x in 1:2) {
  c[[paste("res.", x, "l", sep = "")]] <- heatmap_process(col.clust, "Living", x, FALSE)
  c[[paste("res.", x, "d", sep = "")]] <- heatmap_process(col.clust, "Deceased", x, FALSE)
}

res.c.all <- bind_rows(c)
res.c.all.t <- transpose_u(res.c.all)

pheatmap(
  res.c.all.t,
  cluster_rows = FALSE, cluster_cols = FALSE,
  cellwidth = 10,
  cellheight = 10,
  fontsize = 10,
  filename = "R/test_cols.png" 
)

row.clust <- cutree(map_all$tree_row, k = 10)
res.clust <- rbind(df0.t, cluster = cutree(map_all$tree_row, k = 10))
res.clust.t <- transpose_u(res.clust)
res.clust.t <- cbind(res.clust.t, Survival = df0$`_30DAYST`)

res.clust.all <- rbind(df0.t2, cluster = cutree(map_all$tree_row, k = 10))
res.clust.all <- transpose_u(res.clust.all)

res.clust.all[biomarkers] <- lapply(res.clust.all[biomarkers], as.numeric)
res.clust.all <- pre_process2(res.clust.all)
res.clust.all.mean <- res.clust.all %>%
  group_by(cluster) %>%
  select(all_of(biomarkers)) %>%
  summarise_all(mean, na.rm=T) %>%
  ungroup()
res.clust.all.mean$cluster <- NULL

# Clusters on rows
pheatmap(
  transpose_u(res.clust.all.mean),
  cluster_rows = FALSE, cluster_cols = FALSE,
  cellwidth = 10,
  cellheight = 10,
  fontsize = 10,
  filename = "R/test_sub2.png" 
)

df.list <- list()
for (i in seq_along(1:10)) {
  df.list[[i]] <- res.clust.all %>%
    filter(cluster == i)
}



r <- list()
for (x in 1:10) {
    r[[paste("res.", x, "l", sep = "")]] <- heatmap_process(res.clust.t, "Living", x)
    r[[paste("res.", x, "d", sep = "")]] <- heatmap_process(res.clust.t, "Deceased", x)
}

res.r.all <- bind_cols(r)



clust.all.map.r <- pheatmap(
  res.r.all,
  cluster_rows = FALSE, cluster_cols = FALSE,
  cellwidth = 10,
  cellheight = 10,
  fontsize = 10,
  filename = "R/testr.png" 
)


clust.all.map.c <- pheatmap(
  res.c.all,
  cluster_rows = FALSE, cluster_cols = FALSE,
  cellwidth = 10,
  cellheight = 10,
  fontsize = 10,
  filename = "R/testc.png" 
)

