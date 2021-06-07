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
library("openxlsx")

df = read.csv('df_with_cluster_map.csv', stringsAsFactors = F)

biomarker_cols <- df[52:94] #df0[8:50] 
# biomarker_cols = df0[8:50]
biomarkers <- colnames(biomarker_cols)

for (x in biomarkers) {
  df[paste(x, "_zscore")] <- scale(df[x])
}

col_n = c(252:294)
biomarkers_std <- colnames(df[col_n])

std_cols_cluster1 <- df %>%
  filter(Clusters == "first") %>%
  select(all_of(biomarkers_std))

std_cols_cluster2 <- df %>%
  filter(Clusters == "second") %>%
  select(all_of(biomarkers_std))

std1_mean <- std_cols_cluster1 %>%
  summarise_all(mean, na.rm = T)

std1_mean_t <- transpose(std1_mean)
colnames(std1_mean_t) <- "Cluster 1"
rownames(std1_mean_t) <- colnames(std1_mean)

std2_mean <- std_cols_cluster2 %>%
  summarise_all(mean, na.rm = T)

std2_mean_t <- transpose(std2_mean)
colnames(std2_mean_t) <- "Cluster 2"
rownames(std2_mean_t) <- colnames(std2_mean)

all_mean <- cbind(std1_mean_t, std2_mean_t)

map <- pheatmap(
  all_mean,
  cluster_rows = FALSE, cluster_cols = FALSE,
  cellwidth = 10,
  cellheight = 10,
  fontsize = 10,
  filename = "heatmap.png"
)


