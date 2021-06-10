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
source("pca.R")

df0 = read.xlsx("PROPPR_longitudinal_data_dictionary_edm_5.13.20.xlsx", sheet = "timepoint_0")

dataset_info <- count_data(df0)
# View / print out dataframes.
View(dataset_info$num_gender)
View(dataset_info$num_sepsis)
View(dataset_info$num_ards)
View(dataset_info$num_death)
View(dataset_info$num_gender_inj)

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

heirarch_res <- heirarch_cluster(df0)
print(heirarch_res$cluster_plot)

heirarch_sep_res <- heirarch_cluster_sep(df0)
print(heirarch_sep_res$blunt_plot)
print(heirarch_sep_res$pen_plot)
