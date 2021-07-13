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


res_pre <- pre_process3(df0)
df0.p <- res_pre$df0
std_biomarkers <- res_pre$cols
cols <- c(51:93)
biomarkers <- colnames(df0.p[cols])

df0.p <- remove_na_patients2(df0.p)

df0.p[std_biomarkers] <- lapply(df0.p[std_biomarkers], as.numeric)

df0.p.t <- transpose_u(df0.p %>%
                         select(all_of(std_biomarkers)))

# Remove the 4 outliers
df0.p.t$`576` <- NULL
df0.p.t$`587` <- NULL
df0.p.t$`355` <- NULL
df0.p.t$`579` <- NULL

res.pca <- PCA(df0.p.t, graph = FALSE)
print(get_eig(res.pca))
