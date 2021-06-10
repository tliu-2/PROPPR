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

df0 = read.xlsx("PROPPR_longitudinal_data_dictionary_edm_5.13.20.xlsx", sheet = "timepoint_0")

# PREPROCESSING
biomarker_cols <- df0[51:93] #df0[8:50] 
# biomarker_cols = df0[8:50]
biomarkers = colnames(biomarker_cols)

df0 <- df0 %>%
  filter(INJ_MECH != "Both Types of Injury") 

for (x in biomarkers) {
  if (sum(is.na(df0[x])) / nrow(df0[x]) > 0.2) {
    df0[x] <- NULL
  }
  if(sum(is.na(biomarker_cols[x])) / nrow(biomarker_cols[x]) > 0.2) {
    biomarker_cols[x] <- NULL
  }
}

biomarkers <- colnames(biomarker_cols)

# Impute missing values.
for (x in biomarkers) {
  df0[x] <- with(df0, impute(df0[x], mean))
}

# Divide based on injury mechanism

df_blunt <- df0 %>%
  filter(INJ_MECH == "Blunt Injury Only") %>%
  select(all_of(biomarkers))
df_pen <- df0 %>%
  filter(INJ_MECH == "Penetrating Injury Only")%>%
  select(all_of(biomarkers))

for (x in biomarkers) {
  df_blunt[paste(substr(x, 5, nchar(x)))] <- scale(df_blunt[x])
  df_pen[paste(substr(x, 5, nchar(x)))] <- scale(df_pen[x])
}

std_cols <- c(37:72)
std_biomarkers <- colnames(df_blunt[std_cols])
df_blunt <- df_blunt[std_cols]
df_pen <- df_pen[std_cols]

# K Means on Blunt
fviz_nbclust(df_blunt, kmeans, method = "silhouette")
set.seed(42)
km.res_blunt <- kmeans(df_blunt, 2, nstart = 25)
fviz_cluster(km.res_blunt, data = df_blunt,
             elipse.type = "convex", palette = "jco", repel = TRUE,
             ggtheme = theme_minimal()
)

df_blunt["Cluster"] <- km.res_blunt$cluster
df.b.1 <- df_blunt %>%
  filter(Cluster == "1") %>%
  select(all_of(std_biomarkers)) %>%
  summarise_all(mean, na.rm=T)

df.b.2 <- df_blunt %>%
  filter(Cluster == "2") %>%
  select(all_of(std_biomarkers)) %>%
  summarise_all(mean, na.rm=T)

df.b.all <- rbind(df.b.1, df.b.2)
df.b.t <- transpose(df.b.all)
rownames(df.b.t) <- colnames(df.b.all)
colnames(df.b.t) <- rownames(df.b.all)
pheatmap(
  df.b.t,
  cluster_rows = FALSE, cluster_cols = FALSE,
  cellwidth = 10,
  cellheight = 10,
  fontsize = 10,
)


# K means on Pen
fviz_nbclust(df_blunt, kmeans, method = "silhouette")
set.seed(42)
km.res_pen <- kmeans(df_pen, 2, nstart = 25)
fviz_cluster(km.res_pen, data = df_pen,
             elipse.type = "convex", palette = "jco", repel = TRUE,
             ggtheme = theme_minimal()
)

df_pen["Cluster"] <- km.res_pen$cluster
df.p.1 <- df_pen %>%
  filter(Cluster == "1") %>%
  select(all_of(std_biomarkers)) %>%
  summarise_all(mean, na.rm=T)

df.p.2 <- df_pen %>%
  filter(Cluster == "2") %>%
  select(all_of(std_biomarkers)) %>%
  summarise_all(mean, na.rm=T)

df.p.all <- rbind(df.p.1, df.p.2)
df.p.t <- transpose(df.p.all)
rownames(df.p.t) <- colnames(df.p.all)
colnames(df.p.t) <- rownames(df.p.all)
pheatmap(
  df.p.t,
  cluster_rows = FALSE, cluster_cols = FALSE,
  cellwidth = 10,
  cellheight = 10,
  fontsize = 10,
)
