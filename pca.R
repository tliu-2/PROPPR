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
#df2 = read.xlsx("PROPPR_longitudinal_data_dictionary_edm_5.13.20.xlsx", sheet = "timepoint_2")
#df4 = read.xlsx("PROPPR_longitudinal_data_dictionary_edm_5.13.20.xlsx", sheet = "timepoint_4")
#df6 = read.xlsx("PROPPR_longitudinal_data_dictionary_edm_5.13.20.xlsx", sheet = "timepoint_6")
#df12 = read.xlsx("PROPPR_longitudinal_data_dictionary_edm_5.13.20.xlsx", sheet = "timepoint_12")
#df24 = read.xlsx("PROPPR_longitudinal_data_dictionary_edm_5.13.20.xlsx", sheet = "timepoint_24")
#df48 = read.xlsx("PROPPR_longitudinal_data_dictionary_edm_5.13.20.xlsx", sheet = "timepoint_48")
#df72 = read.xlsx("PROPPR_longitudinal_data_dictionary_edm_5.13.20.xlsx", sheet = "timepoint_72")

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

for (x in biomarkers) {
  df0[x] <- with(df0, impute(df0[x], mean))
}

for (x in biomarkers) {
  df0[paste(x, "_zscore")] <- scale(df0[x])
}


std_cols <- c(227:262)
std_biomarkers <- colnames(df0[std_cols])

df_blunt <- df0 %>%
  filter(INJ_MECH == "Blunt Injury Only") %>%
  select(all_of(std_biomarkers))
df_pen <- df0 %>%
  filter(INJ_MECH == "Penetrating Injury Only")%>%
  select(all_of(std_biomarkers))

df_blunt_mean <- df_blunt %>%
  summarise_all(mean, na.rm = TRUE)
df_pen_mean <- df_pen %>%
  summarise_all(mean, na.rm = TRUE)

all_mean <- rbind(df_blunt_mean, df_pen_mean)
all_mean_t <- transpose(all_mean)
rownames(all_mean_t) <- colnames(all_mean)
colnames(all_mean_t) <- c("Blunt", "Penetrating")

# HEATMAP

map <- pheatmap(
  all_mean_t,
  cluster_rows = FALSE, cluster_cols = FALSE,
  cellwidth = 10,
  cellheight = 10,
  fontsize = 10,
  filename = "heatmap_inj_mech.png"
)


# PERMANOVA
df_c_dist <- dist(df_combined_std, method = "euclidean")
set.seed(42)
dif_c_div <- adonis2(df_c_dist~INJ_MECH, data=df0, permutations = 999, method="bray")
print(dif_c_div)

dispersion <- betadisper(df_c_dist, group=df0$INJ_MECH)
print(dispersion)
permutest(dispersion)
plot(dispersion, hull=FALSE, ellipse = TRUE)


# PCA
df0 = read.xlsx("PROPPR_longitudinal_data_dictionary_edm_5.13.20.xlsx", sheet = "timepoint_0")
df = read.xlsx("PROPPR_longitudinal_data_dictionary_edm_5.13.20.xlsx", sheet = "timepoint_0")
biomarker_cols <- df0[51:93] #df0[8:50] 
# biomarker_cols = df0[8:50]
biomarkers = colnames(biomarker_cols)

for (x in biomarkers) {
  if (sum(is.na(df0[x])) / nrow(df0[x]) > 0.2) {
    df0[x] <- NULL
  }
  if(sum(is.na(biomarker_cols[x])) / nrow(biomarker_cols[x]) > 0.2) {
    biomarker_cols[x] <- NULL
  }
}

biomarkers <- colnames(biomarker_cols)

for (x in biomarkers) {
  df0[paste(x, "_zscore")] <- scale(df0[x])
}

std_cols <- c(227:262)
std_biomarkers <- colnames(df0[std_cols])


df0 <- df0 %>%
  filter(INJ_MECH != "Both Types of Injury") %>%
  select(all_of(std_biomarkers))

df <- df %>%
  filter(INJ_MECH != "Both Types of Injury")

for (x in std_biomarkers) {
  df0[x] <- with(df0, impute(df0[x], mean))
}

combined_pca = PCA(df0, graph = FALSE, ncp = 25)
combined_eig_val <- get_eigenvalue(combined_pca)
print(combined_eig_val) # 15 dimensions account for about 81% of variance.
fviz_eig(combined_pca, addlabels = TRUE, ncp = 25)

var_c <- get_pca_var(combined_pca)
print(var_c$contrib)
fviz_contrib(combined_pca, choice = "var", axes = 1) # Red dashed line indicates the expected average contribution

corrplot(var_c$contrib, is.corr=FALSE)

fviz_pca_ind(combined_pca, 
             geom.ind = "point",
             #col.ind = df$INJ_MECH,
             col.ind = df$INJ_MECH,
             addEllipses = TRUE,
             legend.title = "Injury Mechanism",
             )

pca_graph_d0 <- fviz_pca_ind(combined_pca, geom.ind = "point", col.ind = df$INJ_MECH)

ggpubr::ggpar(pca_graph_d0,
              title = "Principal Component Analysis",
              subtitle = "PROPPR data set",
              xlab = "PC1", ylab = "PC2",
              legend.title = "Injury Mechanism", legend.position = "top",
              ggtheme = theme_gray()
              )

# K MEANS
fviz_nbclust(df0, kmeans, method = "silhouette")
set.seed(42)
km.res <- kmeans(df0, 2, nstart = 25)
fviz_cluster(km.res, data = df0,
             elipse.type = "convex", palette = "jco", repel = TRUE,
             ggtheme = theme_minimal()
             )
print(head(km.res$cluster))
