library("openxlsx")
library("tidyverse")
#library("FactoMineR")
library("factoextra")
library("corrplot")
library(Hmisc)
#library(missMDA)
#library(vegan)
library(pheatmap)
library(data.table)
#library(cluster)
#library(dendextend)
library(ggplot2)
library(gridExtra)
library(viridis)
library(broom)
library(rgl)
library(ComplexHeatmap)

source("./R/functions2.R")

# Misc code

df.inner.join.1 <- inner_join(df.hclust.1, df.kmeans.1, by = "biomarker_key_SUBJECTID")

df.inner.join.2 <- inner_join(df.hclust.2, df.kmeans.2, by = "biomarker_key_SUBJECTID")

df.left.join.1 <- left_join(df.hclust.1, df.kmeans.1, by = "biomarker_key_SUBJECTID")
df.left.join.2 <- left_join(df.hclust.2, df.kmeans.2, by = "biomarker_key_SUBJECTID")

df.clst.list <- list("Cluster1" = df.left.join.1, "Cluster2" = df.left.join.2)
write.xlsx(df.clst.list, file = "Compare_K_Hclust/compare_clusters_left_join.xlsx")

df.old <- read.csv("dataset_clustering_8.4.21.csv", stringsAsFactors = F)
df.new <- read.csv("dataset_clustering_complex_8.4.21.csv", stringsAsFactors = F)

df.check.1 <- inner_join(df.old %>% filter(cluster == 1), df.new %>% filter(cluster == 2), by = c("biomarker_key_PROPPRID"))
df.check.2 <- inner_join(df.old %>% filter(cluster == 2), df.new %>% filter(cluster == 1), by = c("biomarker_key_PROPPRID"))


# Start main analysis code

df0 = read.xlsx("./data/PROPPR_longitudinal_data_dictionary_edm_5.13.20.xlsx", sheet = "timepoint_0")
df0 <- df0 %>%
  filter(INJ_MECH != "Both Types of Injury")

print(df0 %>% select(all_of(outcomes), INJ_MECH, female1, AGE, WHITE_RE, BLACK_RE, OTHER_RE) 
      %>% group_by(INJ_MECH) %>% 
        summarise(mean_age = mean(AGE, na.rm = T),  sd_age = sd(AGE, na.rm = T)))

# Preprocess
res_pre <- pre_process(df0, T)
df0.p <- res_pre$df0
std_biomarkers <- res_pre$std_cols
biomarkers <- res_pre$cols
print(res_pre$removed)

df0.p[std_biomarkers] <- lapply(df0.p[std_biomarkers], as.numeric)

df0.p.t <- transpose_u(df0.p)

subj.num <- colnames(df0.p.t)

df0.p.map <- transpose_u(df0.p.t) %>% select(all_of(std_biomarkers))
df0.p.map <- as.data.frame(lapply(df0.p.map, as.numeric))
rownames(df0.p.map) <- subj.num

table_s1 <- sub("log2_", "", std_biomarkers)
table_s1 <- sub("Hu_", "", table_s1)
table_s1 <- sub("_std", "", table_s1)

df0.p.map.ex <- data.frame(row = subj.num, df0.p.map)
df0.p.map.ex <- as.data.frame(lapply(df0.p.map.ex, as.numeric))
write.xlsx(df0.p.map.ex, file = "test.xlsx")

# End preprocess

# PCA on hclust clusters

map.all.h <- pheatmap(
  transpose_u(df0.p.map),
  cluster_rows = F, cluster_cols = T,
  cellwidth = 5,
  cellheight = 5,
  fontsize = 5,
  cutree_cols = 2,
  #filename = './R_figures/heatmap_hclust_blunt_only.png'
)

# PCA Start
fviz_nbclust(df0.p.map, hcut, method = "silhouette")

clusters <- cutree(map.all.h$tree_col, k = 2)
# Do hclust and assign to subjects.
res.clust.h <- cbind(transpose_u(df0.p.t), cluster = cutree(map.all.h$tree_col, k = 2))
res.clust.h[std_biomarkers] <- lapply(res.clust.h[std_biomarkers], as.numeric)

# Do PCA
res.pca.h <- PCA(res.clust.h[std_biomarkers], graph = F)
fviz_screeplot(res.pca.h, addlabels = T, ylim = c(0, 50))
var_pca <- get_pca_var(res.pca.h)
fviz_contrib(res.pca.h, choice = "var", axes = 1, top = 10)
# PCA Plot
fviz_pca_ind(res.pca.h,
             label = "none",
             col.ind = as.character(res.clust.h$cluster),
             addEllipses = F)
# PCA End

df0.c1.h <- res.clust.h %>%
  filter(cluster == 1)
df0.c2.h <- res.clust.h %>%
  filter(cluster == 2)

print(df0.c1.h %>% count(INJ_MECH))
print(df0.c2.h %>% count(INJ_MECH))

# Summary Heatmap Hclust
res.clust.all.mean <- res.clust.h %>%
  group_by(cluster) %>%
  select(all_of(std_biomarkers)) %>%
  summarise_all(mean, na.rm=T) %>%
  ungroup()
res.clust.all.mean$cluster <- NULL

pheatmap(
  transpose_u(res.clust.all.mean),
  cluster_rows = T, cluster_cols = F,
  cellwidth = 20,
  cellheight = 20,
  fontsize = 10,
  display_numbers = T,
  treeheight_row = 0,
  color = rampcolors,
  breaks = breaks,
  filename = "R_figures/hclust_summary.png" 
)


# K Means

df0.p.k <- transpose_u(df0.p.t)

df0.p2 <- df0.p.k %>%
  select(all_of(std_biomarkers))
df0.p2[std_biomarkers] <- lapply(df0.p2[std_biomarkers], as.numeric)
colnames(df0.p2) <- table_s1
fviz_nbclust(df0.p2, kmeans, method = "silhouette")
fviz_nbclust(df0.p2, kmeans, method = "gap_stat")
fviz_nbclust(df0.p2, kmeans, method = "wss")

color <- turbo(50)
breaks.low <- c(seq(-0.6, -0.41, by = 0.01))
rampcolors.low <- colorRampPalette(color[1:10])(length(breaks.low))
breaks.mid <- c(seq(-0.4, 0.4, by = 0.01))
rampcolors.mid <- colorRampPalette(color[10:40])(length(breaks.mid))
breaks.high <- c(seq(0.41, 0.6, by = 0.01))
rampcolors.high <- colorRampPalette(color[40:50])(length(breaks.high))

breaks <- c(breaks.low, breaks.mid, breaks.high)
rampcolors <- c(rampcolors.low, rampcolors.mid, rampcolors.high)
col.maps <- circlize::colorRamp2(breaks, rampcolors)

set.seed(42)

# Old pheatmap, deprecated. 
#map.all.k <- pheatmap(
#  df0.p2,
#  kmeans_k=2,
#  cluster_rows = F, cluster_cols = T,
#  cellwidth = 20,
#  cellheight = 20,
#  fontsize = 10,
#  cutree_cols = 2,
#  display_numbers = T,
#  treeheight_col = 0,
#  color = rampcolors,
#  breaks = breaks,
#  #filename = "./R_figures/kmeans_heatmap.png"
#)

# K-Means with complex heatmap
set.seed(42)
png(filename = "./R_figures/complex_kmeans_heatmap_v2.png", width=10,height=10,units="in",res=1200)
map.all.k <- Heatmap(
  data.matrix(df0.p2),
  cluster_columns = F,
  name = "z_scores",
  cluster_rows = F,
  row_km = 2,
)
ht <- draw(map.all.k)
ht
dev.off()
# clusterlist = row_order(ht)

# Get clusters and assign to subjects.
# clusters.k <- map.all.k$kmeans$cluster

row.names(df0.p.k) <- df0.p.k$biomarker_key_PROPPRID
rcl.list <- row_order(ht)
lapply(rcl.list, function(x) length(x))
for (i in 1:length(row_order(ht))){
  if (i == 1) {
    clu <- t(t(row.names(data.matrix(df0.p.k[row_order(ht)[[i]],]))))
    out <- cbind(clu, paste(i))
    colnames(out) <- c("biomarker_key_PROPPRID", "cluster")
    } else {
      clu <- t(t(row.names(data.matrix(df0.p.k[row_order(ht)[[i]],]))))
      clu <- cbind(clu, paste(i))
      out <- rbind(out, clu)
    }
}

out

res.clust.k <- cbind(df0.p.k, cluster = map.all.k$kmeans$cluster)
res.clust.k <- inner_join(df0.p.k, out, by = "biomarker_key_PROPPRID", copy = T)

res.clust.k[std_biomarkers] <- lapply(res.clust.k[std_biomarkers], as.numeric)

df0.c1.k.mean <- res.clust.k %>%
  filter(cluster == 1) %>%
  select(all_of(std_biomarkers)) %>%
  summarise_all(mean, na.rm = T)

df0.c2.k.mean <- res.clust.k %>%
  filter(cluster == 2) %>%
  select(all_of(std_biomarkers)) %>%
  summarise_all(mean, na.rm = T)


# Create line plot to visualize differences between the two clusters.
df0.all.k.mean <- rbind(df0.c1.k.mean, df0.c2.k.mean)
df0.all.k.mean.t <- transpose_u(df0.all.k.mean)
df0.all.k.mean.t <- data.frame(biomarkers = rownames(df0.all.k.mean.t), value1 = df0.all.k.mean.t$`1`, value2 = df0.all.k.mean.t$`2`)
df0.all.k.mean.t <- df0.all.k.mean.t %>%
  mutate(diff = value2 - value1) %>%
  arrange(-diff)
df0.all.k.mean.t$biomarkers <- factor(df0.all.k.mean.t$biomarkers, levels = df0.all.k.mean.t$biomarkers[order(-df0.all.k.mean.t$diff)])

p = ggplot() + geom_line(data = df0.all.k.mean.t, aes(x = biomarkers, y = value1, color = "Cluster 1", group = 1)) +
  geom_line(data = df0.all.k.mean.t, aes(x = biomarkers, y = value2, color = "Cluster 2", group = 1)) + 
  theme(axis.text.x = element_text(size = 5, angle = 45, hjust = 1)) + 
  labs(x = "biomarkers", y = "mean z-score of log2")
  
print(p)
ggsave("compare_clusters_line2.png", plot = p)


# Do pheatmaps of individuals clusters.
df0.c1.k <- res.clust.k %>% 
  filter(cluster == 1)
df0.c2.k <- res.clust.k %>%
  filter(cluster == 2)

map.clust.k.1 <- res.clust.k %>%
  select(all_of(std_biomarkers), cluster) %>%
  filter(cluster == 1) %>%
  select(all_of(std_biomarkers))

map.clust.k.1[std_biomarkers] <- lapply(map.clust.k.1[std_biomarkers], as.numeric)

map.clust.k.2 <- res.clust.k %>%
  select(all_of(std_biomarkers), cluster) %>%
  filter(cluster == 2) %>%
  select(all_of(std_biomarkers))
map.clust.k.2[std_biomarkers] <- lapply(map.clust.k.2[std_biomarkers], as.numeric)

color <- turbo(50)
breaks.low <- c(seq(-6, -4, by = 0.01))
rampcolors.low <- colorRampPalette(color[1:10])(length(breaks.low))
breaks.mid <- c(seq(-3.9, 3.9, by = 0.01))
rampcolors.mid <- colorRampPalette(color[10:40])(length(breaks.mid))
breaks.high <- c(seq(4, 6, by = 0.01))
rampcolors.high <- colorRampPalette(color[40:50])(length(breaks.high))

breaks <- c(breaks.low, breaks.mid, breaks.high)
rampcolors <- c(rampcolors.low, rampcolors.mid, rampcolors.high)


set.seed(42)

pheatmap(
  transpose_u(map.clust.k.1),
  cluster_rows = T, cluster_cols = T,
  cellwidth = 1,
  cellheight = 20,
  fontsize = 10,
  display_numbers = F,
  color = rampcolors,
  breaks = breaks,
  filename = "./R_figures/kmeans_heatmap_clust1.png"
)
set.seed(42)
pheatmap(
  transpose_u(map.clust.k.2),
  cluster_rows = T, cluster_cols = T,
  cellwidth = 20,
  cellheight = 20,
  fontsize = 10,
  display_numbers = F,
  color = rampcolors,
  breaks = breaks,
  filename = "./R_figures/kmeans_heatmap_clust2.png"
)

# PCA K-means plot
res.clust.pca.k <- res.clust.k %>%
  select(all_of(std_biomarkers))
res.clust.pca.k[std_biomarkers] <- lapply(res.clust.pca.k[std_biomarkers], as.numeric)

res.pca.k <- PCA(res.clust.pca.k, graph = F)

eig.val.k <- get_eigenvalue(res.pca.k)
eig.val.k

fviz_screeplot(res.pca.k, addlabels = T, ylim = c(0, 50))
var_pca.k <- get_pca_var(res.pca.k)

corrplot(var_pca.k$contrib, is.corr=FALSE)  

fviz_contrib(res.pca.k, choice = "var", axes = 1)
fviz_pca_ind(res.pca.k,
             label = "none",
             col.ind = as.character(res.clust.k$cluster),
             addEllipses = F, axes = c(1,2))
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
write.xlsx(df.export.list, file = "./Compare_K_Hclust/compare_clusters.xlsx")

df.export.list.full <- list("Hclust" = res.clust.h, "Kmeans" = res.clust.k)
write.csv(res.clust.k, file = "dataset_clustering_complex_8.4.21.csv")

# Occurence Output

outcomes <- colnames(df0[173:199])
outcomes <- append(outcomes, c("RBC_10", "RAN_3HRST", "RAN_24HRST"))

df0.c1 <- res.clust.k %>%
  filter(cluster == 1)
df0.c2 <- res.clust.k %>%
  filter(cluster == 2)

cluster.list <- list(df0.c1, df0.c2)#, df0.c3, df0.c4)

res.list <- compare_perc_occur(cluster.list, outcomes)

plots <- make_graphs(res.list$df, res.list$categories)

n <- length(plots)
nCol <- floor(sqrt(n))
g <- arrangeGrob(grobs = plots, ncol = nCol)
ml <- marrangeGrob(grobs = plots, nrow = 1, ncol = 2)
ggsave(file="occurence_kmeans.8.4.21.pdf", ml)


# Check against grouping by injury mech
# K-means
df0.blunt.k <- res.clust.k %>%
  filter(INJ_MECH == "Blunt Injury Only")
df0.pen.k <- res.clust.k %>%
  filter(INJ_MECH == "Penetrating Injury Only")


fviz_pca_ind(res.pca.k,
             label = "none",
             col.ind = as.character(res.clust.k$INJ_MECH),
             addEllipses = F)

df0.p.k[std_biomarkers] <- lapply(df0.p.k[std_biomarkers], as.numeric)
res.clust.all.mean.inj <- df0.p.k %>%
  group_by(INJ_MECH) %>%
  select(all_of(std_biomarkers)) %>%
  summarise_all(mean, na.rm=T) %>%
  ungroup()


pheatmap(
  res.clust.all.mean.inj,
  cluster_rows = FALSE, cluster_cols = FALSE,
  cellwidth = 20,
  cellheight = 20,
  fontsize = 10,
  display_numbers = T,
  color = rampcolors,
  breaks = breaks,
  filename = "R_figures/inj_mech_summary.png" 
)

# Attempt to plot on pc 1, 2, and 3
pca.3d.k <- prcomp(res.clust.pca.k)
summary(pca.3d.k)
scores <- as.data.frame(pca.3d.k$x)
plot3d(scores[,1:3], col=res.clust.k$cluster, size=10, type='p', xlim=c(-10, 10), ylim=c(-10,10), zlim=c(-10,10))
#text3d(scores[,1] + 2, scores[,2] + 10, scores[,3] + 2, texts = c(rownames(scores)), cex=0.7, pos=3)