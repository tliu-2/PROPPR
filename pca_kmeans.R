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
  cellwidth = 1,
  cellheight = 5,
  fontsize = 5,
  #color = rampcolors,
  #breaks = breaks,
  filename = "R/heatmap_alltest.png" 
)

df0.p <- transpose_u(df0.p.t)
df0.p2 <- df0.p %>%
  select(all_of(std_biomarkers))
df0.p2[std_biomarkers] <- lapply(df0.p2[std_biomarkers], as.numeric)

distance  <- get_dist(df0.p2)
fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

fviz_nbclust(df0.p2, kmeans, method = "silhouette")
fviz_nbclust(df0.p2, kmeans, method = "gap_stat")

k5 <- kmeans(df0.p2, centers = 2, nstart = 25)
fviz_cluster(k5, data = df0.p2)

res.clust <- cbind(df0.p, cluster = k5$cluster)

outcomes <- colnames(df0[173:199])
outcomes <- append(outcomes, c("RBC_10", "RAN_3HRST", "RAN_24HRST"))

df0.c1 <- res.clust %>%
  filter(cluster == 1)
df0.c2 <- res.clust %>%
  filter(cluster == 2)
df0.c3 <- res.clust %>%
  filter(cluster == 3)
df0.c4 <- res.clust %>%
  filter(cluster == 4)
df0.c5 <- res.clust %>%
  filter(cluster == 5)

df0.blunt <- res.clust %>%
  filter(INJ_MECH == "Blunt Injury Only")

df0.pen <- res.clust %>%
  filter(INJ_MECH == "Penetrating Injury Only")

df0.blunt.biomarkers <- df0.blunt %>%
  select(all_of(biomarkers))
df0.pen.biomarkers <- df0.pen %>%
  select(all_of(biomarkers))

df0.c1.biomarkers <- df0.c1 %>%
  select(all_of(biomarkers))
df0.c2.biomarkers <- df0.c2 %>%
  select(all_of(biomarkers))
df0.c3.biomarkers <- df0.c3 %>%
  select(all_of(biomarkers))
df0.c4.biomarkers <- df0.c4 %>%
  select(all_of(biomarkers))
df0.c5.biomarkers <- df0.c5 %>%
  select(all_of(biomarkers))

df0.blunt.biomarkers <- as.data.frame(lapply(df0.blunt.biomarkers, as.numeric))
df0.pen.biomarkers <- as.data.frame(lapply(df0.pen.biomarkers, as.numeric))

df0.c1.biomarkers <- as.data.frame(lapply(df0.c1.biomarkers, as.numeric))
df0.c2.biomarkers <- as.data.frame(lapply(df0.c2.biomarkers, as.numeric))
df0.c3.biomarkers <- as.data.frame(lapply(df0.c3.biomarkers, as.numeric))
df0.c4.biomarkers <- as.data.frame(lapply(df0.c4.biomarkers, as.numeric))
df0.c5.biomarkers <- as.data.frame(lapply(df0.c5.biomarkers, as.numeric))


cluster.list <- list(df0.c1, df0.c2)#, df0.c3, df0.c4, df0.c5)

res.list <- compare_perc_occur(cluster.list, outcomes)

plots <- make_graphs(res.list$df, res.list$categories)

n <- length(plots)
nCol <- floor(sqrt(n))
g <- arrangeGrob(grobs = plots, ncol = nCol)
ml <- marrangeGrob(grobs = plots, nrow = 1, ncol = 2)
ggsave(file="occurrencekmeans2clust.pdf", ml)


# Table 3
# Estimate1 = mean of x, Estimate2 = mean of y, Estimate = mean of x + y
# Statistic = t-statistic
test.res <- compare_t_test(df0.c1.biomarkers, df0.c2.biomarkers)
test.res.all <- data.frame(biomarker = rep(names(test.res)), bind_rows(test.res))
print(test.res.all)
print(test.res$log2_Hu_IL_1b__39)
