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

df0 = read.xlsx("PROPPR_longitudinal_data_dictionary_edm_5.13.20.xlsx", sheet = "timepoint_2")
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
  df0[paste(substr(x, 5, nchar(x)))] <- scale(df0[x])
}


std_cols <- c(227:262) # Original range 227:262
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
  filename = "R/heatmap_inj_mech_t2.png"
)
