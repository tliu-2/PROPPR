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
library("logbin")

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

for (x in biomarkers) {
  df0[x] <- with(df0, impute(df0[x], mean))
}

for (x in biomarkers) {
  df0[paste(substr(x, 5, nchar(x)))] <- scale(df0[x])
}

std_cols <- c(227:262)
std_biomarkers <- colnames(df0[std_cols])

df_std <- df0[std_biomarkers]
test_lm <- lm(Death ~ `_Hu_IL_6__19`, data = df0)
print(test_lm)
summary(test_lm)
