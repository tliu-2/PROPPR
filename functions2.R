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
library(ggplot2)
library(rlang)

pre_process <- function(df0) {
  biomarker_cols <- df0[51:93]
  biomarkers <- colnames(biomarker_cols)
  df0[biomarkers] <- lapply(df0[biomarkers], as.numeric)
  for (x in biomarkers) {
    if (sum(is.na(df0[x])) / nrow(df0[x]) > 0.2) {
      df0[x] <- NULL
      element_to_remove <- x
      biomarkers <- biomarkers[!(biomarkers %in% element_to_remove)]
    }
  }
  
  df0 <- remove_na_patients(df0, biomarkers)
  
  
  df0[biomarkers] <- as.data.frame(lapply(df0[biomarkers], as.numeric))
  for (x in biomarkers) {
    df0[x] <- with(df0, impute(df0[x], mean))
  }
  
  for (x in biomarkers) {
    df0[paste(x, "_std", sep = "")] <- scale(df0[x])
  }
  
  
  return(list("df0" = df0, "cols" = colnames(df0[227:length(colnames(df0))])))
}

transpose_u <- function(df) {
  df.t <- transpose(df)
  rownames(df.t) <- colnames(df)
  colnames(df.t) <- rownames(df)
  return(df.t)
}

remove_na_patients <- function(df0, cols) {
  biomarker_cols <- df0[51:93]
  biomarkers <- colnames(biomarker_cols)
  df <- transpose_u(df0)
  for (i in colnames(df)) {
    if (sum(is.na(df[biomarkers, i])) / length(cols) > 0.1) {
      df[i] <- NULL
    }
  }
  return(transpose_u(df))
}