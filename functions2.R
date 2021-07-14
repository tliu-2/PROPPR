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

make_graphs  <- function(df.list, categories) {
  res <- list()
  
  pos <- 1
  for (i in categories) {
    df <- df.list[pos]
    
    y_var <- sym(i)
    res[[i]] <- local({
      i <- i
      rn <- rownames(df)
      p <- ggplot(data = df, aes(x = !! rownames(df), y = !! y_var)) + geom_bar(stat = "identity", fill = "steelblue", na.rm = T) + ggtitle(i) + geom_text(aes(label = !! y_var),  vjust=-0.3, size=3)
      print(p)
    })
    pos <- pos + 1
  }
  return(res)
}

compare_perc_occur <- function(df.list, categories) {
  df.res <- data.frame(matrix(nrow = length(df.list), ncol = 0))
  #res <- list()
  new.categories <- list()
  for (i in categories) {
    occurence.list <- list()
    pos <- 1
    for (d in df.list) {
      total <- nrow(d[i])
      c1 <- nrow(d %>% filter(d[[i]] == 1))
      if (c1 / total < 0.01) {
        break
      }
      occurence.list[paste(pos)] <- c1 / total
      pos <- pos + 1
    }
    if (length(occurence.list) < length(df.list)) {
      next
    }
    df <- bind_rows(occurence.list)
    df <- transpose_u(df)
    colnames(df) <- i
    df.res <- df.res %>%
      add_column(df)
    new.categories <- append(new.categories, i)
    #c1.total <- nrow(df0.c1[i])
    #c1.1 <- nrow(df0.c1 %>%
    #                filter(df0.c1[[i]] == 1))
    #c2.total <- nrow(df0.c2[i])
    #c2.1 <- nrow(df0.c2 %>%
    #               filter(df0.c2[[i]] == 1))
    #df <- data.frame(cluster = "1", occurrence = c1.1/c1.total)
    #df2 <- data.frame(cluster = "2", occurrence = c2.1/c2.total)
    #df <- rbind(df, df2)
    #res[[i]] <- local({
    #  i <- i
    #  p <- ggplot(data = df, aes(x = as.numeric(rownames(df)), y = i)) + geom_bar(stat = "identity", fill = "steelblue", na.rm = T) +
    #    geom_text(aes(label=paste(i, "occurence")), vjust=-0.3, size=3.5) + ggtitle(i)
    #  print(p)
    #})
    df <- NULL
  }
  return(list("df" = df.res, "categories" = new.categories))
}