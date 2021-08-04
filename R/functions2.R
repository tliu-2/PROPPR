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
library(broom)

pre_process <- function(df0, impute = F) {
  biomarker_cols <- df0[51:93]
  biomarkers <- colnames(biomarker_cols)
  df0[biomarkers] <- lapply(df0[biomarkers], as.numeric)
  removed_list <- list()
  for (x in biomarkers) {
    if (sum(is.na(df0[x])) / nrow(df0[x]) > 0.2) {
      df0[x] <- NULL
      element_to_remove <- x
      biomarkers <- biomarkers[!(biomarkers %in% element_to_remove)]
      removed_list[[x]] <- x
    }
  }
  
  df0 <- remove_na_patients(df0, biomarkers)
  
  df0[biomarkers] <- as.data.frame(lapply(df0[biomarkers], as.numeric))
  if (impute) {
    for (x in biomarkers) {
      df0[x] <- with(df0, impute(df0[x], mean))
    }
  }
  
  for (x in biomarkers) {
    df0[paste(x, "_std", sep = "")] <- scale(df0[x])
  }
  
  
  return(list("df0" = df0, "cols" = biomarkers,"std_cols" = colnames(df0[227:length(colnames(df0))]), "removed" = removed_list))
}

transpose_u <- function(df) {
  df.t <- transpose(df)
  rownames(df.t) <- colnames(df)
  colnames(df.t) <- rownames(df)
  return(df.t)
}

remove_na_patients <- function(df0, cols) {
  biomarkers <- cols
  df <- transpose_u(df0)
  pos <- 1
  while (pos <= length(colnames(df))) {
    if (sum(is.na(df[biomarkers, pos])) / length(cols) > 0.1) {
      df[pos] <- NULL
      pos <- 1
    } else {
      pos <- pos + 1
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
      p <- ggplot(data = df, aes(x = !! rownames(df), y = !! y_var))  + ylim(0, 1) + geom_bar(stat = "identity", fill = "steelblue", na.rm = T) + ggtitle(i) + geom_text(aes(label = !! y_var),  vjust=-0.3, size=2.5)
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
    df <- NULL
  }
  return(list("df" = df.res, "categories" = new.categories))
}

compare_t_test <- function(df.c1, df.c2) {
  res <- list()
  cols <- colnames(df.c1)
  for (i in cols) {
    res[[i]] <- broom::tidy(t.test(df.c1[i], df.c2[i]))
  }
  return(res)
}