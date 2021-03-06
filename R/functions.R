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

get_clust <- function(df, s, c) {
  df.res <- df %>%
    filter(cluster == c) %>%
    filter(`_30DAYST` == s)
  return(df.res)
}

heatmap_process <- function(df, s, c, tpose = TRUE) {
  biomarkers <- c(3:38)
  df.t <- df %>%
    filter(cluster == c) %>%
    filter(Survival == s) %>%
    ungroup(cluster, Survival) %>%
    select(all_of(biomarkers))
  if (tpose) {
    df.t.f <- transpose_u(df.t)
    colnames(df.t.f) <- paste("Cluster", c, s)
    return(df.t.f)
  }
  rownames(df.t) <- paste("Cluster", c, s)
  return(df.t)
}

transpose_u <- function(df) {
  df.t <- transpose(df)
  rownames(df.t) <- colnames(df)
  colnames(df.t) <- rownames(df)
  return(df.t)
}

pre_process2 <- function(df0) {
  biomarker_cols <- df0[51:93]
  biomarkers = colnames(biomarker_cols)
  df0[biomarkers] <- lapply(df0[biomarkers], as.numeric)
  for (x in biomarkers) {
   # df0[x] <- with(df0, impute(df0[x], mean))
    if (sum(is.na(df0[x])) / nrow(df0[x] < 0.2)) {
      df0[paste(x, "_std", sep = "")] <- scale(df0[x])
    }
  }
  return(df0)
}

remove_na_patients <- function(df0) {
  for (x in colnames(df0)) {
    if (sum(is.na(df0[x])) / nrow(df0[x]) > 0.1) {
      df0[x] <- NULL
    }
  }
  return(df0)
}

remove_na_patients2 <- function(df0) {
  biomarker_cols <- df0[51:93]
  biomarkers <- colnames(biomarker_cols)
  df <- transpose_u(df0)
  for (i in colnames(df)) {
    if (sum(is.na(df[biomarkers, i])) / length(biomarkers) > 0.1) {
      df[i] <- NULL
    }
  }
  return(transpose_u(df))
}


# Does a general preprocess for the passed dataframe. Standardizes and imputes
# values. Only used in some cases where specific operations are not needed.
# Takes in a dataframe as a parameter
# Returns the modified dataframe.
pre_process <- function(df0) {
  # Pre-processing
  biomarker_cols <- df0[51:93]
  biomarkers = colnames(biomarker_cols)
  df0[biomarkers] <- lapply(df0[biomarkers], as.numeric)
  # Remove columns with more than 20% values missing
  for (x in biomarkers) {
    if (sum(is.na(df0[x])) / nrow(df0[x]) > 0.2) {
      df0[x] <- NULL
      biomarker_cols[x] <- NULL
    }
  }
  
  biomarkers <- colnames(biomarker_cols)
  # Impute missing values.
  for (x in biomarkers) {
    df0[x] <- with(df0, impute(df0[x], mean))
  }
  
  # Standardize data.
  for (x in biomarkers) {
    df0[paste(substr(x, 5, nchar(x)))] <- scale(df0[x])
  }
  
  std_cols <- c(227:262)
  std_biomarkers <- colnames(df0[std_cols])
  
  
  df0 <- df0 %>%
    filter(INJ_MECH != "Both Types of Injury") %>%
    select(all_of(std_biomarkers))
  
  return(df0)
}

# Function that creates variables that hold counts for type of injury, 30 day
# survival, and complications. Perform chi square on living and dead.
# Takes in a dataframe of patient data.
# Returns a list of dataframes containing the number of patients with
# complications, death counts, and gender.
count_data <- function(df0) {
  df0 <- df0 %>%
    filter(INJ_MECH != "Both Types of Injury")
  
  df_inj_type <- df0 %>%
    select(INJ_MECH) %>%
    count(INJ_MECH)
  
  df_res_30_day <- df0 %>%
    select(INJ_MECH, "_30DAYST") %>%
    count(INJ_MECH, df0$`_30DAYST`)
  colnames(df_res_30_day) <- c("INJ MECH", "Survival_30", "n")
  
  df_res_30_day_chi<- df_res_30_day %>%
    filter(Survival_30 != "Unknown/Lost to Follow-Up")
  
  dead <- df_res_30_day_chi %>%
    filter(Survival_30 == "Deceased") %>%
    select(n)
  rownames(dead) <- c("Blunt", "Pen")
  colnames(dead) <- "Deceased"
  
  live <- df_res_30_day_chi %>%
    filter(Survival_30 == "Living") %>%
    select(n)
  
  rownames(live) <- c("Blunt", "Pen")
  colnames(live) <- "Living"
  
  df_30_day <- cbind(live, dead)
  
  chisq <- chisq.test(df_30_day, correct = FALSE)
  
  # 0 = no, 1 = yes for sepsis, ards, sirs
  df_res_sepsis <- df0 %>%
    select(INJ_MECH, Sepsis)%>%
    count(INJ_MECH, Sepsis)
  
  df_res_death <- df0 %>%
    select(INJ_MECH, Death) %>%
    count(INJ_MECH, Death)
  
  df_res_ards <- df0 %>%
    select(INJ_MECH, ARDS) %>%
    count(INJ_MECH, ARDS)
  
  df_res_sirs <- df0 %>%
    select(INJ_MECH, SIRS) %>%
    count(INJ_MECH, SIRS)
  
  # 0 = male, 1 = female
  df_male_female <- df0 %>%
    select(female1) %>%
    count(female1)
  
  df_gender_inj <- df0 %>%
    select(INJ_MECH, female1) %>%
    count(INJ_MECH, female1)
  
  df_demographics <- df0 %>%
    select(AGE, female1) %>%
    group_by(female1) %>%
    summarise(Mean_Age = mean(AGE, na.rm = TRUE))
  
  
  return(list("chi_sq" = chisq, "num_sepsis" = df_res_sepsis, 
              "num_death" = df_res_death, "num_ards" = df_res_ards, 
              "num_res_sirs" = df_res_sirs, "num_gender" = df_male_female,
              "num_gender_inj" = df_gender_inj, 
              "mean_age_gender" = df_demographics, "num_inj_type" = df_inj_type))
  
}

get_table1_data <- function(df0) {
  df0 <- df0 %>%
    filter(INJ_MECH != "Both Types of Injury")
  
  df_age_all <- df0 %>%
    select(AGE)
  
  df_age_mean <- mean(df_age_all$AGE, na.rm=T)
  df_age_sd <- sd(df_age_all$AGE, na.rm=T)
  
  df_blunt_age <- df0 %>%
    filter(INJ_MECH == "Blunt Injury Only") %>%
    select(AGE)
  
  df_blunt_age_mean <- mean(df_blunt_age$AGE, na.rm=T)
  df_blunt_age_sd <- sd(df_blunt_age$AGE, na.rm=T)
  
  df_pen_age <- df0 %>%
    filter(INJ_MECH == "Penetrating Injury Only") %>%
    select(AGE)
  
  df_pen_age_mean <- mean(df_pen_age$AGE, na.rm=T)
  df_pen_age_sd <- sd(df_pen_age$AGE, na.rm=T)
  
  
  df_num_sex <- df0 %>%
    group_by(INJ_MECH) %>%
    count(female1)
  
  return(list("all_age_mean" = df_age_mean, "all_age_sd" = df_age_sd,
              "blunt_age_mean" = df_blunt_age_mean, "blunt_age_sd" = df_blunt_age_sd,
              "pen_age_mean" = df_pen_age_mean, "pen_age_sd" = df_pen_age_sd,
              "num_sex" = df_num_sex))
    
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
      p <- ggplot(data = df, aes(x = !! rownames(df), y = !! y_var)) + geom_bar(stat = "identity", fill = "steelblue", na.rm = T) + ggtitle(i) + geom_text(aes(label = !! y_var),  vjust=-0.3, size=3.5)
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
    if (length(occurence.list) == 0) {
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

compare_perc_occur_r <- function(df.list, categories, df0) {
  for (i in categories) {
    res.list <- get_num_ones(df.list, i)
  }

}

get_num_ones <- function(df.list, category, df0) {
  res.list <- list()
  for (i in seq_along(1:10)) {
      total <- nrow(df0[a])
      nam <- paste("c.", i, ".", category, sep="")
      assign(nam, nrow(df.list[[i]] %>%
                         filter(df[[category]] == 1)))
  }
  return(res.list)
}
# Function that takes two copies of the dataframe that a pca plot will be 
# Generated from. Parameters do not need to be pre-processed. Function will
# standardize and impute values.
# Takes in a dataframe of patient data.
# Returns a list containing the various graphs generated.
pca_plot <- function(df0) {
  df <- df0
  # Pre-processing
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
  
  # Impute missing values.
  for (x in biomarkers) {
    df0[x] <- with(df0, impute(df0[x], mean))
  }
  
  # Standardize data.
  for (x in biomarkers) {
    df0[paste(substr(x, 5, nchar(x)))] <- scale(df0[x])
  }
  
  std_cols <- c(227:262)
  std_biomarkers <- colnames(df0[std_cols])
  
  
  df0 <- df0 %>%
    filter(INJ_MECH != "Both Types of Injury") %>%
    select(all_of(std_biomarkers))
  
  df <- df %>%
    filter(INJ_MECH != "Both Types of Injury")
  
  # Make PCA plot.
  combined_pca = PCA(df0, graph = FALSE, ncp = 25)
  combined_eig_val <- get_eigenvalue(combined_pca)
  print(combined_eig_val) # 15 dimensions account for about 81% of variance.
  scree <- fviz_eig(combined_pca, addlabels = TRUE, ncp = 25)
  
  var_c <- get_pca_var(combined_pca)
  print(var_c$contrib)
  # Red dashed line indicates the expected average contribution
  var_contrib <- fviz_contrib(combined_pca, choice = "var", axes = 1:2)
  
  corr_plot <- corrplot(var_c$contrib, is.corr=FALSE)
  
  plot <- fviz_pca_ind(combined_pca, 
               geom.ind = "point",
               col.ind = df$INJ_MECH,
               addEllipses = TRUE,
               legend.title = "Injury Mechansim",
  )
  
  pca_graph_d0 <- fviz_pca_ind(combined_pca, geom.ind = "point", 
                               col.ind = df$INJ_MECH)
  
  plot_pubr <- ggpubr::ggpar(pca_graph_d0,
                title = "Principal Component Analysis",
                subtitle = "PROPPR data set",
                xlab = "PC1", ylab = "PC2",
                legend.title = "Injury Mechanism", legend.position = "top",
                ggtheme = theme_gray()
  )
  return(list("scree_plot" = scree, "var_contrib" = var_contrib, 
              "plot" = plot_pubr, "plot_ellipse" = plot))
}

# Function that does k-means clustering. Takes in a dataframe that should be
# already preprocessed.
# Returns a list containing the k-means plot and silhouette scores plot.
k_means <- function(df0){
  
  df0 <- pre_process(df0)
  # K MEANS
  sil_plot <- fviz_nbclust(df0, kmeans, method = "silhouette")
  set.seed(42)
  km.res <- kmeans(df0, 2, nstart = 25)
  kmeans_plot <- fviz_cluster(km.res, data = df0,
               elipse.type = "convex", palette = "jco", repel = TRUE,
               ggtheme = theme_minimal()
  )
  return(list("kmeans_plot" = kmeans_plot, "silhouette_plot" = sil_plot))
}


# Function that takes a dataframe and generates heatmaps based on injury
# mechanism and death. Standardizes and imputes data as part of the function.
# Saves heatmaps to disk.
gen_heatmaps <- function(df0) {
  
  # Pre-processing
  biomarker_cols <- df0[51:93]
  biomarkers = colnames(biomarker_cols)
  
  # Remove columns with more than 20% values missing
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
  
  # Standardize data.
  for (x in biomarkers) {
    df0[paste(substr(x, 5, nchar(x)))] <- scale(df0[x])
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
  
  # Compare by injury mechanism
  
  map <- pheatmap(
    all_mean_t,
    cluster_rows = FALSE, cluster_cols = FALSE,
    cellwidth = 10,
    cellheight = 10,
    fontsize = 10,
    filename = "R/heatmap_inj_mech.png"
  )
  
  map <- pheatmap(
    all_mean_t,
    cluster_rows = TRUE, cluster_cols = FALSE,
    cellwidth = 10,
    cellheight = 10,
    fontsize = 10,
    filename = "R/heatmap_inj_mech_cluster.png"
  )
  
  # ARDS, Sepsis, SIRS
  
  df_1_b <- df0 %>%
    filter(INJ_MECH == "Blunt Injury Only") %>%
    filter(Death == 1) %>%
    select(all_of(std_biomarkers))
  
  df_1_p <- df0 %>%
    filter(INJ_MECH == "Penetrating Injury Only") %>%
    filter(Death == 1) %>%
    select(all_of(std_biomarkers))
  
  
  df_0_b <- df0 %>%
    filter(INJ_MECH == "Blunt Injury Only") %>%
    filter(Death == 0) %>%
    select(all_of(std_biomarkers))
  
  df_0_p <- df0 %>%
    filter(INJ_MECH == "Penetrating Injury Only") %>%
    filter(Death == 0) %>%
    select(all_of(std_biomarkers))
  
  df_1_mean_b <- df_1_b %>%
    summarise_all(mean, na.rm = TRUE)
  df_1_mean_p <- df_1_p %>%
    summarise_all(mean, na.rm = TRUE)
  df_0_mean_b <- df_0_b %>%
    summarise_all(mean, na.rm = TRUE)
  df_0_mean_p <- df_0_p %>%
    summarise_all(mean, na.rm = TRUE)
  
  df_1_all <- rbind(df_1_mean_b, df_1_mean_p)
  rownames(df_1_all) <- c("Blunt Death", "Pen Death")
  df_0_all <- rbind(df_0_mean_b, df_0_mean_p)
  rownames(df_0_all) <- c("Blunt No Death", "Pen No Death")
  
  all_mean <- rbind(df_1_all, df_0_all)
  all_mean_t <- transpose(all_mean)
  rownames(all_mean_t) <- colnames(all_mean)
  colnames(all_mean_t) <- c("Blunt Death", "Pen Death", "Blunt No Death", "Pen No Death")
  
  map2 <- pheatmap(
    all_mean_t,
    cluster_rows = FALSE, cluster_cols = FALSE,
    cellwidth = 10,
    cellheight = 10,
    fontsize = 10,
    filename = "R/heatmap_death_v2.png"
  )
  
  # Compare living and death
  
  df_dead <- df0 %>%
    filter(Death == 1) %>%
    select(all_of(std_biomarkers))
  
  df_live <- df0 %>%
    filter(Death == 0) %>%
    select(all_of(std_biomarkers))
  
  df_dead_mean <- df_dead %>%
    summarise_all(mean, na.rm = TRUE)
  df_live_mean <- df_live %>%
    summarise_all(mean, na.rm = TRUE)
  
  all_mean_dl <- rbind(df_live_mean, df_dead_mean)
  all_mean__dl_t <- transpose(all_mean_dl)
  rownames(all_mean__dl_t) <- colnames(all_mean_dl)
  colnames(all_mean__dl_t) <- c("Living", "Deceased")
  
  map2 <- pheatmap(
    all_mean__dl_t,
    cluster_rows = FALSE, cluster_cols = FALSE,
    cellwidth = 10,
    cellheight = 10,
    fontsize = 10,
    filename = "R/heatmap_death_vs_living.png"
  )
  
}

# Function that conducts permanova test on the dataset. Returns the results of
# the betadisper function.
permanova <- function(df0) {
  # PERMANOVA 
  df_c_dist <- dist(df0[std_biomarkers], method = "euclidean")
  set.seed(42)
  dif_c_div <- adonis2(df_c_dist~INJ_MECH, data=df0, permutations = 999, method="bray")
  print(dif_c_div)
  
  dispersion <- betadisper(df_c_dist, group=df0$INJ_MECH)
  print(dispersion)
  permutest(dispersion)
  plot(dispersion, hull=FALSE, ellipse = TRUE)
  
  return(dispersion)
}

# Divides the population into two separate dataframes based on injury type
# and performs k-means clustering on patients of only one injury type.
# Returns a list containing a heatmaps for blunt and one for penetrating
# injury biomarker expression.
k_means_sep <- function(df0) {
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
  b_cluster <- fviz_cluster(km.res_blunt, data = df_blunt,
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
  cluster_heatmap_b <- pheatmap(
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
  p_cluster <- fviz_cluster(km.res_pen, data = df_pen,
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
  cluster_heatmap_p <- pheatmap(
    df.p.t,
    cluster_rows = FALSE, cluster_cols = FALSE,
    cellwidth = 10,
    cellheight = 10,
    fontsize = 10,
  )
  
  return(list("Blunt_Heat" = cluster_heatmap_b, "Pen_Heat" = cluster_heatmap_p,
              "Blunt_Cluster" = b_cluster, "Pen_Cluster" = p_cluster))
}


hierarch_cluster <- function(df0) {
  
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
  
  df0 <- df0 %>%
    select(all_of(biomarkers))
  
  d <- dist(df0, method="euclidean")
  hcl <- hclust(d, method="complete")

  plot(hcl, cex = 0.6, hang = -1)
  
  hc5 <- hclust(d, method="ward.D2")
  sub_grp <- cutree(hc5, k = 2)
  grouped_plot <- plot(hc5, cex = 0.6)
  rect.hclust(hc5, k = 2, border = 2:5)
  plot_cluster <- recordPlot()
  
  hierarch_cluster_plot <- fviz_cluster(list(data = df0, cluster = sub_grp))
  return(list("dendro" = plot_cluster, "cluster_plot" = hierarch_cluster_plot))
}

hierarch_cluster_sep <- function(df0) {
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
  
  # Blunt
  d_b <- dist(df_blunt, method="euclidean")
  hc1_b <- hclust(d_b, method="complete")
  plot(hc1_b, cex = 0.6, hang = -1)
  # rect.hclust(hc1_b, k = 2, border = 2:5)
  sub_grp_b <- cutree(hc1_b, k = 2)
  agg_cluster_b <- recordPlot()
  
  
  hierarch_cluster_plot_b <- fviz_cluster(list(data = df_blunt, cluster = sub_grp_b))
  
  df_b <- df0 %>%
    filter(INJ_MECH == "Blunt Injury Only") %>%
    mutate(cluster = sub_grp_b)
  
  # Penetrating
  d_p <- dist(df_pen, method="euclidean")
  hc1_p <- hclust(d_p, method="complete")
  plot(hc1_p, cex = 0.6, hang = -1)
  # rect.hclust(hc1_p, k = 2, border = 2:5)
  sub_grp_p <- cutree(hc1_p, k = 2)
  agg_cluster_p <- recordPlot()
  
  #rect.hclust(hc5_p, k = 2, border = 2:5)
  
  hierarch_cluster_plot_p <- fviz_cluster(list(data = df_pen, cluster = sub_grp_p))
  
  df_p <- df0 %>%
    filter(INJ_MECH == "Penetrating Injury Only") %>%
    mutate(cluster = sub_grp_p)
  
  return(list("blunt_plot" = hierarch_cluster_plot_b, 
              "pen_plot" = hierarch_cluster_plot_p, "df_b" = df_b, "df_p" = df_p,
              "dendro_b_agg" = agg_cluster_b, "dendro_p_agg" = agg_cluster_p))
  
}

cluster_heatmap <- function(df_b, df_p) {
 
  # Pre-processing
  biomarker_cols <- df0[51:93] #df0[8:50] 
  # biomarker_cols = df0[8:50]
  biomarkers = colnames(biomarker_cols)
  
  for (x in biomarkers) {
    if (sum(is.na(df0[x])) / nrow(df0[x]) > 0.2) {
      df_b[x] <- NULL
      df_p[x] <- NULL
    }
    if(sum(is.na(biomarker_cols[x])) / nrow(biomarker_cols[x]) > 0.2) {
      biomarker_cols[x] <- NULL
    }
  }
  
  biomarkers <- colnames(biomarker_cols)
  
  # Impute missing values.
  for (x in biomarkers) {
    df_b[x] <- with(df_b, impute(df_b[x], mean))
    df_p[x] <- with(df_p, impute(df_p[x], mean))
  } 
  
  # Standardize data.
  for (x in biomarkers) {
    df_b[paste(substr(x, 5, nchar(x)))] <- scale(df_b[x])
    df_p[paste(substr(x, 5, nchar(x)))] <- scale(df_p[x])
  }
  
  std_cols <- c(228:263)
  std_biomarkers <- colnames(df_b[std_cols])
  
  df_b_1 <- df_b %>%
    filter(cluster == 1) %>%
    select(all_of(std_biomarkers)) %>%
    summarise_all(mean, na.rm=T)
  df_b_2 <- df_b %>%
    filter(cluster == 2) %>%
    select(all_of(std_biomarkers)) %>%
    summarise_all(mean, na.rm=T)
  
  df_p_1 <- df_p %>%
    filter(cluster == 1) %>%
    select(all_of(std_biomarkers)) %>%
    summarise_all(mean, na.rm=T)
  df_p_2 <- df_p %>%
    filter(cluster == 2) %>%
    select(all_of(std_biomarkers)) %>%
    summarise_all(mean, na.rm=T)
  
  df_b_all <- rbind(df_b_1, df_b_2)
  df_b_all_t <- transpose(df_b_all)
  rownames(df_b_all_t) <- colnames(df_b_all)
  colnames(df_b_all_t) <- rownames(df_b_all)
  
  df_p_all <- rbind(df_p_1, df_p_2)
  df_p_all_t <- transpose(df_p_all)
  rownames(df_p_all_t) <- colnames(df_p_all)
  colnames(df_p_all_t) <- rownames(df_p_all)
  
  map_b <- pheatmap(
    df_b_all_t,
    cluster_rows = FALSE, cluster_cols = FALSE,
    cellwidth = 10,
    cellheight = 10,
    fontsize = 10
  )
  
  map_p <- pheatmap(
    df_p_all_t,
    cluster_rows = FALSE, cluster_cols = FALSE,
    cellwidth = 10,
    cellheight = 10,
    fontsize = 10
  )
  
  return(list("mapb" = map_b, "mapp" = map_p)) 
  
}