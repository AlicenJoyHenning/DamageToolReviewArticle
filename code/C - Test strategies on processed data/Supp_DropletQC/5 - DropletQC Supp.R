# CONTEXT
# 
# DropletQC with empty droplets included was run by default in the run strategies workflow.
# This script isolates the two DropletQC run (normal damaged run and damaged & empty droplet run) and compares their performance. 

#-------------------------------------------------------------------------------
# PREPARATIONS
#-------------------------------------------------------------------------------

# Load necessary libraries
packages <- c("dplyr", "ggplot2", "tidyr", "PRROC", "patchwork")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    library(pkg)
  }
}

# Set working directory to Zenodo download
setwd("/Users/alicen/Projects/ReviewArticle/Zenodo")


# Load datasets (5 cases)
apoptotic <- read.csv("./C_Test_Strategies/data/benchmark_output/HEK293_apoptotic.csv")
pro_apoptotic <- read.csv("./C_Test_Strategies/data/benchmark_output/HEK293_proapoptotic.csv")
GM18507_dead <- read.csv("./C_Test_Strategies/data/benchmark_output/GM18507_dead.csv")
GM18507_dying <- read.csv("./C_Test_Strategies/data/benchmark_output/GM18507_dying.csv")
PDX <- read.csv("./C_Test_Strategies/data/benchmark_output/PDX_dead.csv")


#-------------------------------------------------------------------------------
# PERFORMANCE METRICS
#-------------------------------------------------------------------------------

# Calculate the performance metrics ----

# List of methods
methods <- c("DropletQC", "DropletQC_empty")

# Corrected precision function
calc_precision_ci <- function(TP_count, FP_count, n_bootstraps = 1000, epsilon = 1e-10) {
  precisions <- numeric(n_bootstraps)
  
  for (i in 1:n_bootstraps) {
    # Bootstrap resampling
    TP_resample <- rbinom(1, TP_count + FP_count, TP_count / (TP_count + FP_count + epsilon))
    FP_resample <- (TP_count + FP_count) - TP_resample
    precisions[i] <- TP_resample / (TP_resample + FP_resample + epsilon)
  }
  precision_mean <- mean(precisions, na.rm = TRUE)
  precision_lower <- quantile(precisions, 0.025, na.rm = TRUE)
  precision_upper <- quantile(precisions, 0.975, na.rm = TRUE)
  return(c(precision_mean, precision_lower, precision_upper))
}

# Corrected TPR (true positive rate) function
calc_tpr_ci <- function(TP_count, FN_count, n_bootstraps = 1000, epsilon = 1e-10) {
  tprs <- numeric(n_bootstraps)
  for (i in 1:n_bootstraps) {
    # Bootstrap resampling
    TP_resample <- rbinom(1, TP_count + FN_count, TP_count / (TP_count + FN_count + epsilon))
    FN_resample <- (TP_count + FN_count) - TP_resample
    tprs[i] <- TP_resample / (TP_resample + FN_resample + epsilon)
  }
  tpr_mean <- mean(tprs, na.rm = TRUE)
  tpr_lower <- quantile(tprs, 0.025, na.rm = TRUE)
  tpr_upper <- quantile(tprs, 0.975, na.rm = TRUE)
  return(c(tpr_mean, tpr_lower, tpr_upper))
}

# Corrected PR-AUC function
calc_pr_auc_ci <- function(scores, labels, n_bootstraps = 1000) {
  aucs <- numeric(n_bootstraps)
  for (i in 1:n_bootstraps) {
    # Bootstrap resampling
    idx <- sample(seq_along(scores), replace = TRUE)
    pr <- pr.curve(scores.class0 = scores[idx], weights.class0 = labels[idx], curve = FALSE)
    aucs[i] <- pr$auc.integral
  }
  auc_mean <- mean(aucs, na.rm = TRUE)
  auc_lower <- quantile(aucs, 0.025, na.rm = TRUE)
  auc_upper <- quantile(aucs, 0.975, na.rm = TRUE)
  return(c(auc_mean, auc_lower, auc_upper))
}

# Loop through datasets and methods
for (item in input_list) {
  
  # Extract components
  df <- item$df
  TP <- item$TP
  TN <- item$TN
  dataset_name <- item$name
  
  # Loop through each method
  for (method in methods) {
    message(method, " ...")
    
    # Define method outcomes
    df$method_outcome <- "na"
    df$method_outcome <- ifelse(df$orig.ident == TP & df[[method]] == "damaged", "true_positive", df$method_outcome)
    df$method_outcome <- ifelse(df$orig.ident == TP & df[[method]] == "cell", "false_negative", df$method_outcome)
    df$method_outcome <- ifelse(df$orig.ident == TN & df[[method]] == "cell", "true_negative", df$method_outcome)
    df$method_outcome <- ifelse(df$orig.ident == TN & df[[method]] == "damaged", "false_positive", df$method_outcome)
    
    # Summarize statistics
    outcome <- df %>% 
      group_by(method_outcome) %>%
      summarise(Count = n(), .groups = 'drop')
    
    # Extract confusion matrix counts 
    TP_count <- outcome %>% filter(method_outcome == "true_positive") %>% pull(Count)
    FP_count <- outcome %>% filter(method_outcome == "false_positive") %>% pull(Count)
    FN_count <- outcome %>% filter(method_outcome == "false_negative") %>% pull(Count)
    TN_count <- outcome %>% filter(method_outcome == "true_negative") %>% pull(Count)
    
    TP_count <- ifelse(length(TP_count) == 0, 0, TP_count)
    FP_count <- ifelse(length(FP_count) == 0, 0, FP_count)
    FN_count <- ifelse(length(FN_count) == 0, 0, FN_count)
    TN_count <- ifelse(length(TN_count) == 0, 0, TN_count)
    
    # Calculate metrics with bootstrapped CIs
    precision_ci <- calc_precision_ci(TP_count, FP_count)
    tpr_ci <- calc_tpr_ci(TP_count, FN_count)
    scores <- c(rep(1, TP_count), rep(0, FN_count), rep(1, FP_count), rep(0, TN_count))
    labels <- c(rep(1, TP_count + FN_count), rep(0, FP_count + TN_count))
    pr_auc_ci <- calc_pr_auc_ci(scores, labels)
    
    # Store results
    results <- rbind(results, 
                     data.frame(dataset = dataset_name, 
                                strategy = method, 
                                precision = precision_ci[1], 
                                tpr = tpr_ci[1], 
                                pr_auc = pr_auc_ci[1], 
                                precision_lower = precision_ci[2], 
                                precision_upper = precision_ci[3],
                                tpr_lower = tpr_ci[2], 
                                tpr_upper = tpr_ci[3],
                                pr_auc_lower = pr_auc_ci[2], 
                                pr_auc_upper = pr_auc_ci[3],
                                stringsAsFactors = FALSE))
  }
}

# Add positive prevalence measures -----

# Before show PR-AUC values, find what the value would be for random guessing 
datasets <- c("apoptotic", "pro_apoptotic", "GM18507_dead", "GM18507_dying", "PDX")


# Find the positive prevalence for each dataset 
positive_prevalence <- list()
positive_prevalence$apoptotic <- table(apoptotic$orig.ident)[1] / (table(apoptotic$orig.ident)[1] + table(apoptotic$orig.ident)[2])
positive_prevalence$pro_apoptotic <- table(pro_apoptotic$orig.ident)[2] / (table(pro_apoptotic$orig.ident)[1] + table(pro_apoptotic$orig.ident)[2])
positive_prevalence$GM18507_dead <- table(GM18507_dead$orig.ident)[2] / (table(GM18507_dead$orig.ident)[2] + table(GM18507_dead$orig.ident)[1])
positive_prevalence$GM18507_dying <- table(GM18507_dying$orig.ident)[2] / (table(GM18507_dying$orig.ident)[2] + table(GM18507_dying$orig.ident)[1])
positive_prevalence$PDX <- table(PDX$orig.ident)[2] / (table(PDX$orig.ident)[2] + table(PDX$orig.ident)[1])

# Convert positive_prevalence to a dataframe for easier plotting
positive_prevalence_df <- data.frame(
  dataset = names(positive_prevalence),
  prevalence = unlist(positive_prevalence)
)

# Add to results 
results <- merge(results, positive_prevalence_df, by = c("dataset"))
View(results)

# Save and view final results
write.csv(results, 
          file = "./C_Test_Strategies/data/groundtruth_performance_metrics/DropletQC_empty_droplet_groundtruth_performance.csv",
          quote = FALSE, 
          row.names = FALSE
)


# Plot 



###


#-------------------------------------------------------------------------------
# PLOT DISTRIBUTION
#-------------------------------------------------------------------------------

# Histogram typical to DropletQC showing nuclear fraction distribution 
plot_histogram <- function(df, project_name){
  
  # Define colour palatte 
  colours <- c(empty_droplet = "blue", "damaged" = "#001E5B", cell = "#F1F1F1")
  
  # Those with traditional nuclear fraction
  if (project_name %in% c("apoptotic", "pro_apoptotic")){
  
  plot <- ggplot(df, aes(x = nf, fill = DropletQC_empty)) +
    geom_histogram(binwidth = 0.01, color = "black", alpha = 0.7, position = "identity") +
    scale_fill_manual(values = colours) +
    labs(x = "Nuclear Fraction", y = "Cells", fill = "Category") +
    theme_classic() +
    theme(axis.title = element_text(face = "bold"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8))
  
  } else {
    
    # Use MALAT1 proxy for nuclear fraction 
    
    plot <- ggplot(df, aes(x = malat1, fill = DropletQC_empty)) +
      geom_histogram(binwidth = 20, color = "black", alpha = 0.7, position = "identity") +
      scale_fill_manual(values = colours) +
      labs(x = "Nuclear Fraction", y = "Cells", fill = "Category") +
      theme_classic() +
      theme(axis.title = element_text(face = "bold"),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8))
    
    
  }
  
  # Save 
  location <- paste0("./C_Test_Strategies/img/", project_name, "_DropletQC_empty_droplets.png")
  ggsave(filename = location,
         plot = plot,
         width = 10, height = 3, units = "in")
  
  return(plot)
  
}

apoptotic_plot <- plot_histogram(apoptotic, "apoptotic_test")
pro_apoptotic_plot <- plot_histogram(pro_apoptotic, "pro_apoptotic")
GM18507_dead_plot <- plot_histogram(GM18507_dead, "GM18507_dead")
GM18507_dying_plot <- plot_histogram(GM18507_dying, "GM18507_dying")
PDX_plot <- plot_histogram(PDX, "PDX")


### End 
