# SCRIPT CONTEXT 
#
# Benchmarking damaged cell detection strategies on ground truth dataset. This script begins with the strategies
# already run on the ground truth cases. The labelled output objects (3 - Run remaining detection strategies.R)
# are loaded at the start of the script. From here, the script uses the labels stored in the objects to 
# calculate confusion metrics (TP, TN, FP, FN) for each strategy, followed by sensitivity and F1 scores.
# The script then plots the sensitivity and F1 scores for each strategy
# 
# Still do : 
# - automated rather than manual confusion metric caculations 
# - automated precision, recall, PR-AUC calculations 
# - finalised plots (new aesthetics, labels, etc)


#-------------------------------------------------------------------------------
# PREPARATIONS
#-------------------------------------------------------------------------------

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(PRROC)
library(patchwork)  # For arranging plots

# Load datasets (5 cases)
apoptotic <- read.csv("/home/alicen/Projects/ReviewArticle/benchmark_results/HEK293_apo.csv")
pro_apoptotic <- read.csv("/home/alicen/Projects/ReviewArticle/benchmark_results/HEK293_pro.csv")
GM18507_dead <- read.csv("/home/alicen/Projects/ReviewArticle/benchmark_results/GM18507_dead.csv")
GM18507_dying <- read.csv("/home/alicen/Projects/ReviewArticle/benchmark_results/GM18507_dying.csv")
PDX <- read.csv("/home/alicen/Projects/ReviewArticle/benchmark_results/PDX_dead.csv")

#-------------------------------------------------------------------------------
# PERFORMANCE METRICS
#-------------------------------------------------------------------------------

# List of methods
methods <- c("ddqc", "DropletQC", "ensembleKQC", "miQC", "valiDrops", 
             "manual_all", "manual_mito_ribo", "manual_mito", "manual_malat1", "manual_mito_isolated")

# Input structure
input_list <- list(
  list(name = "apoptotic", df = apoptotic, TP = "HEK293_apoptotic", TN = "HEK293_control"),
  list(name = "proapoptotic", df = pro_apoptotic, TP = "HEK293_proapoptotic", TN = "HEK293_control"),
  list(name = "GM18507_dead", df = GM18507_dead, TP = "GM18507_dead", TN = "GM18507_control"),
  list(name = "GM18507_dying", df = GM18507_dying, TP = "GM18507_dying", TN = "GM18507_control"),
  list(name = "PDX", df = PDX, TP = "PDX_dead", TN = "PDX_control")
)

# Initialize df for storage
results <- data.frame(dataset = character(), 
                      strategy = character(), 
                      precision = numeric(), 
                      fnr = numeric(), 
                      pr_auc = numeric(),
                      precision_lower = numeric(),
                      precision_upper = numeric(),
                      fnr_lower = numeric(),
                      fnr_upper = numeric(),
                      pr_auc_lower = numeric(),
                      pr_auc_upper = numeric(),
                      stringsAsFactors = FALSE)

# Function to calculate 95% CI ∈ (0, 1) (precision and FNR)
calc_ci <- function(p, n) {
  error <- qnorm(0.975) * sqrt((p * (1 - p)) / n)
  lower <- max(p - error, 0)  
  upper <- min(p + error, 1)  
  return(c(lower, upper))
}

# Function to calculate PR-AUC with bootstrapped CIs (1000 bootstraps)
calc_pr_auc_ci <- function(scores, labels, n_bootstraps = 1000) {
  aucs <- numeric(n_bootstraps)
  for (i in 1:n_bootstraps) {
    # Bootstrap resampling
    idx <- sample(seq_along(scores), replace = TRUE)
    pr <- pr.curve(scores.class0 = scores[idx], weights.class0 = labels[idx], curve = FALSE)
    aucs[i] <- pr$auc.integral
  }
  auc_mean <- mean(aucs)
  auc_lower <- quantile(aucs, 0.025)
  auc_upper <- quantile(aucs, 0.975)
  return(c(auc_mean, auc_lower, auc_upper))
}

# Calculate precision, FNR, and PR-AUC for each strategy for each dataset
for (item in input_list) {
  
  # Extract components of the list 
  df <- item$df
  TP <- item$TP
  TN <- item$TN
  dataset_name <- item$name
  
  # Loop through each method
  for (method in methods) {
    
    df$method_outcome <- "na"
    
    # Assign method_outcome based on conditions
    df$method_outcome <- ifelse(df$orig.ident == TP & df[[method]] == "damaged", "true_positive", df$method_outcome)
    df$method_outcome <- ifelse(df$orig.ident == TP & df[[method]] == "cell", "false_negative", df$method_outcome)
    df$method_outcome <- ifelse(df$orig.ident == TN & df[[method]] == "cell", "true_negative", df$method_outcome)
    df$method_outcome <- ifelse(df$orig.ident == TN & df[[method]] == "damaged", "false_positive", df$method_outcome)
    
    # Summarise for stats
    outcome <- df %>% 
      group_by(method_outcome) %>%
      summarise(Count = n(), .groups = 'drop')
    
    # Extract confusion matrix counts 
    TP_count <- outcome %>% dplyr::filter(method_outcome == "true_positive") %>% pull(Count)
    FP_count <- outcome %>% dplyr::filter(method_outcome == "false_positive") %>% pull(Count)
    FN_count <- outcome %>% dplyr::filter(method_outcome == "false_negative") %>% pull(Count)
    TN_count <- outcome %>% dplyr::filter(method_outcome == "true_negative") %>% pull(Count)
    
    TP_count <- ifelse(length(TP_count) == 0, 0, TP_count)
    FP_count <- ifelse(length(FP_count) == 0, 0, FP_count)
    FN_count <- ifelse(length(FN_count) == 0, 0, FN_count)
    TN_count <- ifelse(length(TN_count) == 0, 0, TN_count)
    
    # Calculate precision and FNR
    precision <- ifelse((TP_count + FP_count) > 0, TP_count / (TP_count + FP_count), NA)
    fnr <- ifelse((TP_count + FN_count) > 0, FN_count / (TP_count + FN_count), NA)
    
    # Calculate confidence intervals
    precision_ci <- ifelse(!is.na(precision), calc_ci(precision, TP_count + FP_count), c(NA, NA))
    fnr_ci <- ifelse(!is.na(fnr), calc_ci(fnr, TP_count + FN_count), c(NA, NA))
    
    # Calculate PR-AUC and its confidence interval using bootstrapping
    scores <- c(rep(1, TP_count), rep(0, FN_count), rep(1, FP_count), rep(0, TN_count))
    labels <- c(rep(1, TP_count + FN_count), rep(0, FP_count + TN_count))
    pr_auc_ci <- calc_pr_auc_ci(scores, labels)
    
    # Store results
    results <- rbind(results, 
                     data.frame(dataset = dataset_name, strategy = method, precision = precision, fnr = fnr, 
                                pr_auc = pr_auc_ci[1], 
                                precision_lower = precision_ci[1], precision_upper = precision_ci[2],
                                fnr_lower = fnr_ci[1], fnr_upper = fnr_ci[2],
                                pr_auc_lower = pr_auc_ci[2], pr_auc_upper = pr_auc_ci[3],
                                stringsAsFactors = FALSE))
  }
}

# View final results
View(results)

#-------------------------------------------------------------------------------
# PLOTTING
#-------------------------------------------------------------------------------

# Convert the strategy column to a factor with the specified order
results$strategy <- factor(results$strategy, levels = methods)

# Pivot the data frame to long format
df_long <- results %>%
  pivot_longer(cols = c(precision, fnr), names_to = "measure", values_to = "value")

# Define the colors for the strategies
strategy_colors <- c("#999FA9", "#586667", "#73966E", "#415741", "#65719E", "#383F54", "#E8E5D6", "#A4A388", "#8B3860", "#522239")

# Create the bar plot for precision
precision_plot <- ggplot(df_long %>% dplyr::filter(measure == "precision"), aes(x = strategy, y = value, fill = strategy)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ dataset, nrow = 1) +
  scale_fill_manual(values = strategy_colors) +
  labs(title = "Precision", y = "Precision") + 
  theme_classic() +  
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = "none",
    legend.justification = "center",
    axis.text = element_text(size = 16),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    strip.background = element_rect(fill = "#8393C6"),
    strip.text = element_text(color = "white", size = 12, face = "bold")
  )

# Create the bar plot for FNR
fnr_plot <- ggplot(df_long %>% dplyr::filter(measure == "fnr"), aes(x = strategy, y = value, fill = strategy)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ dataset, nrow = 1) +
  scale_fill_manual(values = strategy_colors) +
  labs(title = "False Negative Rate (FNR)", y = "FNR") + 
  theme_classic() +  
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = "none",
    legend.justification = "center",
    axis.text = element_text(size = 16),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    strip.background = element_rect(fill = "#8393C6"),
    strip.text = element_text(color = "white", size = 12, face = "bold")
  )

# Create the plot for PR-AUC
pr_auc_plot <- ggplot(results, aes(x = strategy, y = pr_auc, fill = strategy)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ dataset, nrow = 1) +
  scale_fill_manual(values = strategy_colors) +
  labs(title = "PR-AUC", y = "PR-AUC") + 
  theme_classic() +  
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = "none",
    legend.justification = "center",
    axis.text = element_text(size = 16),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    strip.background = element_rect(fill = "#8393C6"),
    strip.text = element_text(color = "white", size = 12, face = "bold")
  )

# Arrange the plots in a single layout
final_plot <- (precision_plot / fnr_plot / pr_auc_plot)

# Display the plot
final_plot
