# SCRIPT CONTEXT 
#
# Benchmarking damaged cell detection strategies on ground truth dataset. This script begins with the strategies
# already run on the ground truth cases. The labelled output objects (3 - Run remaining detection strategies.R)
# are loaded at the start of the script. From here, the script uses the labels stored in the objects to 
# calculate confusion metrics (TP, TN, FP, FN) for each strategy, followed by sensitivity and F1 scores.
# The script then plots the sensitivity and F1 scores for each strategy
# 
# Still do : 
# - confidence intervals for precision and false negative rate 
# - finalize plots (new aesthetics, labels, etc)
# - ranking plot : ranking calculation for each performance metric 
#         - 10 for best, 9 second best ... sum across the datasets


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

# Calculate the performance metrics ----

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


# Define a small epsilon value to avoid division by zero
epsilon <- 1e-10

# Function to calculate precision with bootstrapped CIs (1000 bootstraps)
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

# Function to calculate FNR with bootstrapped CIs (1000 bootstraps)
calc_fnr_ci <- function(TP_count, FN_count, n_bootstraps = 1000, epsilon = 1e-10) {
  fnrs <- numeric(n_bootstraps)
  for (i in 1:n_bootstraps) {
    # Bootstrap resampling
    TP_resample <- rbinom(1, TP_count + FN_count, TP_count / (TP_count + FN_count + epsilon))
    FN_resample <- (TP_count + FN_count) - TP_resample
    fnrs[i] <- FN_resample / (TP_resample + FN_resample + epsilon)
  }
  fnr_mean <- mean(fnrs, na.rm = TRUE)
  fnr_lower <- quantile(fnrs, 0.025, na.rm = TRUE)
  fnr_upper <- quantile(fnrs, 0.975, na.rm = TRUE)
  return(c(fnr_mean, fnr_lower, fnr_upper))
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
  auc_mean <- mean(aucs, na.rm = TRUE)
  auc_lower <- quantile(aucs, 0.025, na.rm = TRUE)
  auc_upper <- quantile(aucs, 0.975, na.rm = TRUE)
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
    
    # Calculate precision and FNR with bootstrapped CIs
    precision_ci <- calc_precision_ci(TP_count, FP_count)
    fnr_ci <- calc_fnr_ci(TP_count, FN_count)
    
    # Calculate PR-AUC and its confidence interval using bootstrapping
    scores <- c(rep(1, TP_count), rep(0, FN_count), rep(1, FP_count), rep(0, TN_count))
    labels <- c(rep(1, TP_count + FN_count), rep(0, FP_count + TN_count))
    pr_auc_ci <- calc_pr_auc_ci(scores, labels)
    
    # Store results
    results <- rbind(results, 
                     data.frame(dataset = dataset_name, strategy = method, 
                                precision = precision_ci[1], fnr = fnr_ci[1], pr_auc = pr_auc_ci[1], 
                                precision_lower = precision_ci[2], precision_upper = precision_ci[3],
                                fnr_lower = fnr_ci[2], fnr_upper = fnr_ci[3],
                                pr_auc_lower = pr_auc_ci[2], pr_auc_upper = pr_auc_ci[3],
                                stringsAsFactors = FALSE))
  }
}

# Save and view final results
write.csv(results, 
          file = "/home/alicen/Projects/ReviewArticle/benchmark_results/performance_metrics/groundtruth_performance.csv",
          quote = FALSE, 
          row.names = FALSE
)

View(results)

#-------------------------------------------------------------------------------
# RANKING PERFORMANCE PER METRIC
#-------------------------------------------------------------------------------

# Rank the tools per performance metric ----

results_ordered <- results[, c("dataset", "strategy", "precision", "fnr", "pr_auc")]

# Function to rank and assign points for precision
rank_and_assign_points <- function(df, metric) {
  df %>%
    arrange(desc(!!sym(metric))) %>%
    mutate(rank = row_number()) %>%
    mutate(points = 11 - rank) %>%
    select(strategy, points)
}

ranked_results <- results_ordered %>% group_by(dataset)

precision_ranked <- rank_and_assign_points(ranked_results, "precision") %>% 
  group_by(strategy) %>% 
  summarise(precision = sum(points))

fnr_ranked <- rank_and_assign_points(ranked_results, "fnr") %>% 
  group_by(strategy) %>% 
  summarise(fnr = sum(points))
fnr_ranked$fnr <- 51 - fnr_ranked$fnr # inverse


pr_auc_ranked <- rank_and_assign_points(ranked_results, "pr_auc") %>% 
  group_by(strategy) %>% 
  summarise(pr_auc = sum(points))


ranked_results_output <- merge(precision_ranked, fnr_ranked, by = c("strategy"))
ranked_results_output <- merge(ranked_results_output, pr_auc_ranked, by = c("strategy"))

# Save and view final results
write.csv(ranked_results_output, 
          file = "/home/alicen/Projects/ReviewArticle/benchmark_results/performance_metrics/performance_ranked.csv",
          quote = FALSE, 
          row.names = FALSE
)

View(ranked_results_output)

#-------------------------------------------------------------------------------
# PLOT RANKS 
#-------------------------------------------------------------------------------

# Plot ranked performance bar graphs ----

# Optional to read in results (if not in global environment)
ranked_results <- read.csv(file = "/home/alicen/Projects/ReviewArticle/benchmark_results/performance_metrics/performance_ranked.csv")

# Add ranking columns
ranked_results <- ranked_results %>%
  mutate(
    precision_rank = rank(-precision, ties.method = "first"),
    fnr_rank = rank(-fnr, ties.method = "first"),
    pr_auc_rank = rank(-pr_auc, ties.method = "first")
  )

# Colour palette
strategy_colours <- c("#999FA9","#586667","#73966E","#415741","#65719E","#383F54","#E8E5D6","#A4A388","#8B3860","#522239")

# Convert the strategy column to a factor with the specified order
ranked_results$strategy <- factor(ranked_results$strategy, levels = methods)

# Define the common theme function
rank_theme <- function() {
  theme_classic() +  
    theme(
      axis.ticks = element_blank(), 
      axis.text = element_blank(), 
      axis.title = element_blank(),
      legend.position = "none",
      panel.border = element_rect(colour = "black", fill = NA),
      plot.title = element_text(vjust = 1, hjust = 0.5, face = "bold", size = 16),
      plot.margin = margin(20, 10, 10, 30)
    )
}

# Create the plots
ranked_precision <- ggplot(ranked_results, aes(x = strategy, y = precision, fill = strategy)) + 
  geom_bar(stat = "identity") + 
  geom_text(aes(label = precision_rank), vjust = -0.5) +
  scale_fill_manual(values = strategy_colours) +  
  labs(title = "Ranked Precision", y = "Precision") + 
  rank_theme() + ylim(0, 45)

ranked_fnr <- ggplot(ranked_results, aes(x = strategy, y = fnr, fill = strategy)) + 
  geom_bar(stat = "identity") + 
  geom_text(aes(label = fnr_rank), vjust = -0.5) +
  scale_fill_manual(values = strategy_colours) +  
  labs(title = "Ranked FNR", y = "FNR") + 
  rank_theme() + ylim(0, 45)

ranked_prauc <- ggplot(ranked_results, aes(x = strategy, y = pr_auc, fill = strategy)) + 
  geom_bar(stat = "identity") + 
  geom_text(aes(label = pr_auc_rank), vjust = -0.5) +
  scale_fill_manual(values = strategy_colours) +  
  labs(title = "Ranked PR-AUC", y = "PR-AUC") + 
  rank_theme() + ylim(0, 49)


#-------------------------------------------------------------------------------
# PLOT PERFORMANCE PER DATASET
#-------------------------------------------------------------------------------

# Bar plot for visualising metrics per dataset (10 bars (tools) x 5 plots (datasets) x 3 rows (metrics)) ----

# Optional to read in results (if not in global environment)
results <- read.csv(file = "/home/alicen/Projects/ReviewArticle/benchmark_results/performance_metrics/groundtruth_performance.csv")


# Convert the strategy column to a factor with the specified order
results$strategy <- factor(results$strategy, levels = methods)

# Pivot the data frame to long format
df_long <- results %>%
  pivot_longer(cols = c(precision, fnr, pr_auc), names_to = "measure", values_to = "value") %>%
  mutate(
    lower = case_when(
      measure == "precision" ~ precision_lower,
      measure == "fnr" ~ fnr_lower,
      measure == "pr_auc" ~ pr_auc_lower
    ),
    upper = case_when(
      measure == "precision" ~ precision_upper,
      measure == "fnr" ~ fnr_upper,
      measure == "pr_auc" ~ pr_auc_upper
    )
  ) %>%
  select(-precision_lower, -fnr_lower, -pr_auc_lower, -precision_upper, -fnr_upper, -pr_auc_upper) %>%
  distinct()

# Precision plot 
precision_plot <- ggplot(df_long %>% dplyr::filter(measure == "precision"), aes(x = strategy, y = value, fill = strategy)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  facet_wrap(~ dataset, nrow = 1, strip.position = "bottom") +
  scale_fill_manual(values = strategy_colors) +
  labs(title = "Precision", y = "Precision") + 
  theme_minimal() +  
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
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey"),
    panel.grid.minor.y = element_line(color = "grey"),
    panel.spacing = unit(1.5, "lines"),
    panel.border = element_blank(),
    plot.title = element_text(vjust = 0.4, hjust = 0.5, face = "bold", size = 16),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(color = "black", size = 12, face = "bold"),
    plot.margin = margin(10, 10, 50, 10)
  )

# FNR plot
fnr_plot <- ggplot(df_long %>% dplyr::filter(measure == "fnr"), aes(x = strategy, y = value, fill = strategy)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  facet_wrap(~ dataset, nrow = 1, strip.position = "bottom") +
  scale_fill_manual(values = strategy_colors) +
  labs(title = "False Negative Rate (FNR)", y = "FNR") + 
  theme_minimal() +  
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
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey"),
    panel.grid.minor.y = element_line(color = "grey"),
    panel.spacing = unit(1.5, "lines"),
    panel.border = element_blank(),
    plot.title = element_text(vjust = 0.4, hjust = 0.5, face = "bold", size = 16),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(color = "black", size = 12, face = "bold"),
    plot.margin = margin(10, 10, 50, 10) 
  )

# PR-AUC plot
pr_auc_plot <- ggplot(df_long %>% dplyr::filter(measure == "pr_auc"), aes(x = strategy, y = value, fill = strategy)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  facet_wrap(~ dataset, nrow = 1, strip.position = "bottom") +
  scale_fill_manual(values = strategy_colors) +
  scale_y_continuous(limits = c(0, 1), labels = scales::number_format(accuracy = 0.01)) +
  labs(title = "Precision-Recall Area Under the Curve (PR-AUC)", y = "Precision") + 
  theme_minimal() +  
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = "bottom",
    legend.justification = "left",
    axis.text = element_text(size = 16),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey"),
    panel.grid.minor.y = element_line(color = "grey"),
    panel.spacing = unit(1.5, "lines"),
    panel.border = element_blank(),
    plot.title = element_text(vjust = 0.4, hjust = 0.5, face = "bold", size = 16),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(color = "black", size = 12, face = "bold"),
    plot.margin = margin(10, 10, 20, 10) 
  )


#-------------------------------------------------------------------------------
# COMBINE PLOTS 
#-------------------------------------------------------------------------------

# Arrange the plots in a single layout with specified widths
final_plot <- (precision_plot + ranked_precision + plot_layout(widths = c(5, 1))) /
  (fnr_plot + ranked_fnr + plot_layout(widths = c(5, 1))) /
  (pr_auc_plot + ranked_prauc + plot_layout(widths = c(5, 1)))


# Save the plot
ggsave(filename = file.path("/home/alicen/Projects/ReviewArticle/benchmark_results/performance_metrics/barplots_ranked.png"), 
       plot = final_plot , width = 18, height = 13, dpi = 300)
ggsave(filename = file.path("/home/alicen/Projects/ReviewArticle/benchmark_results/performance_metrics/barplots_ranked.svg"), 
       plot = final_plot, width = 16, height = 13, dpi = 300)


