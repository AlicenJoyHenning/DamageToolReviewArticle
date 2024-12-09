# SCRIPT CONTEXT 
#
# Benchmarking damaged cell detection strategies on ground truth dataset. This script begins with the strategies
# already run on the ground truth cases. From here, the script uses the labels stored in the objects to 
# calculate confusion metrics (TP, TN, FP, FN) for each strategy, followed by performance metrics: 
#
# 1. Precision 
# 2. False Negative Rate (FNR)
# 3. Precision-recall area under the curve (PR-AUC)


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


# Load datasets (5 cases) - note after SampleQC different to after benchmarking for some 
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
methods <- c("ddqc", "DropletQC", "ensembleKQC", "miQC", "SampleQC",  "scater",  "valiDrops", 
             "manual_fixed_mito", "manual_adaptive_mito", "manual_mito_ribo",  "manual_mito_ribo_library", "manual_library", "manual_malat1", "manual_malat1_mito_ribo")
# Input structure
input_list <- list(
  list(name = "apoptotic", df = apoptotic, TP = "HEK293_apoptotic", TN = "HEK293_control"),
  list(name = "pro_apoptotic", df = pro_apoptotic, TP = "HEK293_proapoptotic", TN = "HEK293_control"),
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
combined_pr_auc_curves <- list()

# Combine items of the same dataset together
for (dataset in datasets) {
  dataset_items <- pr_auc_curves[grep(paste0("^", dataset, "_"), names(pr_auc_curves))]
  combined_pr_auc_curves[[dataset]] <- do.call(rbind, dataset_items)
}

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
          file = "./C_Test_Strategies/data/groundtruth_performance_metrics/groundtruth_performance.csv",
          quote = FALSE, 
          row.names = FALSE
)


#-------------------------------------------------------------------------------
# RANKING PERFORMANCE PER METRIC
#-------------------------------------------------------------------------------

# Rank the tools per performance metric ----

results_ordered <- results[, c("dataset", "strategy", "precision", "tpr", "pr_auc")]

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

tpr_ranked <- rank_and_assign_points(ranked_results, "tpr") %>% 
  group_by(strategy) %>% 
  summarise(tpr = sum(points))


pr_auc_ranked <- rank_and_assign_points(ranked_results, "pr_auc") %>% 
  group_by(strategy) %>% 
  summarise(pr_auc = sum(points))


ranked_results_output <- merge(precision_ranked, tpr_ranked, by = c("strategy"))
ranked_results_output <- merge(ranked_results_output, pr_auc_ranked, by = c("strategy"))

# Save and view final results
write.csv(ranked_results_output, 
          file = "./C_Test_Strategies/data/groundtruth_performance_metrics/performance_ranked.csv",
          quote = FALSE, 
          row.names = FALSE
)

View(ranked_results_output)


#-------------------------------------------------------------------------------
# PLOT PERFORMANCE PER DATASET
#-------------------------------------------------------------------------------

# Bar plot for visualising metrics per dataset (10 bars (tools) x 5 plots (datasets) x 3 rows (metrics)) ----

# Optional to read in results (if not in global environment)
results <- read.csv(file = "./C_Test_Strategies/data/groundtruth_performance_metrics/groundtruth_performance.csv")


# Convert the strategy column to a factor with the specified order

# Order according to set list (colour)
#results$strategy <- factor(results$strategy, levels = methods)
#results$dataset <- factor(results$dataset, levels = datasets)

# Order according to numerical values 
ordered_methods <- ranked_results$strategy[order(ranked_results$precision_rank)]
results$strategy <- factor(results$strategy, levels = ordered_methods)
results$dataset <- ifelse(results$dataset == "pro_apoptotic", "apoptotic_pro", results$dataset)


# Pivot the data frame to long format
df_long <- results %>%
  pivot_longer(cols = c(precision, tpr, pr_auc), names_to = "measure", values_to = "value") %>%
  mutate(
    lower = case_when(
      measure == "precision" ~ precision_lower,
      measure == "tpr" ~ tpr_lower,
      measure == "pr_auc" ~ pr_auc_lower
    ),
    upper = case_when(
      measure == "precision" ~ precision_upper,
      measure == "tpr" ~ tpr_upper,
      measure == "pr_auc" ~ pr_auc_upper
    )
  ) %>%
  select(-precision_lower, -tpr_lower, -pr_auc_lower, -precision_upper, -tpr_upper, -pr_auc_upper) %>%
  distinct()

# Add inverse for stacked bar effect 
df_long$inverse <- 1- df_long$value

# Reshape & ensure performance is plotted first
df_long_melted <- melt(df_long, id.vars = c("dataset", "strategy", "prevalence", "measure", "lower", "upper"), 
                       measure.vars = c("value", "inverse"))

# Ensure the variable column is a factor with the correct levels
df_long_melted$variable <- factor(df_long_melted$variable, levels = c("inverse", "value"))


# Create the stacked bar graph
performance_plot <- ggplot(df_long_melted, aes(x = value, y = strategy, fill = variable)) +
  geom_bar(stat = "identity", width = 0.55) +
  facet_wrap(~ dataset, nrow = 1, strip.position = "bottom") +
  scale_fill_manual(values = c("value" = "#001E5C", "inverse" = "#E6E6E6")) +
  theme_classic() + 
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "none",
        strip.text = element_blank(), # Comment out to verify order 
        strip.background = element_blank())


# Define theme 
performance_bar_theme <- function() {
  theme_minimal() +  
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 16),
      legend.position = "none",
      legend.justification = "center",
      axis.text = element_text(size = 18),
     # axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.y = element_line(color = "black"),
      #axis.line.x = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(1.5, "lines"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(vjust = 3.8, face = "bold", size = 16),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(color = "black", size = 12, face = "bold"),
      plot.margin = margin(10, 10, 50, 10)
    )
}

# Precision plot 
precision_plot <- ggplot(df_long %>% dplyr::filter(measure == "precision"), aes(x = strategy, y = value, fill = strategy)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  facet_wrap(~ dataset, nrow = 1, strip.position = "bottom") +
  #scale_fill_manual(values = strategy_colours) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + 
  labs(title = "", y = "") + 
  performance_bar_theme() + 
  theme(plot.title = element_text(hjust = -0.05),
        axis.title.x = element_blank(),
        axis.text.x = element_text(vjust = -1.5),
        strip.text = element_blank(), # Comment out to verify order 
        strip.background = element_blank())


precision_plot

ggsave(filename = "./C_Test_Strategies/img/precision.png",
       plot = precision_plot,
       width = 22, height = 8.2, units = "in")

# FNR plot

tpr_plot <- ggplot(df_long %>% dplyr::filter(measure == "tpr"), aes(x = strategy, y = value, fill = strategy)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  facet_wrap(~ dataset, nrow = 1, strip.position = "bottom") +
  scale_fill_manual(values = strategy_colours) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + 
  labs(title = "", y = "") + 
  performance_bar_theme() + 
  theme(plot.title = element_text(hjust = -0.05),
        axis.title.x = element_blank(),
        axis.text.x = element_text(vjust = -1.5),
        strip.text = element_blank(), 
        strip.background = element_blank())

tpr_plot

ggsave(filename = "./C_Test_Strategies/img/true_positive_rate.png",
       plot = tpr_plot,
       width = 22, height = 8.2, units = "in")

# PR-AUC plot ----- 
# Before show PR-AUC values, find what the value would be for random guessing 
datasets <- c("apoptotic", "pro_apoptotic", "GM18507_dead", "GM18507_dying", "PDX")
combined_pr_auc_curves <- list()

# Combine items of the same dataset together
for (dataset in datasets) {
  dataset_items <- pr_auc_curves[grep(paste0("^", dataset, "_"), names(pr_auc_curves))]
  combined_pr_auc_curves[[dataset]] <- do.call(rbind, dataset_items)
}

# Find the positive prevalence for each dataset 
positive_prevalence <- list()
positive_prevalence$apoptotic <- table(apoptotic$orig.ident)[1] / (table(apoptotic$orig.ident)[1] + table(apoptotic$orig.ident)[2])
positive_prevalence$apoptotic_pro <- table(pro_apoptotic$orig.ident)[2] / (table(pro_apoptotic$orig.ident)[1] + table(pro_apoptotic$orig.ident)[2])
positive_prevalence$GM18507_dead <- table(GM18507_dead$orig.ident)[2] / (table(GM18507_dead$orig.ident)[2] + table(GM18507_dead$orig.ident)[1])
positive_prevalence$GM18507_dying <- table(GM18507_dying$orig.ident)[2] / (table(GM18507_dying$orig.ident)[2] + table(GM18507_dying$orig.ident)[1])
positive_prevalence$PDX <- table(PDX$orig.ident)[2] / (table(PDX$orig.ident)[2] + table(PDX$orig.ident)[1])

# Convert positive_prevalence to a dataframe for easier plotting
positive_prevalence_df <- data.frame(
  dataset = names(positive_prevalence),
  prevalence = unlist(positive_prevalence)
)


# Bar plot like precision and FNR but for PR-auC with threshold line for random guessing 
pr_auc_plot <- ggplot(df_long %>% dplyr::filter(measure == "pr_auc"), aes(x = strategy, y = value, fill = strategy)) +
  geom_bar(stat = "identity") +
  geom_hline(data = positive_prevalence_df, aes(yintercept = prevalence), linetype = "dashed", color = "black") +
  coord_flip() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  facet_wrap(~ dataset, nrow = 1, strip.position = "bottom") +
  scale_fill_manual(values = strategy_colours) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + 
  labs(title = "", y = "") + 
  performance_bar_theme() + 
  theme(plot.title = element_text(hjust = -0.05),
        axis.title.x = element_blank(),
        axis.text.x = element_text(vjust = -1.5),
        strip.text = element_blank(), 
        strip.background = element_blank())

pr_auc_plot

ggsave(filename = "./C_Test_Strategies/img/precision_recall_AUC.png",
       plot = pr_auc_plot,
       width = 22, height = 8.2, units = "in")


#-------------------------------------------------------------------------------
# COMBINE PLOTS 
#-------------------------------------------------------------------------------

# Arrange the plots in a single layout with specified widths and row labels
final_plot <- (precision_plot + ranked_precision + plot_layout(widths = c(5, 1))) / 
  (fnr_plot + ranked_fnr + plot_layout(widths = c(5, 1))) / 
  (pr_auc_plot + ranked_prauc + plot_layout(widths = c(5, 1)))

# Save the plot
ggsave(filename = file.path("./C_Test_Strategies/img/barplots_ranked.png"), 
       plot = final_plot , width = 25, height = 19, dpi = 300, limitsize = FALSE)



### End

