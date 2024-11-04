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

# Read in damaged-perturbed simulated datasets and labels using loop ----

# Parent directories for the rds objects (housing count matrices)
parent_directory_seurat <- "/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/R_tool_input/"
conditions <- c("control_sim", "stimulated_sim")
percentages <- c("2.5", "5", "10", "15", "20")
reps <- 1:3

# Read in and store all objects in a list 
simulated <- list()
for (condition in conditions) {
  for (percentage in percentages) {
    for (rep in reps) {
      # Construct the file name
      file_name <- paste0(condition, "_", rep, "_", percentage, ".rds")
      file_path <- file.path(parent_directory_seurat, file_name)
      
      # Check if the file exists before reading
      if (file.exists(file_path)) {
        
        # Read the seurat object 
        data <- readRDS(file_path)
        
        # Store the data in the list with a meaningful name
        simulated[[paste0(condition, "_", rep, "_", percentage, "_seurat")]] <- data
      }
    }
  }
}

# Damaged strategy outputs (dfs)
parent_directory_df <- "/home/alicen/Projects/ReviewArticle/damage_perturbation/benchmark_results/"

# Store all resulting output dfs in list 
results <- list()
for (condition in conditions) {
  for (percentage in percentages) {
    for (rep in reps) {
      # Construct the file name
      file_name <- paste0(condition, "_", rep, "_", percentage, ".csv")
      file_path <- file.path(parent_directory_df, file_name)
      
      # Check if the file exists before reading
      if (file.exists(file_path)) {
        
        # Read the seurat object 
        data <- read.csv(file_path)
        
        # Store the data in the list with a meaningful name
        results[[paste0(condition, "_", rep, "_", percentage)]] <- data
      }
    }
  }
}


# Transfer results to seurat objects -----

for (condition in conditions) {
  for (percentage in percentages) {
    for (rep in reps) {
      # Construct the key for the simulated list
      seurat_key <- paste0(condition, "_", rep, "_", percentage, "_seurat")
      result_key <- paste0(condition, "_", rep, "_", percentage)
      
      # Transfer the data frame to the seurat object's meta.data
      simulated[[seurat_key]]@meta.data <- results[[result_key]]
      rownames(simulated[[seurat_key]]@meta.data) <- colnames(simulated[[seurat_key]]@assays$RNA$counts)
      
    }
  }
}


#-------------------------------------------------------------------------------
# PERFORMANCE METRICS
#-------------------------------------------------------------------------------

# Calculate the performance metrics ----

# List of methods
methods <- c("ddqc", "DropletQC", "ensembleKQC", "miQC", "SampleQC",  "scater",  "valiDrops", 
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


# Initialize a list to store PR-AUC curve data
pr_auc_curves <- list()

# Calculate precision, FNR, and PR-AUC for each strategy for each dataset
for (item in input_list) {
  
  # Extract components of the list 
  df <- item$df
  TP <- item$TP
  TN <- item$TN
  dataset_name <- item$name
  
  # Loop through each method
  for (method in methods) {
    
    cat("\n")
    message(method, " ...")
    
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
    
    
    
    # Calculate PR-AUC curve
    if (method == "DropletQC") { next } else {
      
      scores <- c(rep(1, TP_count), rep(0, FN_count), rep(1, FP_count), rep(0, TN_count))
      labels <- c(rep(1, TP_count + FN_count), rep(0, FP_count + TN_count))
      pr <- pr.curve(scores.class0 = scores, weights.class0 = labels, curve = TRUE)
      
      # Store the curve data
      pr_auc_curves[[paste0(dataset_name, "_", method)]] <- data.frame(
        recall = pr$curve[, 1],
        precision = pr$curve[, 2],
        method = method,
        dataset = dataset_name
      )
    
    }
    
  }
}

# Save and view final results
write.csv(results, 
          file = "/home/alicen/Projects/ReviewArticle/benchmark_results/performance_metrics/groundtruth_performance.csv",
          quote = FALSE, 
          row.names = FALSE
)

View(results)


# Plot ranked performance bar graphs ----

# Define the colors
strategy_colours <- c(
  "ddqc" = "#CE9DBA",
  "DropletQC" = "#A799C9",
  "ensembleKQC" = "#E7E4F6", 
  "miQC" = "#808C98", 
  "SampleQC" = "#CED5DB",
  "scater" = "#88A1BD", 
  "valiDrops" = "#D3E2F6",
  "manual_all" = "#4F618F",
  "manual_mito_isolated" = "#A6BEAE", 
  "manual_mito" = "#DCECE2", 
  "manual_mito_ribo" = "#9DBD78",
  "manual_malat1" = "#D7E7BE"
)

# Define the order of strategies
strategy_order <- c(
  "ddqc", "DropletQC", "ensembleKQC", "miQC", "SampleQC", "scater", "valiDrops",
  "manual_all", "manual_mito_isolated", "manual_mito", "manual_mito_ribo", "manual_malat1"
)

# Convert the strategy column to a factor with the specified order
ranked_results$strategy <- factor(ranked_results$strategy, levels = strategy_order )

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
results$strategy <- factor(results$strategy, levels = strategy_order)

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


# Define theme 
performance_bar_theme <- function() {
  theme_minimal() +  
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 16),
      legend.position = "none",
      legend.justification = "center",
      axis.text = element_text(size = 18),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.y = element_line(color = "black"),
      axis.line.x = element_blank(),
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


# PR-AUC plot
pr_auc_plot <- ggplot(df_long %>% dplyr::filter(measure == "pr_auc"), aes(x = strategy, y = value, fill = strategy)) +
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

ggsave(filename = "/home/alicen/Projects/ReviewArticle/benchmark_results/performance_metrics/precision_recall_AUC.png",
       plot = pr_auc_plot,
       width = 17, height = 6.7, units = "in")

#  PR-AUC curve plot ----

create_pr_auc_plot <- function(dataset, data, strategy_colours) {
  p <- ggplot(data %>% filter(dataset == dataset), aes(x = recall, y = precision, color = method)) +
    geom_line(size = 1) +
    scale_color_manual(values = strategy_colours) +
    labs(title = paste("PR-AUC Curve for", dataset),
         x = "Recall",
         y = "Precision") +
    theme_classic() +
    theme(
      legend.position = "right",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  return(p)
}

# Rename list items 
datasets <- c("apoptotic", "proapoptotic", "GM18507_dead", "GM18507_dying", "PDX")

# Initialize a new list to store combined dataframes for each dataset
combined_pr_auc_curves <- list()

# Combine items of the same dataset together
for (dataset in datasets) {
  # Filter items that belong to the current dataset
  dataset_items <- pr_auc_curves[grep(paste0("^", dataset, "_"), names(pr_auc_curves))]
  
  # Combine the items into a single dataframe
  combined_pr_auc_curves[[dataset]] <- do.call(rbind, dataset_items)
}

# Find the positive prevalence for each dataset 
positive_prevalence <- list()

positive_prevalence$apoptotic <- table(apoptotic$orig.ident)[1] / (table(apoptotic$orig.ident)[1] + table(apoptotic$orig.ident)[2])
positive_prevalence$proapoptotic <- table(pro_apoptotic$orig.ident)[2] / (table(pro_apoptotic$orig.ident)[1] + table(pro_apoptotic$orig.ident)[2])
positive_prevalence$GM18507_dead <- table(GM18507_dead$orig.ident)[2] / (table(GM18507_dead$orig.ident)[2] + table(GM18507_dead$orig.ident)[1])
positive_prevalence$GM18507_dying <- table(GM18507_dying$orig.ident)[2] / (table(GM18507_dying$orig.ident)[2] + table(GM18507_dying$orig.ident)[1])
positive_prevalence$PDX <- table(PDX$orig.ident)[2] / (table(PDX$orig.ident)[2] + table(PDX$orig.ident)[1])


# Plot the curves 
plots <- list()
for (dataset in names(combined_pr_auc_curves)) {
  p <- ggplot(combined_pr_auc_curves[[dataset]], aes(x = recall, y = precision, color = method)) +
    geom_line(size = 2) +
    geom_hline(yintercept = positive_prevalence[[dataset]], linetype = "dashed", color = "black") +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +  
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    scale_color_manual(values = strategy_colours) +
    labs(title = "",
         x = "Recall",
         y = "Precision") +
    theme_classic() +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_text(size = 16, hjust = -0.5),
          axis.text.x = element_text(size = 16, vjust = -0.5),
          panel.border = element_rect(fill = NA, color = "black"),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA))
  
  plots[[dataset]] <- p
}

# Name the plots
names(plots) <- names(combined_pr_auc_curves)

PRAUC_plots <- plots$apoptotic | plots$GM18507_dead | plots$GM18507_dying | plots$PDX | plots$proapoptotic

ggsave(filename = "/home/alicen/Projects/ReviewArticle/benchmark_results/performance_metrics/PRAUC_curves.png",
       plot = PRAUC_plots,
       width = 15, height = 5.3, units = "in")


### End

