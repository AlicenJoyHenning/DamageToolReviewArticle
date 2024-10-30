# SCRIPT CONTEXT 
#
# Benchmarking damaged cell detection strategies. The labelled output objects are loaded at the start of the script. 
#   (3 - Run remaining detection strategies.R) From here, the script uses the labels stored in the objects to calculate: 
#
# 1. The proportion of damaged cells for each tool - violin plot 
# 2. The similarity between the tools (Cohen's Kappa)-  PCA plots
# 3. Consistency score : combines proportion damaged & similarity - bar plot 
# 4. Dimensionality reduction plots - UMAPs or tSNEs using mt-rb genes 



#-------------------------------------------------------------------------------
# PREPARATIONS
#-------------------------------------------------------------------------------

# Load necessary libraries
library(cowplot)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(purrr)
library(irr)

# Load datasets  -----
# Ground truth 
apoptotic <- read.csv("/home/alicen/Projects/ReviewArticle/benchmark_results/HEK293_apo.csv")
pro_apoptotic <- read.csv("/home/alicen/Projects/ReviewArticle/benchmark_results/HEK293_pro.csv")
GM18507_dead <- read.csv("/home/alicen/Projects/ReviewArticle/benchmark_results/GM18507_dead.csv")
GM18507_dying <- read.csv("/home/alicen/Projects/ReviewArticle/benchmark_results/GM18507_dying.csv")
PDX <- read.csv("/home/alicen/Projects/ReviewArticle/benchmark_results/PDX_dead.csv")

# Non-ground truth 
A549 <- read.csv("/home/alicen/Projects/ReviewArticle/benchmark_results/A549.csv")
dLiver <- read.csv("/home/alicen/Projects/ReviewArticle/benchmark_results/dLiver.csv")
dLung <- read.csv("/home/alicen/Projects/ReviewArticle/benchmark_results/dLung.csv")
dPBMC <- read.csv("/home/alicen/Projects/ReviewArticle/benchmark_results/dPBMC.csv")
ductal <- read.csv("/home/alicen/Projects/ReviewArticle/benchmark_results/ductal.csv")
glio <- read.csv("/home/alicen/Projects/ReviewArticle/benchmark_results/glio.csv")
HCT116 <- read.csv("/home/alicen/Projects/ReviewArticle/benchmark_results/HCT116.csv")
hLiver <- read.csv("/home/alicen/Projects/ReviewArticle/benchmark_results/hLiver.csv")
hLung <- read.csv("/home/alicen/Projects/ReviewArticle/benchmark_results/hLung.csv")
hodgkin <- read.csv("/home/alicen/Projects/ReviewArticle/benchmark_results/hodgkin.csv")
hPBMC <- read.csv("/home/alicen/Projects/ReviewArticle/benchmark_results/hPBMC.csv")
Jurkat <- read.csv("/home/alicen/Projects/ReviewArticle/benchmark_results/Jurkat.csv")
mLiver <- read.csv("/home/alicen/Projects/ReviewArticle/benchmark_results/mLiver.csv")
mLung <- read.csv("/home/alicen/Projects/ReviewArticle/benchmark_results/mLung.csv")
mPBMC <- read.csv("/home/alicen/Projects/ReviewArticle/benchmark_results/mPBMC.csv")


# Simulated 
# Damaged strategy outputs (dfs)
parent_directory_df <- "/home/alicen/Projects/ReviewArticle/damage_perturbation/benchmark_results/"

# Store all resulting output dfs in list 
simulated <- list()
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
        simulated[[paste0(condition, "_", rep, "_", percentage)]] <- data
      }
    }
  }
}

# Check same 
View(simulated$control_sim_1_2.5)

# Define palette

# Define the colors
strategy_colours <- c(
  "ddqc" = "#A799C9",
  "DropletQC" = "#E7E4F6",
  "ensembleKQC" = "#808C98", 
  "miQC" = "#CED5DB", 
  "scater" = "#88A1BD", 
  "valiDrops" = "#D3E2F6",
  "manual_all" = "#4F618F",
  "manual_mito_isolated" = "#A6BEAE", 
  "manual_mito" = "#DCECE2", 
  "manual_mito_ribo" = "#9DBD78",
  "manual_malat1" = "#D7E7BE"
)


# Define the desired order
desired_order <-  c("ddqc", "DropletQC", "ensembleKQC", "miQC", "scater", "valiDrops",
                    "manual_all", "manual_mito_isolated", "manual_mito", "manual_mito_ribo", "manual_malat1"
)


# Define theme 
barplot_theme <- function() {
  theme_minimal() +  
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 16),
      legend.position = "none",
      legend.justification = "center",
      axis.text = element_text(size = 18),
      # axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.spacing = unit(1.5, "lines"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      plot.title = element_text(vjust = 3.8, hjust = 0, face = "bold", size = 16),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(color = "black", size = 12, face = "bold"),
      plot.margin = margin(50, 20, 50, 20) # top right bottom left 
    )
}


#-------------------------------------------------------------------------------
# PROPORTION DAMAGED 
#-------------------------------------------------------------------------------

# Proportion damaged as bar plot ----

find_proportion_damaged <- function(data, project_name){
  
  data <- data[, c("ddqc", "DropletQC", "ensembleKQC", "miQC", "scater", "valiDrops",
                   "manual_all", "manual_mito_isolated", "manual_mito", "manual_mito_ribo", "manual_malat1")]
  
  # Calculate damaged proportion for each column
  damaged_proportion <- sapply(data, function(column) {
    damaged_count <- sum(column == "damaged")
    total_count <- length(column)
    (damaged_count / total_count) * 100
  })
  
  # Convert to data frame for better readability
  damaged_proportion_df <- data.frame(strategy = names(damaged_proportion), damaged_proportion = damaged_proportion)
  damaged_proportion_df <- setNames(damaged_proportion_df, c("strategy", project_name))
  
  return(damaged_proportion_df)
  
}

# Run through all datasets
# Ground truth 
apoptotic_df <- find_proportion_damaged(apoptotic, "apoptotic")
pro_apoptotic_df <- find_proportion_damaged(pro_apoptotic, "pro_apoptotic")
GM18507_dead_df <- find_proportion_damaged(GM18507_dead, "GM18507_dead")
GM18507_dying_df <- find_proportion_damaged(GM18507_dying, "GM18507_dying")
PDX_df <- find_proportion_damaged(PDX, "PDX")

# Non ground truth 
A549_df <- find_proportion_damaged(A549, "A549")
dLiver_df <- find_proportion_damaged(dLiver, "dLiver")
dLung_df <- find_proportion_damaged(dLung, "dLung")
dPBMC_df <- find_proportion_damaged(dPBMC, "dPBMC")
ductal_df <- find_proportion_damaged(ductal, "ductal")
glio_df <- find_proportion_damaged(glio, "glio")
HCT116_df <- find_proportion_damaged(HCT116, "HCT116")
hLiver_df <- find_proportion_damaged(hLiver, "hLiver")
hLung_df <- find_proportion_damaged(hLung, "hLung")
hodgkin_df <- find_proportion_damaged(hodgkin, "hodgkin")
hPBMC_df <- find_proportion_damaged(hPBMC, "hPBMC")
Jurkat_df <- find_proportion_damaged(Jurkat, "Jurkat")
mLiver_df <- find_proportion_damaged(mLiver, "mLiver")
mLung_df <- find_proportion_damaged(mLung, "mLung") 
mPBMC_df <- find_proportion_damaged(mPBMC, "mPBMC")

# Simulated
simulated_df <- list()
for (name in names(simulated)){
  
  simulated_df[[name]] <-  find_proportion_damaged(simulated[[name]], as.character(name))

}


# List of data frames
non_groundtruth_dfs <- list(A549_df, dLiver_df, dLung_df, dPBMC_df, ductal_df, 
                            glio_df, HCT116_df, hLiver_df, hLung_df, hodgkin_df, 
                            hPBMC_df, Jurkat_df, mLiver_df, mLung_df, mPBMC_df)

groundtruth_dfs <- list(apoptotic_df, pro_apoptotic_df, 
                        GM18507_dead_df, GM18507_dying_df, 
                        PDX_df)

simulated_dfs <- list(simulated_df$control_sim_1_2.5, simulated_df$control_sim_1_5, simulated_df$control_sim_1_10, simulated_df$control_sim_1_15, simulated_df$control_sim_1_20,
                      simulated_df$control_sim_2_2.5, simulated_df$control_sim_2_5, simulated_df$control_sim_2_10, simulated_df$control_sim_2_15, simulated_df$control_sim_2_20,
                      simulated_df$control_sim_3_2.5, simulated_df$control_sim_3_5, simulated_df$control_sim_3_10, simulated_df$control_sim_3_15, simulated_df$control_sim_3_20,
                      simulated_df$stimulated_sim_1_2.5, simulated_df$stimulated_sim_1_5, simulated_df$stimulated_sim_1_10, simulated_df$stimulated_sim_1_15, simulated_df$stimulated_sim_1_20,
                      simulated_df$stimulated_sim_2_2.5, simulated_df$stimulated_sim_2_5, simulated_df$stimulated_sim_2_10, simulated_df$stimulated_sim_2_15, simulated_df$stimulated_sim_2_20,
                      simulated_df$stimulated_sim_3_2.5, simulated_df$stimulated_sim_3_5, simulated_df$stimulated_sim_3_10, simulated_df$stimulated_sim_3_15, simulated_df$stimulated_sim_3_20
)


# Merge all data frames by 'strategy' column
non_groundtruth_df <- reduce(non_groundtruth_dfs, full_join, by = "strategy")
groundtruth_df <- reduce(groundtruth_dfs, full_join, by = "strategy")
simulated_full_df <- reduce(simulated_dfs, full_join, by = "strategy")

# Find median values 
groundtruth_df <- groundtruth_df %>%
  rowwise() %>%
  mutate(
    lowest_value = min(c_across(apoptotic:PDX)),
    highest_value = max(c_across(apoptotic:PDX)),
    median_value = median(c_across(apoptotic:PDX))
  ) %>%
  ungroup() %>% 
  mutate(strategy = factor(strategy, levels = desired_order))

non_groundtruth_df <- non_groundtruth_df %>%
  rowwise() %>%
  mutate(
    lowest_value = min(c_across(A549:mLung)),
    highest_value = max(c_across(A549:mLung)),
    median_value = median(c_across(A549:mLung))
  ) %>%
  ungroup() %>%
  mutate(strategy = factor(strategy, levels = desired_order))

simulated_full_df <- simulated_full_df %>%
  rowwise() %>%
  mutate(
    lowest_value = min(c_across(control_sim_1_2.5:stimulated_sim_3_20)),
    highest_value = max(c_across(control_sim_1_2.5:stimulated_sim_3_20)),
    median_value = median(c_across(control_sim_1_2.5:stimulated_sim_3_20))
  ) %>%
  ungroup() %>%
  mutate(strategy = factor(strategy, levels = desired_order))



# Damaged proportion plots
proportion_damaged_groundtruth_plot <- ggplot(groundtruth_df, aes(x = strategy, y = median_value, fill = strategy)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = lowest_value, ymax = highest_value), width = 0.2) +
  scale_fill_manual(values = strategy_colours) +
  coord_flip() + 
  scale_y_reverse(limits = c(60, 0)) +  # Reverse the y-axis
  labs(title = "a      Damaged proportion", y = "") + 
  barplot_theme() + 
  theme(plot.title = element_blank(),
        plot.margin = margin(50, 10, 50, 30), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 34, vjust = -0.5),
        legend.text = element_text(size = 22),
        legend.spacing.x = unit(1.0, 'cm'), 
        legend.position = "none",
        legend.justification = c(-0.05, -0.2), 
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)) + 
  guides(fill = guide_legend(nrow = 1, keywidth=0.4, keyheight=0.1, default.unit="inch"))



proportion_damaged_non_groundtruth_plot <- ggplot(non_groundtruth_df, aes(x = strategy, y = median_value, fill = strategy)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = lowest_value, ymax = highest_value), width = 0.2) +
  scale_fill_manual(values = strategy_colours) +
  coord_flip() + 
  scale_y_reverse(limits = c(60, 0)) +  # Reverse the y-axis
  labs(title = "a      Damaged proportion", y = "") + 
  barplot_theme() + 
  theme(plot.title = element_blank(),
        plot.margin = margin(50, 10, 50, 30), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 34, vjust = -0.5),
        legend.text = element_text(size = 22),
        legend.spacing.x = unit(1.0, 'cm'), 
        legend.position = "none",
        legend.justification = c(-0.05, -0.2), 
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)) + 
  guides(fill = guide_legend(nrow = 1, keywidth=0.4, keyheight=0.1, default.unit="inch"))


 
proportion_damaged_simulated_plot <- ggplot(simulated_full_df, aes(x = strategy, y = median_value, fill = strategy)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = lowest_value, ymax = highest_value), width = 0.2) +
  scale_fill_manual(values = strategy_colours) +
  coord_flip() + 
  scale_y_reverse(limits = c(60, 0)) +  # Reverse the y-axis
  labs(title = "a      Damaged proportion", y = "") + 
  barplot_theme() + 
  theme(plot.title = element_blank(),
        plot.margin = margin(50, 10, 50, 30), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 34, vjust = -0.5),
        legend.text = element_text(size = 22),
        legend.spacing.x = unit(1.0, 'cm'), 
        legend.position = "none",
        legend.justification = c(-0.05, -0.2), 
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)) + 
  guides(fill = guide_legend(nrow = 1, keywidth=0.4, keyheight=0.1, default.unit="inch"))



proportion_damaged_simulated_plot | proportion_damaged_non_groundtruth_plot | proportion_damaged_groundtruth_plot

#-------------------------------------------------------------------------------
# PROPORTION UNIQUE
#-------------------------------------------------------------------------------

# Proportion of uniquely labelled damaged cells per tool across datasets -----

find_proportion_unique <- function(data, project_name, groundtruth = NULL) {
  
  # Define tools and columns of interest for the data frame 
  tools <- c("ddqc", "DropletQC", "ensembleKQC", "miQC", "scater", "valiDrops", 
             "manual_all", "manual_mito_isolated", "manual_mito", "manual_mito_ribo", "manual_malat1")
  
  columns <- c("orig.ident", tools)
  data <- data[, columns]
  
  # Determine uniqueness by label from none or at most one other tool
  is_unique <- function(tool_label, other_labels) {
    if (tool_label == "damaged" && sum(other_labels == "damaged") <= 1) {
      return("unique")
    } else {
      return("common")
    }
  }
  
  for (tool in tools) {
    other_tools <- tools[tools != tool]
    data[[tool]] <- apply(data, 1, function(row) {
      is_unique(row[[tool]], row[other_tools])
    })
  }
  
  # Calculate unique proportion for each column
  unique_proportion <- sapply(data[, tools], function(column) {
    unique_count <- sum(column == "unique")
    total_count <- length(column)
    round((unique_count / total_count) * 100, 0)
  })
  
  # Convert to data frame for better readability
  unique_proportion_df <- data.frame(strategy = names(unique_proportion), unique_proportion = unique_proportion)
  unique_proportion_df <- setNames(unique_proportion_df, c("strategy", project_name))
  
  
  # If the input has ground truth labels, calculate the proportion of the unique labels that are correct 
  if (!is.null(groundtruth)) {
    
    # Function to determine if unique AND correct
    is_unique_and_correct <- function(tool_label, orig_ident, groundtruth) {
      if (tool_label == "unique" && orig_ident == groundtruth) {
        return("correct")
      } else if (tool_label == "unique") {
        return("incorrect")
      } else {
        return(tool_label)
      }
    }
    
    # Create object with isolated unique labels to store correctness
    data_correct <- data[apply(data[, tools], 1, function(row) any(row == "unique")), ]
    
    for (tool in tools) { 
      data_correct[[tool]] <- apply(data_correct, 1, function(row) {
        is_unique_and_correct(row[[tool]], row[["orig.ident"]], groundtruth)
      })
    }
    
    # Calculate the proportion of correct labels for every unique label
    correct_unique_proportion <- sapply(tools, function(tool) {
      correct_count <- sum(data_correct[[tool]] == "correct")
      unique_count <- sum(data_correct[[tool]] %in% c("correct", "incorrect"))
      if (unique_count == 0) {
        return(0)  
      } else {
        return(round((correct_count / unique_count) * 100, 0))
      }
    })
    
    # Convert to data frame for better readability
    correct_unique_proportion_df <- data.frame(strategy = tools, correct_unique_proportion = correct_unique_proportion)
    correct_unique_proportion_df <- setNames(correct_unique_proportion_df, c("strategy", project_name))
    
    return(list(unique = unique_proportion_df, 
                correct = correct_unique_proportion_df))
    
  } else {
    
    return(unique_proportion_df)
    
  }
}

# Run through all datasets to find unique damaged labels (udf)
A549_udf <- find_proportion_unique(A549, "A549")
dLiver_udf <- find_proportion_unique(dLiver, "dLiver")
dLung_udf <- find_proportion_unique(dLung, "dLung")
dPBMC_udf <- find_proportion_unique(dPBMC, "dPBMC")
ductal_udf <- find_proportion_unique(ductal, "ductal")
glio_udf <- find_proportion_unique(glio, "glio")
HCT116_udf <- find_proportion_unique(HCT116, "HCT116")
hLiver_udf <- find_proportion_unique(hLiver, "hLiver")
hLung_udf <- find_proportion_unique(hLung, "hLung")
hodgkin_udf <- find_proportion_unique(hodgkin, "hodgkin")
hPBMC_udf <- find_proportion_unique(hPBMC, "hPBMC")
Jurkat_udf <- find_proportion_unique(Jurkat, "Jurkat")
mLiver_udf <- find_proportion_unique(mLiver, "mLiver")
mLung_udf <- find_proportion_unique(mLung, "mLung") #
mPBMC_udf <- find_proportion_unique(mPBMC, "mPBMC")


simulated_udf <- list()
for (name in names(simulated)) {
  
  simulated_udf[[name]] <- find_proportion_unique(simulated[[name]], as.character(name))
  
}


# Groundtruth specify damaged cell label 
apoptotic_udf <- find_proportion_unique(apoptotic, "apoptotic", "HEK293_apoptotic")
pro_apoptotic_udf <- find_proportion_unique(pro_apoptotic, "pro_apoptotic", "HEK293_proapoptotic")
GM18507_dead_udf <- find_proportion_unique(GM18507_dead, "GM18507_dead", "GM18507_dead")
GM18507_dying_udf <- find_proportion_unique(GM18507_dying, "GM18507_dying", "GM18507_dying")
PDX_udf <- find_proportion_unique(PDX, "PDX", "PDX_dead")


# List of proportion Unique Data Frames (udf) and proportion Correct Unique Data Frames (cudf)
non_groundtruth_udfs <- list(A549_udf, dLiver_udf, dLung_udf, dPBMC_udf, ductal_udf, 
                            glio_udf, HCT116_udf, hLiver_udf, hLung_udf, hodgkin_udf, 
                            hPBMC_udf, Jurkat_udf, mLiver_udf, mLung_udf, mPBMC_udf)

groundtruth_udfs <- list(apoptotic_udf$unique, pro_apoptotic_udf$unique, 
                        GM18507_dead_udf$unique, GM18507_dying_udf$unique, 
                        PDX_udf$unique)

# Simulated
simulated_udfs <- list(simulated_udf$control_sim_1_2.5, simulated_udf$control_sim_1_5, simulated_udf$control_sim_1_10, simulated_udf$control_sim_1_15, simulated_udf$control_sim_1_20,
                       simulated_udf$control_sim_2_2.5, simulated_udf$control_sim_2_5, simulated_udf$control_sim_2_10, simulated_udf$control_sim_2_15, simulated_udf$control_sim_2_20,
                       simulated_udf$control_sim_3_2.5, simulated_udf$control_sim_3_5, simulated_udf$control_sim_3_10, simulated_udf$control_sim_3_15, simulated_udf$control_sim_3_20,
                       simulated_udf$stimulated_sim_1_2.5, simulated_udf$stimulated_sim_1_5, simulated_udf$stimulated_sim_1_10, simulated_udf$stimulated_sim_1_15, simulated_udf$stimulated_sim_1_20,
                       simulated_udf$stimulated_sim_2_2.5, simulated_udf$stimulated_sim_2_5, simulated_udf$stimulated_sim_2_10, simulated_udf$stimulated_sim_2_15, simulated_udf$stimulated_sim_2_20,
                       simulated_udf$stimulated_sim_3_2.5, simulated_udf$stimulated_sim_3_5, simulated_udf$stimulated_sim_3_10, simulated_udf$stimulated_sim_3_15, simulated_udf$stimulated_sim_3_20
)

groundtruth_cudfs <- list(apoptotic_udf$correct, pro_apoptotic_udf$correct, 
                         GM18507_dead_udf$correct, GM18507_dying_udf$correct, 
                         PDX_udf$correct)


# Merge all data frames by 'strategy' column
non_groundtruth_udf <- reduce(non_groundtruth_udfs, full_join, by = "strategy")
groundtruth_udf <- reduce(groundtruth_udfs, full_join, by = "strategy")
simulated_udf <- reduce(simulated_udfs, full_join, by = "strategy")
groundtruth_cudf <- reduce(groundtruth_cudfs, full_join, by = "strategy")


# Find median values
non_groundtruth_udf <- non_groundtruth_udf %>%
  rowwise() %>%
  mutate(
    lowest_value = min(c_across(A549:mLung)),
    highest_value = max(c_across(A549:mLung)),
    median_value = median(c_across(A549:mLung))
  ) %>%
  ungroup() %>%
  mutate(strategy = factor(strategy, levels = desired_order))

groundtruth_udf <- groundtruth_udf %>%
  rowwise() %>%
  mutate(
    lowest_value = min(c_across(apoptotic:PDX)),
    highest_value = max(c_across(apoptotic:PDX)),
    median_value = median(c_across(apoptotic:PDX))
  ) %>%
  ungroup() %>%
  mutate(strategy = factor(strategy, levels = desired_order))

simulated_udf <- simulated_udf %>%
  rowwise() %>%
  mutate(
    lowest_value = min(c_across(control_sim_1_2.5:stimulated_sim_3_20)),
    highest_value = max(c_across(control_sim_1_2.5:stimulated_sim_3_20)),
    median_value = median(c_across(control_sim_1_2.5:stimulated_sim_3_20))
  ) %>%
  ungroup() %>%
  mutate(strategy = factor(strategy, levels = desired_order))

groundtruth_cudf <- groundtruth_cudf %>%
  rowwise() %>%
  mutate(
    correct_percent = median(c_across(apoptotic:PDX))
  ) %>%
  ungroup() %>%
  mutate(strategy = factor(strategy, levels = desired_order))

 # Add the percent correctness to ground truth data frame
groundtruth_udf$correct_percent <- groundtruth_cudf$correct_percent


# Proportion unique plots 
proportion_unique_groundtruth_plot <- ggplot(groundtruth_udf, aes(x = strategy, y = median_value, fill = strategy)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = lowest_value, ymax = highest_value), width = 0.2) +
  scale_fill_manual(values = strategy_colours) +
  coord_flip() + 
  barplot_theme() + 
  theme(plot.title = element_blank(),
        plot.margin = margin(50, 20, 50, 30), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 34, vjust = -0.5),
        legend.text = element_text(size = 22),
        legend.spacing.x = unit(1.0, 'cm'), 
        legend.position = "none",
        legend.justification = c(-0.05, -0.2), 
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)) + 
  guides(fill = guide_legend(nrow = 1, keywidth=0.4, keyheight=0.1, default.unit="inch"))

ground_prop_plots <- proportion_damaged_groundtruth_plot | proportion_unique_groundtruth_plot 
ggsave(filename = file.path("/home/alicen/Projects/ReviewArticle/benchmark_results/comparison_metrics/proportions_groundtruth_plot.png"), 
       plot = ground_prop_plots, width = 14, height = 12, dpi = 300)


proportion_unique_non_groundtruth_plot <- ggplot(non_groundtruth_udf, aes(x = strategy, y = median_value, fill = strategy)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = lowest_value, ymax = highest_value), width = 0.2) +
  scale_fill_manual(values = strategy_colours) +
  coord_flip() + 
  barplot_theme() + 
  theme(plot.title = element_blank(),
        plot.margin = margin(50, 20, 50, 30), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 34, vjust = -0.5),
        legend.text = element_text(size = 22),
        legend.spacing.x = unit(1.0, 'cm'), 
        legend.position = "none",
        legend.justification = c(-0.05, -0.2), 
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)) + 
  guides(fill = guide_legend(nrow = 1, keywidth=0.4, keyheight=0.1, default.unit="inch"))

non_ground_prop_plots <- proportion_damaged_non_groundtruth_plot | proportion_unique_non_groundtruth_plot 
ggsave(filename = file.path("/home/alicen/Projects/ReviewArticle/benchmark_results/comparison_metrics/proportions_non_groundtruth_plot.png"), 
       plot = non_ground_prop_plots, width = 14, height = 12, dpi = 300)


proportion_unique_simulated_plot <- ggplot(simulated_udf, aes(x = strategy, y = median_value, fill = strategy)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = lowest_value, ymax = highest_value), width = 0.2) +
  scale_fill_manual(values = strategy_colours) +
  coord_flip() + 
  barplot_theme() + 
  theme(plot.title = element_blank(),
        plot.margin = margin(50, 20, 50, 10), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 34, vjust = -0.5),
        legend.text = element_text(size = 22),
        legend.spacing.x = unit(1.0, 'cm'), 
        legend.position = "none",
        legend.justification = c(-0.05, -0.2), 
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)) + 
  guides(fill = guide_legend(nrow = 1, keywidth=0.4, keyheight=0.1, default.unit="inch"))


simulated_prop_plots <- proportion_damaged_simulated_plot | proportion_unique_simulated_plot 
ggsave(filename = file.path("/home/alicen/Projects/ReviewArticle/benchmark_results/comparison_metrics/proportions_simulated_plot.png"), 
       plot = simulated_prop_plots, width = 14, height = 12, dpi = 300)




#-------------------------------------------------------------------------------
# SIMILATIRY PCA  
#-------------------------------------------------------------------------------

# Similarity of tools in PCA  ----

# Calculate pairwise cohen's kappa scores for each dataset
calculate_similarity <- function(data){
  
  # Isolate columns of interest 
  data <- data[, c("ddqc", "DropletQC", "ensembleKQC", "miQC", "scater", "valiDrops",
                   "manual_all", "manual_mito_isolated", "manual_mito", "manual_mito_ribo", "manual_malat1")]
  
  # Convert damaged/cell labels to binary 1/0
  data <- as.data.frame(lapply(data, function(x) ifelse(x == 'damaged', 1, 0)))
  
  # Extract the number of tools (matrix dimensions defining)
  tools <- colnames(data)
  num_tools <- length(data)
  
  
  # Calculate cohen kappa ------
  
  # Initialize an empty matrix for Cohen's Kappa
  kappa_matrix <- matrix(NA, nrow = num_tools, ncol = num_tools, 
                         dimnames = list(tools, tools))
  
  for (i in 1:num_tools) {
    for (j in 1:num_tools) {
      kappa_result <- kappa2(data.frame(data[, i], data[, j]))$value
      kappa_matrix[i, j] <- kappa_result
    }
  }
  
  # Convert any NaNs in the matrix to 0
  kappa_matrix[is.nan(kappa_matrix)] <- 0
  
  # Name the rows & columns 
  tool_names <- c("ddqc", "DropletQC", "ensembleKQC", "miQC", "scater", "valiDrops",
  "manual_all", "manual_mito_isolated", "manual_mito", "manual_mito_ribo", "manual_malat1")

  rownames(kappa_matrix) <- colnames(kappa_matrix) <- tool_names
  
  return(kappa_matrix)
  
}

# Apply the calculation
# Ground truth 
apoptotic_matrix <- calculate_similarity(apoptotic)
pro_apoptotic_matrix <- calculate_similarity(pro_apoptotic)
GM18507_dead_matrix <- calculate_similarity(GM18507_dead)
GM18507_dying_matrix <- calculate_similarity(GM18507_dying)
PDX_matrix <- calculate_similarity(PDX)

# Non-ground truth
A549_matrix <- calculate_similarity(A549)
dLiver_matrix <- calculate_similarity(dLiver)
dLung_matrix <- calculate_similarity(dLung)
dPBMC_matrix <- calculate_similarity(dPBMC)
ductal_matrix <- calculate_similarity(ductal)
glio_matrix <- calculate_similarity(glio)
HCT116_matrix <- calculate_similarity(HCT116)
hLiver_matrix <- calculate_similarity(hLiver)
hLung_matrix <- calculate_similarity(hLung)
hodgkin_matrix <- calculate_similarity(hodgkin)
hPBMC_matrix <- calculate_similarity(hPBMC)
Jurkat_matrix <- calculate_similarity(Jurkat)
mLiver_matrix <- calculate_similarity(mLiver)
mLung_matrix <- calculate_similarity(mLung)
mPBMC_matrix <- calculate_similarity(mPBMC)

# Simulated
simulated_matrix <- list()
for (name in names(simulated)){
  
  simulated_matrix[[name]] <- calculate_similarity(simulated[[name]])
  
}


# Calculate the median similarity scores across input datasets
calculate_median <- function(input_list) {
  
  # Initialize an empty list to store the matrices
  matrices <- list()
  
  # Extract matrix
  for (item in input_list) {
    matrices <- append(matrices, list(item))
  }
  
  # Combine into array
  matrices_array <- array(unlist(matrices), dim = c(11, 11, length(matrices)))
  
  # Calculate the median for each position across the matrices
  result <- apply(matrices_array, c(1, 2), median, na.rm = TRUE)
  

  return(result) # 11 x 11 matrix 
  
}

# Define lists for input
groundtruth_matrices <- list(A549_matrix, dLiver_matrix, dLung_matrix, dPBMC_matrix, ductal_matrix,
                             glio_matrix, HCT116_matrix, hLiver_matrix, hLung_matrix, hodgkin_matrix,
                             hPBMC_matrix, Jurkat_matrix, mLiver_matrix, mLung_matrix, mPBMC_matrix)

non_groundtruth_matrices <- list(apoptotic_matrix, pro_apoptotic_matrix, GM18507_dead_matrix, GM18507_dying_matrix, PDX_matrix)

# Simulated
simulated_matrices <- list(simulated_matrix$control_sim_1_2.5, simulated_matrix$control_sim_1_5, simulated_matrix$control_sim_1_10, simulated_matrix$control_sim_1_15, simulated_matrix$control_sim_1_20,
                           simulated_matrix$control_sim_2_2.5, simulated_matrix$control_sim_2_5, simulated_matrix$control_sim_2_10, simulated_matrix$control_sim_2_15, simulated_matrix$control_sim_2_20,
                           simulated_matrix$control_sim_3_2.5, simulated_matrix$control_sim_3_5, simulated_matrix$control_sim_3_10, simulated_matrix$control_sim_3_15, simulated_matrix$control_sim_3_20,
                           simulated_matrix$stimulated_sim_1_2.5, simulated_matrix$stimulated_sim_1_5, simulated_matrix$stimulated_sim_1_10, simulated_matrix$stimulated_sim_1_15, simulated_matrix$stimulated_sim_1_20,
                           simulated_matrix$stimulated_sim_2_2.5, simulated_matrix$stimulated_sim_2_5, simulated_matrix$stimulated_sim_2_10, simulated_matrix$stimulated_sim_2_15, simulated_matrix$stimulated_sim_2_20,
                           simulated_matrix$stimulated_sim_3_2.5, simulated_matrix$stimulated_sim_3_5, simulated_matrix$stimulated_sim_3_10, simulated_matrix$stimulated_sim_3_15, simulated_matrix$stimulated_sim_3_20
)

# Find the median 
groundtruth_similarity <- calculate_median(groundtruth_matrices)
non_groundtruth_similarity <- calculate_median(non_groundtruth_matrices)
simulated_similarity <- calculate_median(simulated_matrices)

PlotSimilarity <- function(matrix, title, metric = "Cohen's Kappa") {
  
  # Add the tool names -----
  similarity_matrix <- matrix
  rownames(similarity_matrix) <- colnames(similarity_matrix) <- c("ddqc", "DropletQC", "ensembleKQC", "miQC", "scater", "valiDrops",
  "manual_all", "manual_mito_isolated", "manual_mito", "manual_mito_ribo", "manual_malat1")
  
  # Add a small constant to the diagonal elements to avoid zero variance issues
  diag(similarity_matrix) <- diag(similarity_matrix) + 1e-6
  
  # Create a PCA plot -----
  
  # Perform PCA on the similarity matrix
  pca_result <- prcomp(similarity_matrix, scale. = TRUE)
  
  # Extract PCA coordinates for the first two principal components
  pca_coords <- as.data.frame(pca_result$x[, 1:2])
  colnames(pca_coords) <- c("PC1", "PC2")
  
  # Add a column for the tool names
  pca_coords$tool <- rownames(similarity_matrix)
  

  # Plot without legend (clustered viewing)
  pca_plot <- ggplot(pca_coords, aes(x = PC1, y = PC2, fill = tool, label = tool)) +
    geom_point(size = 10, shape = 21, color = "white") +  
    # geom_text_repel(aes(label = tool, nudge_y = nudge_y), size = 8, color = "black", segment.color = "grey", segment.size = 0.5, direction = "both") +
    scale_fill_manual(values = strategy_colours) +
    theme_classic() +
    labs(title = title,
         x = "PC 1",
         y = "PC 2") +
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),  
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      axis.text = element_blank(),
      axis.ticks =  element_blank(),
      axis.title.y = element_text(size = 18, vjust = 0.5),
      axis.title.x = element_text(size = 18, vjust = -0.5),
      plot.title = element_text(hjust = -0.1, vjust = 3.8, face = "bold", size = 22), 
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(50, 50, 68, 50))
  

  # Return output -----
  return(pca_plot)
  
} 

# Run the function on each covariate case
groundtruth_similarity_PCA <- PlotSimilarity(matrix = groundtruth_similarity, title = "")
ggsave(filename = file.path("/home/alicen/Projects/ReviewArticle/benchmark_results/comparison_metrics/groundtruth_PCA.png"), 
       plot = groundtruth_similarity_PCA, width = 8, height = 7, dpi = 300)


non_groundtruth_similarity_PCA <- PlotSimilarity(non_groundtruth_similarity, title = "")
ggsave(filename = file.path("/home/alicen/Projects/ReviewArticle/benchmark_results/comparison_metrics/non_groundtruth_PCA.png"), 
       plot = non_groundtruth_similarity_PCA, width = 8, height = 7, dpi = 300)


simulated_similarity_PCA <- PlotSimilarity(simulated_similarity, title = "")
ggsave(filename = file.path("/home/alicen/Projects/ReviewArticle/benchmark_results/comparison_metrics/simulated_PCA.png"), 
       plot = simulated_similarity_PCA, width = 8, height = 7, dpi = 300)



#-------------------------------------------------------------------------------
# CONSISTENCY SCORES  
#-------------------------------------------------------------------------------

# Consistency scoring metrics  ----

# consistency = (weighting constant x similarity deviation) + (weighting constant x proportion damaged detection deviation

calculate_deviation_scores <- function(kappa_matrices, damaged_df) {
  
  
  tool_names <- c("ddqc", "DropletQC", "ensembleKQC", "miQC", "scater", "valiDrops", 
                  "manual_all", "manual_mito_isolated", "manual_mito", "manual_mito_ribo", "manual_malat1")

    # Initialize an empty data frame to store the standard deviation of Kappa scores for each tool
  deviation_scores <- data.frame(tool = character(),
                                 similarity_deviation = numeric(), 
                                 stringsAsFactors = FALSE)
  
  # For each strategy/tool (10 total)
  for (tool in 1:length(tool_names)) {
    
    pairwise_deviations <- numeric()  # To store deviations for each pairwise comparison
    
    # Calculate the pairwise standard deviations across datasets
    for (other_tool in 1:length(tool_names)) {
      if (tool != other_tool) {
        
        # Collect the Kappa scores between tool and other_tool across datasets
        kappa_scores <- sapply(kappa_matrices, function(mat) mat[tool, other_tool])
        
        # Calculate the standard deviation of these Kappa scores
        deviation <- sd(kappa_scores, na.rm = TRUE)
        pairwise_deviations <- c(pairwise_deviations, deviation)
      }
    }
    
    # Mean deviation for this tool (across all its pairwise tool comparisons)
    mean_deviation <- mean(pairwise_deviations, na.rm = TRUE)
    
    # Store the result in the data frame
    deviation_scores <- rbind(deviation_scores, data.frame(tool = tool_names[tool], similarity_deviation = mean_deviation))
  }
  
  # Remove repeated rows 
  deviation_scores <- unique(deviation_scores)
  
  # Deviations of proportion damaged detected for each tool across datasets 
  deviation_scores$damaged_deviation <- apply(groundtruth_df, 1, sd, na.rm = TRUE)
  
  # Rescale similarity deviation to [0, 1] range
  # sim_min <- min(deviation_scores$similarity_deviation, na.rm = TRUE)
  # sim_max <- max(deviation_scores$similarity_deviation, na.rm = TRUE)
  # deviation_scores$similarity_deviation_rescaled <- (deviation_scores$similarity_deviation - sim_min) / (sim_max - sim_min)
  
  # Rescale damaged deviation to [0, 1] range
  damaged_min <- min(deviation_scores$damaged_deviation, na.rm = TRUE)
  damaged_max <- max(deviation_scores$damaged_deviation, na.rm = TRUE)
  deviation_scores$damaged_deviation_rescaled <- (deviation_scores$damaged_deviation - damaged_min) / (damaged_max - damaged_min)
  
  # Combine the rescaled deviations to calculate the final consistency score for each tool
  # Weights for similarity deviation and damaged proportion deviation
  w_similarity <- 0.5
  w_damaged <- 0.5
  
  deviation_scores$consistency <- w_similarity * deviation_scores$similarity_deviation + 
    w_damaged * deviation_scores$damaged_deviation_rescaled
  
  deviation_scores$consistency <- (1 - deviation_scores$consistency)
  
  # Convert the 'tool' column to a factor with levels specified in the desired order
  deviation_scores$tool <- factor(deviation_scores$tool, levels = tool_names)
  
  return(deviation_scores)
}


groundtruth_kappa_matrices <- list(apoptotic_matrix, apoptotic_matrix, GM18507_dead_matrix, GM18507_dying_matrix, PDX_matrix)
non_groundtruth_kappa_matrices <- list(A549_matrix, dLiver_matrix, dLung_matrix, 
                                       dPBMC_matrix, ductal_matrix, glio_matrix, HCT116_matrix,
                                       hLiver_matrix, hLung_matrix, hodgkin_matrix, 
                                       hPBMC_matrix, Jurkat_matrix, mLiver_matrix,
                                       mLung_matrix, mPBMC_matrix 
)
simulated_kappa_matrices <- list(simulated_matrix$control_sim_1_2.5, simulated_matrix$control_sim_1_5, simulated_matrix$control_sim_1_10, simulated_matrix$control_sim_1_15, simulated_matrix$control_sim_1_20,
                           simulated_matrix$control_sim_2_2.5, simulated_matrix$control_sim_2_5, simulated_matrix$control_sim_2_10, simulated_matrix$control_sim_2_15, simulated_matrix$control_sim_2_20,
                           simulated_matrix$control_sim_3_2.5, simulated_matrix$control_sim_3_5, simulated_matrix$control_sim_3_10, simulated_matrix$control_sim_3_15, simulated_matrix$control_sim_3_20,
                           simulated_matrix$stimulated_sim_1_2.5, simulated_matrix$stimulated_sim_1_5, simulated_matrix$stimulated_sim_1_10, simulated_matrix$stimulated_sim_1_15, simulated_matrix$stimulated_sim_1_20,
                           simulated_matrix$stimulated_sim_2_2.5, simulated_matrix$stimulated_sim_2_5, simulated_matrix$stimulated_sim_2_10, simulated_matrix$stimulated_sim_2_15, simulated_matrix$stimulated_sim_2_20,
                           simulated_matrix$stimulated_sim_3_2.5, simulated_matrix$stimulated_sim_3_5, simulated_matrix$stimulated_sim_3_10, simulated_matrix$stimulated_sim_3_15, simulated_matrix$stimulated_sim_3_20
)

groundtruth_deviation_scores <- calculate_deviation_scores(groundtruth_kappa_matrices, groundtruth_df)
non_groundtruth_deviation_scores <- calculate_deviation_scores(non_groundtruth_kappa_matrices, non_groundtruth_df)
simulated_deviation_scores <- calculate_deviation_scores(simulated_kappa_matrices, simulated_df)


# Plot 
consistency_groundtruth_plot <- ggplot(groundtruth_deviation_scores, aes(x = tool, y = consistency, fill = tool)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = strategy_colours) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(0, 1)) +
  coord_flip() + 
  labs(title = "", y = "") + 
  barplot_theme() + 
  theme(plot.title = element_text(hjust = -0.3, size = 22),
        plot.margin = margin(50, 20, 90, 20),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 30, vjust = -0.7),
        plot.background = element_rect(fill = "white", color = NA)
        ) 
ggsave(filename = file.path("/home/alicen/Projects/ReviewArticle/benchmark_results/comparison_metrics/consistency_groundtruth_plot.png"), 
       plot = consistency_groundtruth_plot, width = 7, height = 12, dpi = 300)


consistency_non_groundtruth_plot <- ggplot(non_groundtruth_deviation_scores, aes(x = tool, y = consistency, fill = tool)) +
  geom_bar(stat = "identity") +
  coord_flip() + 
  scale_fill_manual(values = strategy_colours) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(0, 1)) +
  labs(title = "", y = "") + 
  barplot_theme() + 
  theme(plot.title = element_text(hjust = -0.3, size = 22),
        plot.margin = margin(50, 20, 90, 20), # top right bottom left 
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 30, vjust = -0.7),
        plot.background = element_rect(fill = "white", color = NA)
        ) 
ggsave(filename = file.path("/home/alicen/Projects/ReviewArticle/benchmark_results/comparison_metrics/consistency_non_groundtruth_plot.png"), 
       plot = consistency_non_groundtruth_plot , width = 7, height = 12, dpi = 300)

consistency_simulated_plot <- ggplot(simulated_deviation_scores, aes(x = tool, y = consistency, fill = tool)) +
  geom_bar(stat = "identity") +
  coord_flip() + 
  scale_fill_manual(values = strategy_colours) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(0, 1)) +
  labs(title = "", y = "") + 
  barplot_theme() + 
  theme(plot.title = element_text(hjust = -0.3, size = 22),
        plot.margin = margin(50, 20, 90, 20), # top right bottom left 
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 30, vjust = -0.7),
        plot.background = element_rect(fill = "white", color = NA)
  ) 

ggsave(filename = file.path("/home/alicen/Projects/ReviewArticle/benchmark_results/comparison_metrics/consistency_simulated_plot.png"), 
       plot = consistency_simulated_plot, width = 7, height = 12, dpi = 300)

consistency_groundtruth_plot | consistency_non_groundtruth_plot | consistency_simulated_plot


#-------------------------------------------------------------------------------
# COMBINE
#-------------------------------------------------------------------------------

# Combine plots ----

# # # Quick calculation of total cell numbers for labels 
# data_frames <- list(A549, dLiver, dLung, dPBMC, ductal, 
#                     glio, HCT116, hLiver, hLung, hodgkin, 
#                     hPBMC, Jurkat, mLiver, mLung, mPBMC,
#                     apoptotic, pro_apoptotic, GM18507_dead, GM18507_dying, PDX)
# 
# 
# # Total cell number
# total_entries <- 0
# for (df in data_frames) {
#   total_entries <- total_entries + dim(df)[1]
# }
# 
# # Simulated
# total_entries <- 0
# for (name in names(simulated)) {
#    total_entries <- total_entries + dim(simulated[[name]])[1]
# }
# 
# # Median cell number from all datasets 
# cell_numbers <- data.frame(cell_number = c(
#   dim(A549)[1], dim(dLiver)[1], dim(dLung)[1], dim(dPBMC)[1], dim(ductal)[1], 
#   dim(glio)[1], dim(HCT116)[1], dim(hLiver)[1], dim(hLung)[1], dim(hodgkin)[1], 
#   dim(hPBMC)[1], dim(Jurkat)[1], dim(mLiver)[1], dim(mLung)[1], dim(mPBMC)[1],
#   dim(apoptotic)[1], dim(pro_apoptotic)[1], dim(GM18507_dead)[1], dim(GM18507_dying)[1], dim(PDX)[1]))
# 
# 
# median_value <- median(cell_numbers$cell_number)
# print(median_value)




# Create titles
groundtruth_title <- ggdraw() + 
  draw_label("Ground truth datasets (n = 6, cells = 40 604)", 
             fontface = 'bold', 
             size = 24,
             x = 0, 
             hjust = -0.12) + 
  theme(plot.margin = margin(0, 0, -50, 2))

non_groundtruth_title <- ggdraw() + 
  draw_label("Non-ground truth datasets (n = 15, cells = 94 603)", 
             fontface = 'bold', 
             size = 24,
             x = 0, 
             hjust = -0.11) + 
  theme(plot.margin = margin(0, 0, -50, 2))

# Arrange the groundtruth plots
groundtruth_plots <- plot_grid(
  proportion_damaged_groundtruth_plot,
  proportion_unique_groundtruth_plot, 
  consistency_groundtruth_plot,
  groundtruth_similarity_PCA,
  ncol = 4, 
  rel_widths = c(1, 1, 1, 2.2)
)

# Combine the groundtruth title and plots
groundtruth_combined <- plot_grid(
  groundtruth_title, 
  groundtruth_plots, 
  ncol = 1, 
  rel_heights = c(0.1, 1)
)

# Arrange the non-groundtruth plots
non_groundtruth_plots <- plot_grid(
  proportion_damaged_non_groundtruth_plot,
  proportion_unique_non_groundtruth_plot, 
  consistency_non_groundtruth_plot,
  non_groundtruth_similarity_PCA,
  ncol = 4, 
  rel_widths = c(1, 1, 1, 2.2)
)

# Combine the non-groundtruth title and plots
non_groundtruth_combined <- plot_grid(
  non_groundtruth_title, 
  non_groundtruth_plots, 
  ncol = 1, 
  rel_heights = c(0.1, 1)
)

# Combine the groundtruth and non-groundtruth plots
final_plot <- plot_grid(
  groundtruth_combined, 
  non_groundtruth_combined, 
  ncol = 1, 
  rel_heights = c(1, 1)
) + theme(plot.background = element_rect(fill = "white", color = NA))

# Save the final plot as a PNG with appropriate dimensions
ggsave("/home/alicen/Projects/ReviewArticle/benchmark_results/comparison_metrics/comparison_metrics.png", plot = final_plot, width = 38, height = 18, units = "in")


### End 
