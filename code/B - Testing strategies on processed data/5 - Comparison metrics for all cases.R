# SCRIPT CONTEXT 
#
# Benchmarking damaged cell detection strategies. The labelled output objects are loaded at the start of the script. 
# (3 - Run remaining detection strategies.R) From here, the script uses the labels stored in the objects to calculate: 
#
# 1. The proportion of damaged cells for each tool - violin plot 
# 2. The similarity between the tools (Cohen's Kappa)-  PCA plots
# 3. Consistency score : combines proportion damaged & similarity - bar plot 
# 4. Dimensionality reduction plots - UMAPs or tSNEs using mt-rb genes 



#-------------------------------------------------------------------------------
# PREPARATIONS
#-------------------------------------------------------------------------------

# Load necessary libraries
library(ggplot2)
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


# Define palette
strategy_colours <- c("#999FA9","#586667","#73966E","#415741","#65719E","#383F54","#E8E5D6","#A4A388","#8B3860","#522239")


#-------------------------------------------------------------------------------
# PROPORTION DAMAGED 
#-------------------------------------------------------------------------------

# Proportion damaged as bar plot ----

find_proportion_damaged <- function(data, project_name){
  
  data <- data[, c("ddqc", "DropletQC", "ensembleKQC", "miQC", "valiDrops", 
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
mLung_df <- find_proportion_damaged(mLung, "mLung") #
mPBMC_df <- find_proportion_damaged(mPBMC, "mPBMC")
apoptotic_df <- find_proportion_damaged(apoptotic, "apoptotic")
pro_apoptotic_df <- find_proportion_damaged(pro_apoptotic, "pro_apoptotic")
GM18507_dead_df <- find_proportion_damaged(GM18507_dead, "GM18507_dead")
GM18507_dying_df <- find_proportion_damaged(GM18507_dying, "GM18507_dying")
PDX_df <- find_proportion_damaged(PDX, "PDX")

# List of data frames
non_groundtruth_dfs <- list(A549_df, dLiver_df, dLung_df, dPBMC_df, ductal_df, 
                            glio_df, HCT116_df, hLiver_df, hLung_df, hodgkin_df, 
                            hPBMC_df, Jurkat_df, mLiver_df, mLung_df, mPBMC_df)

groundtruth_dfs <- list(apoptotic_df, pro_apoptotic_df, 
                        GM18507_dead_df, GM18507_dying_df, 
                        PDX_df)


# Merge all data frames by 'strategy' column
non_groundtruth_df <- reduce(non_groundtruth_dfs, full_join, by = "strategy")
groundtruth_df <- reduce(groundtruth_dfs, full_join, by = "strategy")

# Find median values 
groundtruth_df <- groundtruth_df %>%
  rowwise() %>%
  mutate(
    lowest_value = min(c_across(apoptotic:PDX)),
    highest_value = max(c_across(apoptotic:PDX)),
    median_value = median(c_across(apoptotic:PDX))
  ) %>%
  ungroup()

non_groundtruth_df <- non_groundtruth_df %>%
  rowwise() %>%
  mutate(
    lowest_value = min(c_across(A549:mLung)),
    highest_value = max(c_across(A549:mLung)),
    median_value = median(c_across(A549:mLung))
  ) %>%
  ungroup()

# Plot as bar plots 

# Define theme 
proportion_damaged_theme <- function() {
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
      panel.border = element_rect(colour = "black", fill = NA),
      plot.title = element_text(vjust = 3.8, hjust = 0, face = "bold", size = 16),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(color = "black", size = 12, face = "bold"),
      plot.margin = margin(10, 10, 50, 10)
    )
}

# Precision plot 
proportion_damaged_groundtruth_plot <- ggplot(groundtruth_df, aes(x = strategy, y = median_value, fill = strategy)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = lowest_value, ymax = highest_value), width = 0.2) +
  scale_fill_manual(values = strategy_colours) +
  labs(title = "Damaged proportion", y = "") + 
  proportion_damaged_theme() + theme(plot.title = element_text(hjust = -0.05)) + ylim(0, 100)

proportion_damaged_non_groundtruth_plot <- ggplot(non_groundtruth_df, aes(x = strategy, y = median_value, fill = strategy)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = lowest_value, ymax = highest_value), width = 0.2) +
  scale_fill_manual(values = strategy_colours) +
  labs(title = "Damaged proportion", y = "") + 
  proportion_damaged_theme() + theme(plot.title = element_text(hjust = -0.05)) + ylim(0, 100)


proportion_plots <- (proportion_damaged_groundtruth_plot / proportion_damaged_non_groundtruth_plot)


#-------------------------------------------------------------------------------
# PROPORTION UNIQUE
#-------------------------------------------------------------------------------

# Proportion of uniquely labelled damaged cells per tool across datasets -----
data <- apoptotic
find_proportion_unique <- function(data, project_name, grountruth = "HEK293_apoptotic"){
  
  tools <- c("ddqc", "DropletQC", "ensembleKQC", "miQC", "valiDrops", 
             "manual_all", "manual_mito_isolated", "manual_mito", "manual_mito_ribo", "manual_malat1")
  
  
  if (!is.na(groundtruth)){
    
    # Define columns 
    columns <- c("orig.ident", "ddqc", "DropletQC", "ensembleKQC", "miQC", "valiDrops", 
               "manual_all", "manual_mito_isolated", "manual_mito", "manual_mito_ribo", "manual_malat1")
    
    data <- data[, columns]
    
    
    # Determine uniqueness by label from none or at most one other tool
    is_unique <- function(tool_label, other_labels) {
      if (tool_label == "damaged" && sum(other_labels == "damaged") <= 1) {
        return("unique")
      } else {
        return("common")
      }
    }
    
    # Determine if unique and correct
    is_unique_and_correct <- function(tool_label, other_labels, groundtruth) {
      if (tool_label == "damaged" && sum(other_labels == "damaged") <= 1 && tool_label == groundtruth) {
        return("correct")
      } else {
        return("incorrect")
      }
    }
    
    # Apply the functions to each tool
    for (tool in tools) {
      other_tools <- tools[tools != tool]
      data[[paste0(tool, "_unique")]] <- apply(data, 1, function(row) {
        is_unique(row[[tool]], row[other_tools])
      })
      data[[paste0(tool, "_correct")]] <- apply(data, 1, function(row) {
        is_unique_and_correct(row[[tool]], row[other_tools], row["groundtruth"])
      })
    }
    
    # Calculate the proportion of unique labels per tool
    unique_proportion <- sapply(tools, function(tool) {
      unique_count <- sum(data[[paste0(tool, "_unique")]] == "unique")
      total_count <- nrow(data)
      (unique_count / total_count) * 100
    })
    
    # Calculate the proportion of unique labels that are correct per tool
    correct_proportion <- sapply(tools, function(tool) {
      correct_count <- sum(data[[paste0(tool, "_correct")]] == "correct")
      unique_count <- sum(data[[paste0(tool, "_unique")]] == "unique")
      if (unique_count == 0) {
        return(0)
      } else {
        return((correct_count / unique_count) * 100)
      }
    })
    
    
    
    # Convert to data frame for better readability
    unique_proportion_df <- data.frame(strategy = names(unique_proportion), unique_proportion = unique_proportion)
    unique_proportion_df <- setNames(unique_proportion_df, c("strategy", project_name))
    
  }
  
  
  if (groundtruth == FALSE){
  
  data <- data[, tools]
  
  
  # Determine uniqueness by label from none or at most one other tool
  is_unique <- function(tool_label, other_labels) {
    if (tool_label == "damaged" && sum(other_labels == "damaged") == 0) {
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
  
  # Calculate damaged proportion for each column
  unique_proportion <- sapply(data, function(column) {
    unique_count <- sum(column == "unique")
    total_count <- length(column)
    (unique_count / total_count) * 100
  })
  
  # Convert to data frame for better readability
  unique_proportion_df <- data.frame(strategy = names(unique_proportion), unique_proportion = unique_proportion)
  unique_proportion_df <- setNames(unique_proportion_df, c("strategy", project_name))
  
  }
  
  
  
  return(unique_proportion_df)
  
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
apoptotic_udf <- find_proportion_unique(apoptotic, "apoptotic")
pro_apoptotic_udf <- find_proportion_unique(pro_apoptotic, "pro_apoptotic")
GM18507_dead_udf <- find_proportion_unique(GM18507_dead, "GM18507_dead")
GM18507_dying_udf <- find_proportion_unique(GM18507_dying, "GM18507_dying")
PDX_udf <- find_proportion_unique(PDX, "PDX")

# List of data frames
non_groundtruth_udfs <- list(A549_udf, dLiver_udf, dLung_udf, dPBMC_udf, ductal_udf, 
                            glio_udf, HCT116_udf, hLiver_udf, hLung_udf, hodgkin_udf, 
                            hPBMC_udf, Jurkat_udf, mLiver_udf, mLung_udf, mPBMC_udf)

groundtruth_udfs <- list(apoptotic_udf, pro_apoptotic_udf, 
                        GM18507_dead_udf, GM18507_dying_udf, 
                        PDX_udf)


# Merge all data frames by 'strategy' column
non_groundtruth_udf <- reduce(non_groundtruth_udfs, full_join, by = "strategy")
groundtruth_udf <- reduce(groundtruth_udfs, full_join, by = "strategy")

# Find median values 
groundtruth_udf <- groundtruth_udf %>%
  rowwise() %>%
  mutate(
    lowest_value = min(c_across(apoptotic:PDX)),
    highest_value = max(c_across(apoptotic:PDX)),
    median_value = median(c_across(apoptotic:PDX))
  ) %>%
  ungroup()

non_groundtruth_udf <- non_groundtruth_udf %>%
  rowwise() %>%
  mutate(
    lowest_value = min(c_across(A549:mLung)),
    highest_value = max(c_across(A549:mLung)),
    median_value = median(c_across(A549:mLung))
  ) %>%
  ungroup()

# Using same theme 

# Precision plot 
proportion_unique_groundtruth_plot <- ggplot(groundtruth_udf, aes(x = strategy, y = median_value, fill = strategy)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = lowest_value, ymax = highest_value), width = 0.2) +
  scale_fill_manual(values = strategy_colours) +
  labs(title = "Unique proportion", y = "") + 
  proportion_damaged_theme() + theme(plot.title = element_text(hjust = -0.05)) + ylim(0, 100)

proportion_unique_non_groundtruth_plot <- ggplot(non_groundtruth_udf, aes(x = strategy, y = median_value, fill = strategy)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = lowest_value, ymax = highest_value), width = 0.2) +
  scale_fill_manual(values = strategy_colours) +
  labs(title = "Unique proportion", y = "") + 
  proportion_damaged_theme() + theme(plot.title = element_text(hjust = -0.05)) + ylim(0, 100)


unique_proportion_plots <- (proportion_unique_groundtruth_plot / proportion_unique_non_groundtruth_plot)



#-------------------------------------------------------------------------------
# SIMILATIRY PCA  
#-------------------------------------------------------------------------------

# Similarity of tools in PCA  ----

# Calculate pairwise cohen's kappa scores for each dataset
calculate_similarity <- function(data){
  
  # Isolate columns of interest 
  data <- data[, c("ddqc", "DropletQC", "ensembleKQC", "miQC", "valiDrops", 
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
  
  return(kappa_matrix)
  
}

# Apply the calculation
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
apoptotic_matrix <- calculate_similarity(apoptotic)
pro_apoptotic_matrix <- calculate_similarity(pro_apoptotic)
GM18507_dead_matrix <- calculate_similarity(GM18507_dead)
GM18507_dying_matrix <- calculate_similarity(GM18507_dying)
PDX_matrix <- calculate_similarity(PDX)



# Calculate the median similarity scores across input datasets
calculate_median <- function(input_list) {
  
  # Initialize an empty list to store the matrices
  matrices <- list()
  
  # Extract matrix
  for (item in input_list) {
    matrices <- append(matrices, list(item))
  }
  
  # Combine into array
  matrices_array <- array(unlist(matrices), dim = c(10, 10, length(matrices)))
  
  # Calculate the median for each position across the matrices
  result <- apply(matrices_array, c(1, 2), median, na.rm = TRUE)
  
  
  
  return(result) # 10 x 10 matrix 
}

# Define lists for input
groundtruth_matrices <- list(A549_matrix, dLiver_matrix, dLung_matrix, dPBMC_matrix, ductal_matrix,
                             glio_matrix, HCT116_matrix, hLiver_matrix, hLung_matrix, hodgkin_matrix,
                             hPBMC_matrix, Jurkat_matrix, mLiver_matrix, mLung_matrix, mPBMC_matrix)

non_groundtruth_matrices <- list(apoptotic_matrix, pro_apoptotic_matrix, GM18507_dead_matrix, GM18507_dying_matrix, PDX_matrix)


# Find the median 
groundtruth_similarity <- calculate_median(groundtruth_matrices)
non_groundtruth_similarity <- calculate_median(non_groundtruth_matrices)

PlotSimilarity <- function(matrix, title, metric = "Cohen's Kappa") {
  
  # Add the tool names -----
  similarity_matrix <- matrix
  rownames(similarity_matrix) <- colnames(similarity_matrix) <- c("ddqc", "DropletQC", "ensembleKQC", "miQC", "valiDrops", 
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
  
  if (title == "Lung") {
    
    # Create PCA plot with legend & inverted axes 
    pca_coords_rotated <- pca_coords %>%
      mutate(PC1 = -PC1,  # Invert PC1
             PC2 = -PC2)  # Invert PC2
    
    
    # Create the PCA plot using ggplot2 with rotated points
    pca_plot <- ggplot(pca_coords_rotated, aes(x = PC1, y = PC2, fill = tool, label = tool)) +
      geom_point(size = 5, shape = 21, color = "black") +  # Use shape 21 for filled circles with an outline
      scale_fill_manual(values = tool_colors) +
      theme_classic() +
      labs(title = title,
           x = "PC 1",  # Keep the axis labels correct after rotation
           y = "PC 2") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            panel.background = element_rect(fill = "white", color = NA),  
            panel.border = element_rect(color = "black", fill = NA, size = 1),
            plot.background = element_rect(fill = "white", color = NA),
            axis.text = element_blank(),
            axis.ticks =  element_blank(),
            legend.position = "none")
    
    
  } 
  
  if (title == "Liver"){
    
    # Flip axes to match other PCA
    pca_coords_rotated <- pca_coords %>%
      mutate(PC1 = -PC1,  # Invert PC1
             PC2 = -PC2)  # Invert PC2
    
    
    # Create the PCA plot using ggplot2 with rotated points
    pca_plot <- ggplot(pca_coords_rotated, aes(x = PC1, y = PC2, fill = tool, label = tool)) +
      geom_point(size = 5, shape = 21, color = "black") +  # Use shape 21 for filled circles with an outline
      scale_fill_manual(values = tool_colors) +
      theme_classic() +
      labs(title = title,
           x = "PC 1",  # Keep the axis labels correct after rotation
           y = "PC 2") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            panel.background = element_rect(fill = "white", color = NA),  
            panel.border = element_rect(color = "black", fill = NA, size = 1),
            plot.background = element_rect(fill = "white", color = NA),
            axis.text = element_blank(),
            axis.ticks =  element_blank(),
            legend.position = "none")
    
  }
  
  # Plot without legend (clustered viewing)
  pca_plot <- ggplot(pca_coords, aes(x = PC1, y = PC2, fill = tool, label = tool)) +
      geom_point(size = 10, shape = 21, color = "white") +  
      scale_fill_manual(values = strategy_colours) +
      theme_classic() +
      labs(title = title,
           x = "PC 1",
           y = "PC 2") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            # legend.position = "bottom",
            legend.position = "none",
            legend.title = element_blank(),
            panel.background = element_rect(fill = "white", color = NA),  
            panel.border = element_rect(color = "black", fill = NA, size = 1),
            axis.text = element_blank(),
            axis.ticks =  element_blank(),
            plot.background = element_rect(fill = "white", color = NA)) 
  

  # Return output -----
  return(pca_plot)
  
}

# Run the function on each covariate case
groundtruth_similarity_PCA <- PlotSimilarity(groundtruth_similarity, title = "Cohen's Kappa")
non_groundtruth_similarity_PCA <- PlotSimilarity(non_groundtruth_similarity, title = "Cohen's Kappa")

similarity_plots <- groundtruth_similarity_PCA / non_groundtruth_similarity_PCA
proportion_plots | similarity_plots

#-------------------------------------------------------------------------------
# CONSISTENCY SCORES  
#-------------------------------------------------------------------------------

# Consistency scores as bar plot  ----




#-------------------------------------------------------------------------------
# VISUALISE LABELS 
#-------------------------------------------------------------------------------

# Compile UMAP plots with labels ----




