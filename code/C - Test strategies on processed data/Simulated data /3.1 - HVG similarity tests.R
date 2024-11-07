# SCRIPT CONTEXT 
# 
# After gathering the damaged strategy results from the simulated data, 
# the performance of the strategies is evaluated relative to the unperturbed/ 
# damage-free data. 
#
# 1. HVG similarity 
# The set of highly variable genes that define the damage-free data will be compared 
# to those of the same data that has undergone filtering by the damaged strategies after 
# been perturbed to include damaged. For instance,  
#   - HVGs of replicate 1 of the control data will be compared to the data
#     resulting from the filtered output of the 11 damaged strategies. This is repeated for each 
#     level of damage that was added (2.5 %, 5 %, 10 %, 15 %, and 20 %). Specifically: 
#     control_sim_1 ~ control_sim_1_2.5, control_sim_1_5, control_sim_1_10, control_sim_1_15, control_sim_1_20 




#-------------------------------------------------------------------------------
# PREPARATIONS
#-------------------------------------------------------------------------------

# Load libraries -------

packages <- c("cowplot", "devtools", "dplyr", "ggplot2", "glmGamPoi", "Matrix", "robustbase",
              "png", "Seurat", "tidyr", "SingleCellExperiment", "scuttle", "presto")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    library(pkg)
  }
}

# Read in the 6 unperturbed simulated ("parent") datasets ----

control_sim_1 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/simulated_references/control_sim_1_seurat.rds")
control_sim_2 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/simulated_references/control_sim_2_seurat.rds")
control_sim_3 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/simulated_references/control_sim_3_seurat.rds")
stimulated_sim_1 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/simulated_references/stimulated_sim_1_seurat.rds")
stimulated_sim_2 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/simulated_references/stimulated_sim_2_seurat.rds")
stimulated_sim_3 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/simulated_references/stimulated_sim_3_seurat.rds")

# Sizes for confirmation
dim(control_sim_1)[2]    # 4377
dim(control_sim_2)[2]    # 5981
dim(control_sim_3)[2]    # 5443
dim(stimulated_sim_1)[2] # 4279
dim(stimulated_sim_2)[2] # 4967
dim(stimulated_sim_3)[2] # 5463

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

# Confirm
# View(simulated$control_sim_1_2.5_seurat@meta.data)


# Filter controls -----

# Filter control cases to be more like reality, control = no damaged (small), damaged = added damage (larger), filtering = taking away from damage (reduction in damaged size, closer to control)

filter_controls <- function(control, project_name){
  
  filtered_controls <- list()
  
  percentages <- c(2.5, 5, 10, 15, 20)
  
  for (percent in percentages){
    
    # Retrieve matching perturbed object
    perturbed <- simulated[[paste0(project_name, "_", percent, "_seurat")]]
    unperturbed_cells <- rownames(perturbed@meta.data)[perturbed$orig.ident == "cell"]
    filtered <- subset(control, cells = unperturbed_cells) 
    
    # Add to list 
    filtered_controls[[paste0(project_name, "_", percent)]] <- filtered
                  
  }

  return(filtered_controls)
  
}

control_sim_1 <- filter_controls(control_sim_1, "control_sim_1")
control_sim_2 <- filter_controls(control_sim_2, "control_sim_2")
control_sim_3 <- filter_controls(control_sim_3, "control_sim_3")
stimulated_sim_1 <- filter_controls(stimulated_sim_1, "stimulated_sim_1")
stimulated_sim_2 <- filter_controls(stimulated_sim_2, "stimulated_sim_2")
stimulated_sim_3 <- filter_controls(stimulated_sim_3, "stimulated_sim_3")


#-------------------------------------------------------------------------------
# HVGs
#-------------------------------------------------------------------------------

# HVG similarity scoring ----
# Using the Jaccard Index, the HVG sets will be compared 

# Function to calculate the Jaccard Index (ratio of the intersection over union)
calculate_jaccard <- function(set_A, set_B){
  
  # Calculate intersection size
  intersection <- length(intersect(set_A, set_B))
  union <- length(union(set_A, set_B))
  
  # Calculate sizes of each set
  size_A <- length(set_A)
  size_B <- length(set_B)
  
  # Jaccard ----
  jaccard_index <- intersection / union
  
  return(jaccard_index)
  
}


# Parent HVGs ----

# List of parent, undamaged Seurat objects
seurat_objects <- list(
  
  # 2.5 % 
  control_sim_1_2.5 = control_sim_1$control_sim_1_2.5, 
  control_sim_2_2.5 = control_sim_2$control_sim_2_2.5, 
  control_sim_3_2.5 = control_sim_3$control_sim_3_2.5,
  stimulated_sim_1_2.5 = stimulated_sim_1$stimulated_sim_1_2.5, 
  stimulated_sim_2_2.5 = stimulated_sim_2$stimulated_sim_2_2.5, 
  stimulated_sim_3_2.5 = stimulated_sim_3$stimulated_sim_3_2.5, 
  
  # 5 % 
  control_sim_1_5 = control_sim_1$control_sim_1_5, 
  control_sim_2_5 = control_sim_2$control_sim_2_5, 
  control_sim_3_5 = control_sim_3$control_sim_3_5,
  stimulated_sim_1_5 = stimulated_sim_1$stimulated_sim_1_5, 
  stimulated_sim_2_5 = stimulated_sim_2$stimulated_sim_2_5, 
  stimulated_sim_3_5 = stimulated_sim_3$stimulated_sim_3_5, 
  
  # 10 % 
  control_sim_1_10 = control_sim_1$control_sim_1_10, 
  control_sim_2_10 = control_sim_2$control_sim_2_10, 
  control_sim_3_10 = control_sim_3$control_sim_3_10,
  stimulated_sim_1_10 = stimulated_sim_1$stimulated_sim_1_10, 
  stimulated_sim_2_10 = stimulated_sim_2$stimulated_sim_2_10, 
  stimulated_sim_3_10 = stimulated_sim_3$stimulated_sim_3_10, 
  
  # 15 % 
  control_sim_1_15 = control_sim_1$control_sim_1_15, 
  control_sim_2_15 = control_sim_2$control_sim_2_15, 
  control_sim_3_15 = control_sim_3$control_sim_3_15,
  stimulated_sim_1_15 = stimulated_sim_1$stimulated_sim_1_15, 
  stimulated_sim_2_15 = stimulated_sim_2$stimulated_sim_2_15, 
  stimulated_sim_3_15 = stimulated_sim_3$stimulated_sim_3_15,
  
  # 20 % 
  control_sim_1_20 = control_sim_1$control_sim_1_20, 
  control_sim_2_20 = control_sim_2$control_sim_2_20, 
  control_sim_3_20 = control_sim_3$control_sim_3_20,
  stimulated_sim_1_20 = stimulated_sim_1$stimulated_sim_1_20, 
  stimulated_sim_2_20 = stimulated_sim_2$stimulated_sim_2_20, 
  stimulated_sim_3_20 = stimulated_sim_3$stimulated_sim_3_20 
  
)

# Calculate and store variable features (vf) for each parent dataset
variable_features <- list()

for (name in names(seurat_objects)) {
  seurat_objects[[name]] <- NormalizeData(seurat_objects[[name]]) %>% FindVariableFeatures(nfeatures = 5000)
  variable_features[[name]] <- VariableFeatures(seurat_objects[[name]])
}


# Compare control (parent) to damaged (unfiltered) -----

# Define the control and stimulated simulations
simulations <- c("control_sim_1", "control_sim_2", "control_sim_3", 
                 "stimulated_sim_1", "stimulated_sim_2", "stimulated_sim_3")

# Define the damage percentages
damage_percentages <- c("2.5", "5", "10", "15", "20")

# Initialize an empty list to store the test cases
test_list <- list()

for (sim in simulations) {
  for (damage_percent in damage_percentages) {
    project_name <- paste0(sim, "_", damage_percent)
    test_list[[project_name]] <- list(
      test = simulated[[paste0(sim, "_", damage_percent, "_seurat")]],
      damaged_percent = damage_percent,
      project_name = project_name
    )
  }
}

# Convert the list to a format that matches the original structure
test_list <- unname(test_list)
  
compare_sets <- function(test,             # seurat of the damaged case of interest 
                         damaged_percent,  # 2.5, 5, 10, 15, or 20
                         project_name      # string to save this column as 
){
  # Initialize empty df
  results <- data.frame(case = character(),
                        jaccard_index = numeric(),
                        intersection = numeric(), 
                        control_unique = numeric(),
                        damage_unique = numeric(),
                        stringsAsFactors = FALSE)
  
  # Find variable features (vf) for isolated cells
  test <- NormalizeData(test) %>% FindVariableFeatures(nfeatures = 5000)
  vf <- VariableFeatures(test)
  
  # Retrieve parent vf 
  parent_vf <- variable_features[[paste0(project_name)]]
  
  # Calculate the overlap and unique genes to each set
  intersection <- length(intersect(vf, parent_vf))
  union <- length(union(vf, parent_vf))
  damage_unique <- length(vf) - intersection
  control_unique <- length(parent_vf) - intersection 
  
  # Find the actual genes unique to each case 
  damage_unique_HVGs <- setdiff(vf, parent_vf)
  control_unique_HVGs <- setdiff(parent_vf, vf)
  
  # Calculate the similarity (sm) of this vf set to that of the specified parent 
  sm <- calculate_jaccard(vf, parent_vf)
  
  # Create a new row with method, set info & jaccard index (named after input string)
  new_row <- data.frame(case = project_name, 
                        jaccard_index = sm,
                        intersection = intersection,
                        control_unique = control_unique,
                        damage_unique = damage_unique,
                        stringsAsFactors = FALSE)
  
  # Append the new row to results
  results <- rbind(results, new_row)
  
  # Return both results and HVGs
  return(list(results = results, 
              HVGs = list(damage_HVGs = damage_unique_HVGs, 
                          control_HVGs = control_unique_HVGs)))
}

# Run through all cases with a loop that appends results and HVGs
hvg_parent_v_unfiltered <- data.frame(case = character(),
                                      jaccard_index = numeric(),
                                      intersection = numeric(), 
                                      control_unique = numeric(),
                                      damage_unique = numeric(),
                                      stringsAsFactors = FALSE)


HVGs <- list()
for (item in test_list) {
  
  cat("\n")
  message(item$project_name)
  
  result <- compare_sets(test = item$test, 
                         damaged_percent = item$damaged_percent,
                         project_name = item$project_name)
  
  # Append the results to the dataframe
  hvg_parent_v_unfiltered <- rbind(hvg_parent_v_unfiltered, result$results)
  
  # Update the HVGs list
  HVGs[[item$project_name]] <- result$HVGs
  
}

# Add damage percent 
# Extract the number after the second underscore and add it to a new column
hvg_parent_v_unfiltered$damage <- sub(".*_.*_(\\d+)$", "\\1", hvg_parent_v_unfiltered$case)
hvg_parent_v_unfiltered$damage <- ifelse(!hvg_parent_v_unfiltered$damage %in% c(5, 10, 15, 20),
  sub(".*_.*_([0-9]+\\.[0-9]+)$", "\\1", hvg_parent_v_unfiltered$case), 
  hvg_parent_v_unfiltered$damage)

# Save outputs 
saveRDS(HVGs, "/home/alicen/Projects/ReviewArticle/damage_perturbation/analysis_results/HVGs_unique.rds")
write.csv(hvg_parent_v_unfiltered, 
          "/home/alicen/Projects/ReviewArticle/damage_perturbation/analysis_results/parent_v_unfiltered.csv",
          quote = FALSE, row.names = FALSE)




# Function for comparing hvg of parent to tests ----

compare_sets <- function(test,             # seurat of the damaged case of interest 
                         damaged_percent,  # 2.5, 5, 10, 15, or 20
                         project_name      # string to save this column as 
){
  
  # Initialize empty df
  results = data.frame(method = character(),
                       jaccard_index = numeric(),
                       propotion_damaged = numeric(),
                       stringsAsFactors = FALSE)
  
  # Define each damaged detection method (as exists in column of Seurat meta data)
  methods <-  c("ddqc", "DropletQC", "ensembleKQC", "miQC", "SampleQC", "scater", "valiDrops",
                "manual_all", "manual_mito_ribo", "manual_mito", "manual_malat1", "manual_mito_isolated")

  for (method in methods) {
    
    # Identify & isolate true cells for the method
    cells <- rownames(test@meta.data)[test[[method]] == "cell"]
    matrix <- test@assays$RNA$counts[, cells]
    object <- CreateSeuratObject(counts = matrix)
    
    # Also calculate the proportion damaged 
    proportion_damaged <- 1 - (length(cells) / dim(test)[2])
    
    # Find variable features (vf) for isolated cells
    object <- NormalizeData(object) %>% FindVariableFeatures(nfeatures = 5000)
    vf <- VariableFeatures(object)
    
    # Retrieve parent vf 
    parent_vf <- variable_features[[paste0(project_name)]]
    
    # Calculate the similarity (sm) of this vf set to that of the specified parent 
    sm <- calculate_jaccard(vf, parent_vf)
      
    # Create a new row with method, set info & jaccard index (named after input string)
    new_row <- data.frame(method = method, stringsAsFactors = FALSE)
    new_row[[project_name]] <- sm  
    new_row[[paste0(project_name, "_", "proportion")]] <- proportion_damaged
    
    # Append the new row to results
    results <- rbind(results, new_row)
    
  } 

  return(results)
  
}

# Function to run the comparisons for a specific damage percentage & calculate median
run_compare_sets <- function(damage_percentage) {

  # Run for each of the 6 datasets (3 control and 3 stimulated)
  control_sim_1_sm <- compare_sets(test = simulated[[paste0("control_sim_1_", damage_percentage, "_seurat")]], 
                                   damaged_percent = damage_percentage, 
                                   project_name = paste0("control_sim_1_", damage_percentage))
  
  control_sim_2_sm <- compare_sets(simulated[[paste0("control_sim_2_", damage_percentage, "_seurat")]], 
                                   damaged_percent = damage_percentage, 
                                   paste0("control_sim_2_", damage_percentage))
  
  control_sim_3_sm <- compare_sets(simulated[[paste0("control_sim_3_", damage_percentage, "_seurat")]], 
                                   damaged_percent = damage_percentage, 
                                   paste0("control_sim_3_", damage_percentage))
  
  stimulated_sim_1_sm <- compare_sets(simulated[[paste0("stimulated_sim_1_", damage_percentage, "_seurat")]], 
                                      damaged_percent = damage_percentage, 
                                      paste0("stimulated_sim_1_", damage_percentage))
  
  stimulated_sim_2_sm <- compare_sets(simulated[[paste0("stimulated_sim_2_", damage_percentage, "_seurat")]], 
                                      damaged_percent = damage_percentage, 
                                      paste0("stimulated_sim_2_", damage_percentage))
  
  stimulated_sim_3_sm <- compare_sets(simulated[[paste0("stimulated_sim_3_", damage_percentage, "_seurat")]], 
                                      damaged_percent = damage_percentage, 
                                      paste0("stimulated_sim_3_", damage_percentage))
  
  # Merge and find median across the 6 datasets 
  median_sim <- merge(control_sim_1_sm, control_sim_2_sm, by = "method", all = TRUE) %>%
    merge(control_sim_3_sm, by = "method", all = TRUE) %>%
    merge(stimulated_sim_1_sm, by = "method", all = TRUE) %>%
    merge(stimulated_sim_2_sm, by = "method", all = TRUE) %>%
    merge(stimulated_sim_3_sm, by = "method", all = TRUE)
  

  # Identify columns specific to the similarity index value
  columns_to_median <- grep(".*_proportion$", colnames(median_sim), invert = TRUE, value = TRUE)
  columns_to_median <- columns_to_median[columns_to_median != "method"]
  
  # Calculate the median for these columns only
  median_col_name <- paste0("median_", damage_percentage)
  median_sim[[median_col_name]] <- apply(median_sim[, columns_to_median], 1, median, na.rm = TRUE)
  
  # Identify columns that end in _damaged
  columns_damaged <- grep(".*_proportion$", colnames(median_sim), value = TRUE)
  
  # Calculate the median for the _damaged columns
  median_sim[["median_proportion"]] <- apply(median_sim[, columns_damaged], 1, median, na.rm = TRUE)
  
  # Calculate the lowest and highest values for columns not ending in "_proportion"
  lowest_col_name <- paste0("lowest_", damage_percentage)
  highest_col_name <- paste0("highest_", damage_percentage)
  
  # Add columns for the lowest and highest values (for plotting intervals)
  median_sim[[lowest_col_name]] <- apply(median_sim[, columns_to_median], 1, min, na.rm = TRUE)
  median_sim[[highest_col_name]] <- apply(median_sim[, columns_to_median], 1, max, na.rm = TRUE)
  
  
  return(median_sim)
  
}

# Run for each damage percentage
median_sim_2.5 <- run_compare_sets("2.5")
median_sim_5 <- run_compare_sets("5")
median_sim_10 <- run_compare_sets("10")
median_sim_15 <- run_compare_sets("15")
median_sim_20 <- run_compare_sets("20")

# Merge median values for plotting 
combined_medians <- median_sim_2.5 %>%
  dplyr::select(method, median_2.5, lowest_2.5, highest_2.5, median_proportion) %>%
  full_join(median_sim_5 %>% dplyr::select(method, median_5, lowest_5, highest_5, median_proportion), by = "method") %>%
  full_join(median_sim_10 %>% dplyr::select(method, median_10, lowest_10, highest_10, median_proportion), by = "method") %>%
  full_join(median_sim_15 %>% dplyr::select(method, median_15, lowest_15, highest_15, median_proportion), by = "method") %>%
  full_join(median_sim_20 %>% dplyr::select(method, median_20, lowest_20, highest_20, median_proportion), by = "method")

View(combined_medians)

write.csv(combined_medians, "/home/alicen/Projects/ReviewArticle/damage_perturbation/analysis_results/hvg_similarity.csv", quote = FALSE, row.names = FALSE)


### End 