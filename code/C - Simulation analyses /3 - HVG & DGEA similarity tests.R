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
# 
# 2. DEG correctness 
# For differential gene expression testing, each dataset will be matched to its corresponding case, control ~ stimulated. 
# From here, the DEGs between the cells will be found and compared to the DEGs between its undamged "parent" datasets. For instance, 
#   - control_sim_1 vs stimulated_sim_1 -> DEGs, compared to 
#   - control_sim_1_2.5 vs stimulated_sim_2.5 -> DEGs. 



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
  seurat_objects[[name]] <- NormalizeData(seurat_objects[[name]]) %>% FindVariableFeatures(nfeatures = 100)
  variable_features[[name]] <- VariableFeatures(seurat_objects[[name]])
}



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
  methods <-  c("ddqc", "DropletQC", "ensembleKQC", "miQC", "scater", "valiDrops",
                "manual_all", "manual_mito_ribo", "manual_mito", "manual_malat1", "manual_mito_isolated")

  for (method in methods) {
    
    # Identify & isolate true cells for the method
    cells <- rownames(test@meta.data)[test[[method]] == "cell"]
    matrix <- test@assays$RNA$counts[, cells]
    object <- CreateSeuratObject(counts = matrix)
    
    # Also calculate the proportion damaged 
    proportion_damaged <- 1 - (length(cells) / dim(test)[2])
    
    # Find variable features (vf) for isolated cells
    object <- NormalizeData(object) %>% FindVariableFeatures(nfeatures = 100)
    vf <- VariableFeatures(object)
    
    # Retrieve parent vf 
    parent_vf <- variable_features[[paste0(project_name)]]
    
    # Calculate the similarity (sm) of this vf set to that of the specified parent 
    sm <- calculate_jaccard(vf, parent_vf)
      
    # Create a new row with method & jaccard index (named after input string)
    new_row <- data.frame(method = method, stringsAsFactors = FALSE)
    new_row[[project_name]] <- sm  # Assign similarity value to the dynamically named column
    new_row[[paste0(project_name, "_", "damaged")]] <- proportion_damaged
      
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
  

  # Identify columns that do not end in _damaged
  columns_to_median <- grep(".*_damaged$", colnames(median_sim), invert = TRUE, value = TRUE)
  columns_to_median <- columns_to_median[columns_to_median != "method"]
  
  # Calculate the median for these columns only
  median_col_name <- paste0("median_", damage_percentage)
  median_sim[[median_col_name]] <- apply(median_sim[, columns_to_median], 1, median, na.rm = TRUE)
  
  # Identify columns that end in _damaged
  columns_damaged <- grep(".*_damaged$", colnames(median_sim), value = TRUE)
  
  # Calculate the median for the _damaged columns
  median_sim[["median_proportion"]] <- apply(median_sim[, columns_damaged], 1, median, na.rm = TRUE)
  
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
  dplyr::select(method, median_2.5, median_proportion) %>%
  full_join(median_sim_5 %>% dplyr::select(method, median_5, median_proportion), by = "method") %>%
  full_join(median_sim_10 %>% dplyr::select(method, median_10, median_proportion), by = "method") %>%
  full_join(median_sim_15 %>% dplyr::select(method, median_15, median_proportion), by = "method") %>%
  full_join(median_sim_20 %>% dplyr::select(method, median_20, median_proportion), by = "method")

write.csv(combined_medians, "/home/alicen/Projects/ReviewArticle/damage_perturbation/analysis_results/hvg_similarity.csv", quote = FALSE, row.names = FALSE)


#-------------------------------------------------------------------------------
# DEGs
#-------------------------------------------------------------------------------

# F1 scoring ----

calculate_f1 <- function(set_A, set_B) {
  
  # Calculate confusion metrics
  TP <- length(intersect(set_A, set_B))
  FN <- length(setdiff(set_A, set_B))
  FP <- length(setdiff(set_B, set_A))
  
  # Calculate precision
  precision <- TP / (TP + FP)
  
  # Calculate recall
  recall <- TP / (TP + FN)
  
  # Calculate F1 score
  f1_score <- 2 * (precision * recall) / (precision + recall)
  
  return(f1_score)
}

# Psuedo-bulk DEGs between cell type in control and stimulated case (5 cell types)
# Compute the F1 score for each set of cell-type specific DEGs 

# Calculate parent DEGs ----
# Edit merge all together, column with dataset, add this to psuedobulking to get 3 replicates for each 
parents_merged <-  merge(control_sim_1, 
                         y = list(stimulated_sim_1, 
                                  control_sim_2,  stimulated_sim_2, 
                                  control_sim_3,  stimulated_sim_3), 
                        add.cell.ids = c("control_1", "stimulated_1", 
                                         "control_2", "stimulated_2",
                                         "control_3", "stimulated_3"), 
                        project = "sim")

# Pseudobulk counts (average expression of cell types, not cells, are used)
parents_merged$stim <- gsub("_[0-9]+", "", parents_merged$orig.ident)
parents_psuedo <- AggregateExpression(parents_merged, assays = "RNA", return.seurat = T, group.by = c("orig.ident", "celltype", "stim"))

# Add meta data column with celltype & stimulated status & set this to the Idents
parents_psuedo$celltype.stim <- paste0(parents_psuedo$celltype, "_", parents_psuedo$stim)
Idents(parents_psuedo) <- "celltype.stim"

# DEG 
sets <- list(
  list(name = "B", ident_1 = "B_control", ident_2 = "B_stimulated"), 
  list(name = "DC", ident_1 = "DC_control", ident_2 = "DC_stimulated"), 
  list(name = "Monocyte", ident_1 = "Monocyte_control", ident_2 = "Monocyte_stimulated"), 
  list(name = "NK", ident_1 = "NK_control", ident_2 = "NK_stimulated"), 
  list(name = "T", ident_1 = "T_control", ident_2 = "T_stimulated"))
  

parent_DEGs <- list()

for (set in sets){
  
  # Run DESeq2 
  bulk.deg <- FindMarkers(object = parents_psuedo, 
                          ident.1 = set$ident_1, 
                          ident.2 = set$ident_2,
                          test.use = "DESeq2")
  
  # Keep only significantly differentially expressed genes 
  bulk.deg <- subset(bulk.deg, (p_val_adj <= 0.05) & (p_val_adj != 0))
  bulk.deg <- rownames(bulk.deg)

  parent_DEGs[[set$name]] <- bulk.deg
  
}

# Calculate unfiltered damaged DEGs ----

compare_deg_sets <- function(damage_percent,  # Case of interest (2.5, 5, 10, 15)
                             output_path = "/home/alicen/Projects/ReviewArticle/damage_perturbation/analysis_results/"
){
  
  # Retrieve 6 datasets of interest 
  control_sim_1 <- simulated[[paste0("control_sim_1_", damage_percent, "_seurat")]]
  control_sim_2 <- simulated[[paste0("control_sim_2_", damage_percent, "_seurat")]]
  control_sim_3 <- simulated[[paste0("control_sim_3_", damage_percent, "_seurat")]]
  stimulated_sim_1 <- simulated[[paste0("stimulated_sim_1_", damage_percent, "_seurat")]]
  stimulated_sim_2 <- simulated[[paste0("stimulated_sim_2_", damage_percent, "_seurat")]]
  stimulated_sim_3 <- simulated[[paste0("stimulated_sim_3_", damage_percent, "_seurat")]]
  
  # Merge together for later aggregation 
  merged <-  merge(control_sim_1, 
                   y = list(stimulated_sim_1, 
                            control_sim_2,  stimulated_sim_2, 
                            control_sim_3,  stimulated_sim_3), 
                   add.cell.ids = c("control_1", "stimulated_1", 
                                    "control_2", "stimulated_2",
                                    "control_3", "stimulated_3"), 
                   project = paste0("test", damage_percent))
  
  # Ensure all counts are used 
  merged <- JoinLayers(merged)
  
  # Sample names into column
  merged$dataset <- sub("(_[^_]+).*", "\\1", rownames(merged@meta.data))
  merged$stim <- gsub("_[0-9]+", "", merged$dataset)

    
  # Identify & isolate true cells for the method
  cells <- rownames(merged@meta.data)[merged[[method]] == "cell"]
  meta_data <- merged@meta.data[cells,]
  matrix <- merged@assays$RNA$counts[, cells]
  matrix <- round(matrix)
  object <- CreateSeuratObject(counts = matrix)
  object@meta.data <-  meta_data
    
  # Create psuedobulk counts =
  object_psuedo <- AggregateExpression(object, assays = "RNA", return.seurat = T, group.by = c("celltype", "stim", "dataset"))
  
  # Add meta data column with celltype & stimulated status & set this to the Idents
  object_psuedo$celltype.stim <- paste0(object_psuedo$celltype, "_", object_psuedo$stim)
  Idents(object_psuedo) <- "celltype.stim"
    
  # Create a table of the counts of each identifier
  ident_counts <- table(Idents(object_psuedo))
    
  # DEG   
  sets <- list(
    list(name = "B", ident_1 = "B_control", ident_2 = "B_stimulated"), 
    list(name = "DC", ident_1 = "DC_control", ident_2 = "DC_stimulated"), 
    list(name = "Monocyte", ident_1 = "Monocyte_control", ident_2 = "Monocyte_stimulated"), 
    list(name = "NK", ident_1 = "NK_control", ident_2 = "NK_stimulated"), 
    list(name = "T", ident_1 = "T_control", ident_2 = "T_stimulated"))
    
  for (set in sets) {
      
      # Check if both identifiers are present in the data and have a size of at least 3
      if (set$ident_1 %in% names(ident_counts) & set$ident_2 %in% names(ident_counts) & 
          ident_counts[set$ident_1] >= 3 & ident_counts[set$ident_2] >= 3) {
        
        # Run DESeq2 
        bulk.deg <- FindMarkers(object = object_psuedo, 
                                ident.1 = set$ident_1, 
                                ident.2 = set$ident_2,
                                test.use = "DESeq2")
        
        # Keep only significantly differentially expressed genes 
        bulk.deg <- subset(bulk.deg, (p_val_adj <= 0.05) & (p_val_adj != 0))
        bulk.deg <- rownames(bulk.deg)
        
        object_DEGs[[set$name]] <- bulk.deg
        
      } else {
        
        cat("\n")
        message("Skipped ", set)
        cat("\n")
        
        next
      }
    }
    
    # Compare the DEGs of the test, filtered object to the parent DEGs for each cell type
    B_F1 <- calculate_f1(parent_DEGs$B, object_DEGs$B)
    DC_F1 <- calculate_f1(parent_DEGs$DC, object_DEGs$DC)
    Monocyte_F1 <- calculate_f1(parent_DEGs$Monocyte, object_DEGs$Monocyte)
    NK_F1 <- calculate_f1(parent_DEGs$NK, object_DEGs$NK)
    T_F1 <- calculate_f1(parent_DEGs$T, object_DEGs$T)
    
    median_F1 <- median(B_F1, DC_F1, Monocyte_F1,  NK_F1, T_F1)
  
  return(median_F1)
  
}

# Run for each case of interest 
uf_damaged_2.5 <- compare_deg_sets(damage_percent = "2.5")
uf_damaged_5 <- compare_deg_sets(damage_percent = "5")
uf_damaged_10 <- compare_deg_sets(damage_percent = "10") 
uf_damaged_15 <- compare_deg_sets(damage_percent = "15") 
uf_damaged_20 <- compare_deg_sets(damage_percent = "20") 


# Calculate DEGs of the test cases and compare to parent ----

compare_deg_sets <- function(damage_percent,  # Case of interest (2.5, 5, 10, 15)
                             output_path = "/home/alicen/Projects/ReviewArticle/damage_perturbation/analysis_results/"
){
  
  # Retrieve 6 datasets of interest 
  control_sim_1 <- simulated[[paste0("control_sim_1_", damage_percent, "_seurat")]]
  control_sim_2 <- simulated[[paste0("control_sim_2_", damage_percent, "_seurat")]]
  control_sim_3 <- simulated[[paste0("control_sim_3_", damage_percent, "_seurat")]]
  stimulated_sim_1 <- simulated[[paste0("stimulated_sim_1_", damage_percent, "_seurat")]]
  stimulated_sim_2 <- simulated[[paste0("stimulated_sim_2_", damage_percent, "_seurat")]]
  stimulated_sim_3 <- simulated[[paste0("stimulated_sim_3_", damage_percent, "_seurat")]]
  
  # Merge together for later aggregation 
  merged <-  merge(control_sim_1, 
                   y = list(stimulated_sim_1, 
                            control_sim_2,  stimulated_sim_2, 
                            control_sim_3,  stimulated_sim_3), 
                   add.cell.ids = c("control_1", "stimulated_1", 
                                    "control_2", "stimulated_2",
                                    "control_3", "stimulated_3"), 
                   project = paste0("test", damage_percent))
  
  # Ensure all counts are used 
  merged <- JoinLayers(merged)
  
  # Sample names into column
  merged$dataset <- sub("(_[^_]+).*", "\\1", rownames(merged@meta.data))
  merged$stim <- gsub("_[0-9]+", "", merged$dataset)
  
  
  # Filter based on method of interest ---- 
  
  # Initialize empty df
  results = data.frame(method = character(),
                       F1_score = numeric(),
                       stringsAsFactors = FALSE)
  
  # Define each damaged detection method (as exists in column of Seurat meta data)
  methods <-  c("ddqc", "DropletQC", "ensembleKQC", "miQC", "scater", "valiDrops",
                "manual_all", "manual_mito_ribo", "manual_mito", "manual_malat1", "manual_mito_isolated")
  
  # method = "manual_all"
  for (method in methods) {
    
    # Print new lines
    cat("\n\n")
    message(method)
    
    # Identify & isolate true cells for the method
    cells <- rownames(merged@meta.data)[merged[[method]] == "cell"]
    meta_data <- merged@meta.data[cells,]
    matrix <- merged@assays$RNA$counts[, cells]
    matrix <- round(matrix)
    object <- CreateSeuratObject(counts = matrix)
    object@meta.data <-  meta_data

    # Create psuedobulk counts =
    object_psuedo <- AggregateExpression(object, assays = "RNA", return.seurat = T, group.by = c("celltype", "stim", "dataset"))
    
    # Add meta data column with celltype & stimulated status & set this to the Idents
    object_psuedo$celltype.stim <- paste0(object_psuedo$celltype, "_", object_psuedo$stim)
    Idents(object_psuedo) <- "celltype.stim"
    
    # Create a table of the counts of each identifier
    ident_counts <- table(Idents(object_psuedo))
    
    # DEG 
    sets <- list(
      list(name = "B", ident_1 = "B_control", ident_2 = "B_stimulated"), 
      list(name = "DC", ident_1 = "DC_control", ident_2 = "DC_stimulated"), 
      list(name = "Monocyte", ident_1 = "Monocyte_control", ident_2 = "Monocyte_stimulated"), 
      list(name = "NK", ident_1 = "NK_control", ident_2 = "NK_stimulated"), 
      list(name = "T", ident_1 = "T_control", ident_2 = "T_stimulated"))
    
    for (set in sets) {
      
      # Check if both identifiers are present in the data and have a size of at least 3
      if (set$ident_1 %in% names(ident_counts) & set$ident_2 %in% names(ident_counts) & 
          ident_counts[set$ident_1] >= 3 & ident_counts[set$ident_2] >= 3) {
          
          # Run DESeq2 
          bulk.deg <- FindMarkers(object = object_psuedo, 
                                  ident.1 = set$ident_1, 
                                  ident.2 = set$ident_2,
                                  test.use = "DESeq2")
          
          # Keep only significantly differentially expressed genes 
          bulk.deg <- subset(bulk.deg, (p_val_adj <= 0.05) & (p_val_adj != 0))
          bulk.deg <- rownames(bulk.deg)
          
          object_DEGs[[set$name]] <- bulk.deg
          
      } else {
        
        cat("\n")
        message("Skipped ", set)
        cat("\n")
        
        next
      }
    }
    
    # Compare the DEGs of the test, filtered object to the parent DEGs for each cell type
    B_F1 <- calculate_f1(parent_DEGs$B, object_DEGs$B)
    DC_F1 <- calculate_f1(parent_DEGs$DC, object_DEGs$DC)
    Monocyte_F1 <- calculate_f1(parent_DEGs$Monocyte, object_DEGs$Monocyte)
    NK_F1 <- calculate_f1(parent_DEGs$NK, object_DEGs$NK)
    T_F1 <- calculate_f1(parent_DEGs$T, object_DEGs$T)
    
    median_F1 <- median(B_F1, DC_F1, Monocyte_F1,  NK_F1, T_F1)
    
    # Create a new row with method & F1 score 
    new_row <- data.frame(method = method, stringsAsFactors = FALSE)
    project_name <- paste0("damaged_", damage_percent)
    new_row[[project_name]] <- median_F1 
    
    # Append the new row to results
    results <- rbind(results, new_row)
    
  } 

  write.csv(results, 
            paste0(output_path, project_name, ".csv"), 
            quote = FALSE, 
            row.names = FALSE)
  
  return(results)
  
}

# Run for each case of interest 
damaged_2.5 <- compare_deg_sets(damage_percent = "2.5")
damaged_5 <- compare_deg_sets(damage_percent = "5")
damaged_10 <- compare_deg_sets(damage_percent = "10") 
damaged_15 <- compare_deg_sets(damage_percent = "15") 
damaged_20 <- compare_deg_sets(damage_percent = "20") 











