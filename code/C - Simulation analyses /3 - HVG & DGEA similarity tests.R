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
# 2. Cell type preservation 
# Using the cell type labels of the data,
#
# 
# 3. DEG correctness 
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
seurat_objects <- list(control_sim_1 = control_sim_1, 
                       control_sim_2 = control_sim_2, 
                       control_sim_3 = control_sim_3,
                       stimulated_sim_1 = stimulated_sim_1, 
                       stimulated_sim_2 = stimulated_sim_2, 
                       stimulated_sim_3 = stimulated_sim_3)

# Calculate and store variable features (vf) for each parent dataset
variable_features <- list()
for (name in names(seurat_objects)) {
  seurat_objects[[name]] <- NormalizeData(seurat_objects[[name]]) %>% FindVariableFeatures(nfeatures = 100)
  variable_features[[name]] <- VariableFeatures(seurat_objects[[name]])
}


# Function for comparing hvg of parent to tests ----

compare_sets <- function(parent_vf,   # Variable features of seurat object without damage 
                         test,        # seurat of the damaged case of interest (2.5, 5, 10, 15, 20)
                         project_name # string to save this column as 
){
  
  # Initialize empty set 
  results = data.frame(method = character(),
                       jaccard_index = numeric(),
                       stringsAsFactors = FALSE)
  
  # Define each damaged detection method (as exists in column of Seurat meta data)
  methods <-  c("ddqc", "DropletQC", "ensembleKQC", "miQC", "scater", "valiDrops",
                "manual_all", "manual_mito_ribo", "manual_mito", "manual_malat1", "manual_mito_isolated")
  
  for (method in methods) {
    
    # Identify & isolate true cells for the method
    cells <- rownames(test@meta.data)[test[[method]] == "cell"]
    matrix <- test@assays$RNA$counts[, cells]
    object <- CreateSeuratObject(counts = matrix)
    
    # Find variable features (vf) for isolated cells
    object <- NormalizeData(object) %>% FindVariableFeatures(nfeatures = 100)
    vf <- VariableFeatures(object)
    
    # Calculate the similarity (sm) of this vf set to that of the specified parent 
    sm <- calculate_jaccard(vf, parent_vf)
    
    # Create a new row with method & jaccard index (named after input string)
    new_row <- data.frame(method = method, stringsAsFactors = FALSE)
    new_row[[project_name]] <- sm  # Assign similarity value to the dynamically named column
    
    # Append the new row to results
    results <- rbind(results, new_row)
    
  } 

  return(results)
  
}

# Function to run the comparisons for a specific damage percentage & calculate median
run_compare_sets <- function(damage_percentage) {
  
  # Run for each of the 6 datasets (3 control and 3 stimulated)
  control_sim_1_sm <- compare_sets(variable_features$control_sim_1, 
                                   simulated[[paste0("control_sim_1_", damage_percentage, "_seurat")]], 
                                   paste0("control_sim_1_", damage_percentage))
  
  control_sim_2_sm <- compare_sets(variable_features$control_sim_2, 
                                   simulated[[paste0("control_sim_2_", damage_percentage, "_seurat")]], 
                                   paste0("control_sim_2_", damage_percentage))
  
  control_sim_3_sm <- compare_sets(variable_features$control_sim_3, 
                                   simulated[[paste0("control_sim_3_", damage_percentage, "_seurat")]], 
                                   paste0("control_sim_3_", damage_percentage))
  
  stimulated_sim_1_sm <- compare_sets(variable_features$stimulated_sim_1, 
                                      simulated[[paste0("stimulated_sim_1_", damage_percentage, "_seurat")]], 
                                      paste0("stimulated_sim_1_", damage_percentage))
  
  stimulated_sim_2_sm <- compare_sets(variable_features$stimulated_sim_2, 
                                      simulated[[paste0("stimulated_sim_2_", damage_percentage, "_seurat")]], 
                                      paste0("stimulated_sim_2_", damage_percentage))
  
  stimulated_sim_3_sm <- compare_sets(variable_features$stimulated_sim_3, 
                                      simulated[[paste0("stimulated_sim_3_", damage_percentage, "_seurat")]], 
                                      paste0("stimulated_sim_3_", damage_percentage))
  
  # Merge and find median across the 6 datasets 
  median_sim <- merge(control_sim_1_sm, control_sim_2_sm, by = "method", all = TRUE) %>%
    merge(control_sim_3_sm, by = "method", all = TRUE) %>%
    merge(stimulated_sim_1_sm, by = "method", all = TRUE) %>%
    merge(stimulated_sim_2_sm, by = "method", all = TRUE) %>%
    merge(stimulated_sim_3_sm, by = "method", all = TRUE)
  

  median_col_name <- paste0("median_", damage_percentage)
  median_sim[[median_col_name]] <- apply(median_sim[, -1], 1, median, na.rm = TRUE)
  
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
  dplyr::select(method, median_2.5) %>%
  full_join(median_sim_5 %>% dplyr::select(method, median_5), by = "method") %>%
  full_join(median_sim_10 %>% dplyr::select(method, median_10), by = "method") %>%
  full_join(median_sim_15 %>% dplyr::select(method, median_15), by = "method") %>%
  full_join(median_sim_20 %>% dplyr::select(method, median_20), by = "method")

write.csv(combined_medians, "/home/alicen/Projects/ReviewArticle/damage_perturbation/analysis_results/hvg_similarity_100.csv", quote = FALSE, row.names = FALSE)


#-------------------------------------------------------------------------------
# DEGs
#-------------------------------------------------------------------------------

# DEG F1 scoring ----


# Psuedo-bulk DEGs between cell type in control and stimulated case (5 cell types)
# Compute the F1 score for each set of cell-type specific DEGs 

# Calculate parent DEGs all together 
# Edit merge all together, column with dataset, add this to psuedobulking to get 3 rreplicates for each 
sim_1 <-  merge(control_sim_1, y = stimulated_sim_1, add.cell.ids = c("control", "stimulated"), project = "sim_1")
sim_2 <-  merge(control_sim_2, y = stimulated_sim_2, add.cell.ids = c("control", "stimulated"), project = "sim_2")
sim_3 <-  merge(control_sim_3, y = stimulated_sim_1, add.cell.ids = c("control", "stimulated"), project = "sim_3")

# Pseudobulk counts (average expression of cell types, not cells, are used)

sim_1_psuedo <- AggregateExpression(sim_1, assays = "RNA", return.seurat = T, group.by = c("orig.ident", "celltype"))
sim_2_psuedo <- AggregateExpression(sim_2, assays = "RNA", return.seurat = T, group.by = c("orig.ident", "celltype"))
sim_3_psuedo <- AggregateExpression(sim_3, assays = "RNA", return.seurat = T, group.by = c("orig.ident", "celltype"))

# Add meta data column with celltype & stimulated status & set this to the Idents
sim_1_psuedo$orig.ident <- gsub("-[0-9]+_", "_", sim_1_psuedo$orig.ident)
sim_2_psuedo$orig.ident <- gsub("-[0-9]+_", "_", sim_2_psuedo$orig.ident)
sim_3_psuedo$orig.ident <- gsub("-[0-9]+_", "_", sim_3_psuedo$orig.ident)

Idents(sim_1_psuedo) <- "orig.ident"
Idents(sim_2_psuedo) <- "orig.ident"
Idents(sim_3_psuedo) <- "orig.ident"

# DEG 
bulk.mono.de <- FindMarkers(object = sim_1_psuedo, 
                            ident.1 = "control_T", 
                            ident.2 = "stimulated_T",
                            test.use = "DESeq2")


table(sim_1$celltype)

# Calculate and store variable features (vf) for each parent dataset
variable_features <- list()
for (name in names(seurat_objects)) {
  seurat_objects[[name]] <- NormalizeData(seurat_objects[[name]]) %>% FindVariableFeatures(nfeatures = 100)
  variable_features[[name]] <- VariableFeatures(seurat_objects[[name]])
}

cell_types <- c("Monocyte", "DC", "T", "B", "NK")


