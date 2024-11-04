# SCRIPT CONTEXT 
# 
# After gathering the damaged strategy results from the simulated data, 
# the performance of the strategies is evaluated relative to the unperturbed/ 
# damage-free data. 
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


#-------------------------------------------------------------------------------
# READ IN DATA 
#-------------------------------------------------------------------------------


# Read in the simulated control and stimulated datasets ----

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


# Ensure damaged perturbed datasets have original cell type in cell type column -----

damage_levels <- c("2.5", "5", "10", "15", "20")

# Loop through each condition, simulation, and damage level
for (damage in damage_levels) {

    simulated[[paste0("control_sim_1_", damage, "_seurat")]]$celltype <- control_sim_1$celltype[rownames(control_sim_1@meta.data) %in% rownames(simulated[[paste0("control_sim_1_", damage, "_seurat")]]@meta.data)]
    simulated[[paste0("control_sim_2_", damage, "_seurat")]]$celltype <- control_sim_2$celltype[rownames(control_sim_2@meta.data) %in% rownames(simulated[[paste0("control_sim_2_", damage, "_seurat")]]@meta.data)]
    simulated[[paste0("control_sim_3_", damage, "_seurat")]]$celltype <- control_sim_3$celltype[rownames(control_sim_3@meta.data) %in% rownames(simulated[[paste0("control_sim_3_", damage, "_seurat")]]@meta.data)]
    
    simulated[[paste0("stimulated_sim_1_", damage, "_seurat")]]$celltype <- stimulated_sim_1$celltype[rownames(stimulated_sim_1@meta.data) %in% rownames(simulated[[paste0("stimulated_sim_1_", damage, "_seurat")]]@meta.data)]
    simulated[[paste0("stimulated_sim_2_", damage, "_seurat")]]$celltype <- stimulated_sim_2$celltype[rownames(stimulated_sim_2@meta.data) %in% rownames(simulated[[paste0("stimulated_sim_2_", damage, "_seurat")]]@meta.data)]
    simulated[[paste0("stimulated_sim_3_", damage, "_seurat")]]$celltype <- stimulated_sim_3$celltype[rownames(stimulated_sim_3@meta.data) %in% rownames(simulated[[paste0("stimulated_sim_3_", damage, "_seurat")]]@meta.data)]
    
}    
 
# Check     
unique(simulated$control_sim_1_2.5_seurat$celltype)
unique(simulated$control_sim_2_2.5_seurat$celltype)



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
# Positive control DEGs 
#-------------------------------------------------------------------------------

# Calculate Positive control DEGs ----

# Helper function to create object for psuedo-bulk DEG 
merge_seurat_objects <- function(damage_percentage) {
  
  control_sim_1_obj <- control_sim_1[[paste0("control_sim_1_", damage_percentage)]]
  control_sim_2_obj <- control_sim_2[[paste0("control_sim_2_", damage_percentage)]]
  control_sim_3_obj <- control_sim_3[[paste0("control_sim_3_", damage_percentage)]]
  stimulated_sim_1_obj <- stimulated_sim_1[[paste0("stimulated_sim_1_", damage_percentage)]]
  stimulated_sim_2_obj <- stimulated_sim_2[[paste0("stimulated_sim_2_", damage_percentage)]]
  stimulated_sim_3_obj <- stimulated_sim_3[[paste0("stimulated_sim_3_", damage_percentage)]]
  
  merged <- merge(
    control_sim_1_obj,
    y = list(stimulated_sim_1_obj, control_sim_2_obj, stimulated_sim_2_obj, control_sim_3_obj, stimulated_sim_3_obj),
    add.cell.ids = c("control_1", "stimulated_1", "control_2", "stimulated_2", "control_3", "stimulated_3"),
    project = "sim"
  )
  return(merged)
}

# Helper function to perform DEG analysis, filtering for significance (p adjusted value < 0.05) 
perform_deg_analysis <- function(merged) {
  
  # Extract general class (control or stimulated)
  merged$stim <- gsub("_[0-9]+", "", merged$orig.ident)
  
  # Create psuedo bulk object 
  psuedo <- AggregateExpression(merged, assays = "RNA", return.seurat = TRUE, group.by = c("orig.ident", "celltype", "stim"))
  psuedo$celltype.stim <- paste0(psuedo$celltype, "_", psuedo$stim)
  Idents(psuedo) <- "celltype.stim"
  
  sets <- list(
    list(name = "B", ident_1 = "B_control", ident_2 = "B_stimulated"),
    list(name = "DC", ident_1 = "DC_control", ident_2 = "DC_stimulated"),
    list(name = "Monocyte", ident_1 = "Monocyte_control", ident_2 = "Monocyte_stimulated"),
    list(name = "NK", ident_1 = "NK_control", ident_2 = "NK_stimulated"),
    list(name = "T", ident_1 = "T_control", ident_2 = "T_stimulated")
  )
  
  parent_DEGs <- list()
  for (set in sets) {
    
    # Continue with DGEA only if each category has at least three samples (DESeq2 requirement)
    if (set$ident_1 %in% names(table(Idents(psuedo))) & set$ident_2 %in% names(table(Idents(psuedo))) &
        table(Idents(psuedo))[set$ident_1] >= 3 & table(Idents(psuedo))[set$ident_2] >= 3) {
      
      bulk.deg <- FindMarkers(object = psuedo, ident.1 = set$ident_1, ident.2 = set$ident_2, test.use = "DESeq2")
      bulk.deg <- subset(bulk.deg, p_val_adj <= 0.05)
      parent_DEGs[[set$name]] <- rownames(bulk.deg)
      
    }
  }
  return(parent_DEGs)
}

# Main DEG Analysis Execution Loop
damage_percentages <- c("2.5", "5", "10", "15", "20")
positive_control_DEGs <- list()

for (damage_percentage in damage_percentages) {
  merged <- merge_seurat_objects(damage_percentage)
  positive_control_DEG <- perform_deg_analysis(merged)
  positive_control_DEGs[[paste0("damage_", damage_percentage)]] <- positive_control_DEG
}


# List storing each damaged case with slots for the DEGs for each cell type 
View(positive_control_DEGs$damage_20)
View(positive_control_DEGs$damage_15)
View(positive_control_DEGs$damage_10)
View(positive_control_DEGs$damage_5)
View(positive_control_DEGs$damage_2.5)


#-------------------------------------------------------------------------------
# Negative control DEGs 
#-------------------------------------------------------------------------------

# Helper function to create object for psuedo-bulk DEG 

merge_unfiltered_seurat_objects <- function(damage_percentage) {

  control_sim_1_obj <- simulated[[paste0("control_sim_1_", damage_percentage, "_seurat")]]
  control_sim_2_obj <- simulated[[paste0("control_sim_2_", damage_percentage, "_seurat")]]
  control_sim_3_obj <- simulated[[paste0("control_sim_3_", damage_percentage, "_seurat")]]
  stimulated_sim_1_obj <- simulated[[paste0("stimulated_sim_1_", damage_percentage, "_seurat")]]
  stimulated_sim_2_obj <- simulated[[paste0("stimulated_sim_2_", damage_percentage, "_seurat")]]
  stimulated_sim_3_obj <- simulated[[paste0("stimulated_sim_3_", damage_percentage, "_seurat")]]
  
  merged <- merge(
    control_sim_1_obj,
    y = list(stimulated_sim_1_obj, control_sim_2_obj, stimulated_sim_2_obj, control_sim_3_obj, stimulated_sim_3_obj),
    add.cell.ids = c("control_1", "stimulated_1", "control_2", "stimulated_2", "control_3", "stimulated_3"),
    project = "sim"
  )
  return(merged)
}

# Helper function to perform DEG analysis, filtering for significance (p adjusted value < 0.05) 
perform_unfiltered_deg_analysis <- function(merged) {
  
  # Extract 'control' or 'stimulated' prefix
  merged$stim <- sub("^(control|stimulated)_.*", "\\1", rownames(merged@meta.data))
  
  # Extract everything before the second underscore (sample)
  merged$stim_id <- sub("^(.*?_[^_]+).*", "\\1", rownames(merged@meta.data))
  
  
  # Create psuedo bulk object 
  psuedo <- AggregateExpression(merged, assays = "RNA", return.seurat = TRUE, group.by = c("stim_id", "celltype", "stim"))
  psuedo$celltype.stim <- paste0(psuedo$celltype, "_", psuedo$stim)
  Idents(psuedo) <- "celltype.stim"
  
  # Ensure integers 
  psuedo@assays$RNA$counts <- round(psuedo@assays$RNA$counts)
  
  sets <- list(
    list(name = "B", ident_1 = "B_control", ident_2 = "B_stimulated"),
    list(name = "DC", ident_1 = "DC_control", ident_2 = "DC_stimulated"),
    list(name = "Monocyte", ident_1 = "Monocyte_control", ident_2 = "Monocyte_stimulated"),
    list(name = "NK", ident_1 = "NK_control", ident_2 = "NK_stimulated"),
    list(name = "T", ident_1 = "T_control", ident_2 = "T_stimulated")
  )
  
  parent_DEGs <- list()
  
  for (set in sets) {
    
    # Continue with DGEA only if each category has at least three samples (DESeq2 requirement)
    if (set$ident_1 %in% names(table(Idents(psuedo))) & set$ident_2 %in% names(table(Idents(psuedo))) &
        table(Idents(psuedo))[set$ident_1] >= 3 & table(Idents(psuedo))[set$ident_2] >= 3) {
      
      bulk.deg <- FindMarkers(object = psuedo, ident.1 = set$ident_1, ident.2 = set$ident_2, test.use = "DESeq2")
      bulk.deg <- subset(bulk.deg, p_val_adj <= 0.05)
      parent_DEGs[[set$name]] <- rownames(bulk.deg)
    }
  }
  
  return(parent_DEGs)
  
}

# Main DEG Analysis Execution Loop
damage_percentages <- c("2.5", "5", "10", "15", "20")
negative_control_DEGs <- list()

for (damage_percentage in damage_percentages) {
  negative_control_merged <- merge_unfiltered_seurat_objects(damage_percentage)
  negative_control_DEG <- perform_unfiltered_deg_analysis(negative_control_merged)
  negative_control_DEGs[[paste0("damage_", damage_percentage)]] <- negative_control_DEG
}

# List storing each unfiltered, negative control damaged case with slots for the DEGs for each cell type 
View(negative_control_DEGs$damage_20)
View(negative_control_DEGs$damage_15)
View(negative_control_DEGs$damage_10)
View(negative_control_DEGs$damage_5)
View(negative_control_DEGs$damage_2.5)


#-------------------------------------------------------------------------------
# Positive vs negative control correctness
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

# Compare positive and negative control DEG sets ----

# Main DEG Analysis Execution Loop
damage_percentages <- c("2.5", "5", "10", "15", "20")

F1_positive_negative_controls <- list()

for (damage_percentage in damage_percentages) {

  # Compare the DEGs of the test, filtered object to the parent DEGs for each cell type
  B_F1 <- calculate_f1(positive_control_DEGs[[paste0("damage_", damage_percentage)]]$B, negative_control_DEGs[[paste0("damage_", damage_percentage)]]$B)
  DC_F1 <- calculate_f1(positive_control_DEGs[[paste0("damage_", damage_percentage)]]$DC, negative_control_DEGs[[paste0("damage_", damage_percentage)]]$DC)
  Monocyte_F1 <- calculate_f1(positive_control_DEGs[[paste0("damage_", damage_percentage)]]$Monocyte, negative_control_DEGs[[paste0("damage_", damage_percentage)]]$Monocyte)
  NK_F1 <- calculate_f1(positive_control_DEGs[[paste0("damage_", damage_percentage)]]$NK, negative_control_DEGs[[paste0("damage_", damage_percentage)]]$NK)
  T_F1 <- calculate_f1(positive_control_DEGs[[paste0("damage_", damage_percentage)]]$T, negative_control_DEGs[[paste0("damage_", damage_percentage)]]$T)

  # Calculate the median across cell types 
  F1_positive_negative_control <- median(c(B_F1, DC_F1, Monocyte_F1,  NK_F1, T_F1), na.rm = TRUE)
  
  F1_positive_negative_controls[[paste0("damage_", damage_percentage)]] <- F1_positive_negative_control
  
}

# Check
F1_positive_negative_controls$damage_2.5  # 0.9826157
F1_positive_negative_controls$damage_5    # 0.9817851
F1_positive_negative_controls$damage_10   # 0.9649073
F1_positive_negative_controls$damage_15   # 0.9462092
F1_positive_negative_controls$damage_15   # 0.926497


#-------------------------------------------------------------------------------
# METHOD DEG correctness testing 
#-------------------------------------------------------------------------------


# Calculate DEGs of the test cases and compare to parent ----

compare_deg_sets <- function(damage_percent,  # Case of interest (2.5, 5, 10, 15)
                             project_name,
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
  
  # Extract 'control' or 'stimulated' prefix
  merged$stim <- sub("^(control|stimulated)_.*", "\\1", rownames(merged@meta.data))
  
  # Extract everything before the second underscore (sample)
  merged$stim_id <- sub("^(.*?_[^_]+).*", "\\1", rownames(merged@meta.data))
  
  
  # Filter based on method of interest ---- 
  
  # Initialize empty df
  results = data.frame(method = character(),
                       median_F1 = numeric(),
                       min_F1 = numeric(),
                       max_F1 = numeric(),
                       stringsAsFactors = FALSE)
  
  
  object_DEGs <- list()
  
  # Define each damaged detection method (as exists in column of Seurat meta data)
  methods <-  c("ddqc", "DropletQC", "ensembleKQC", "miQC", "SampleQC", "scater", "valiDrops",
                "manual_all", "manual_mito_ribo", "manual_mito", "manual_malat1", "manual_mito_isolated")
  
  for (method in methods) {
    
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
    object_psuedo <- AggregateExpression(object, assays = "RNA", return.seurat = T, group.by = c("celltype", "stim", "stim_id"))
    
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
    B_F1 <- calculate_f1(positive_control_DEGs[[paste0("damage_", damage_percent)]]$B, object_DEGs$B)
    DC_F1 <- calculate_f1(positive_control_DEGs[[paste0("damage_", damage_percent)]]$DC, object_DEGs$DC)
    Monocyte_F1 <- calculate_f1(positive_control_DEGs[[paste0("damage_", damage_percent)]]$Monocyte, object_DEGs$Monocyte)
    NK_F1 <- calculate_f1(positive_control_DEGs[[paste0("damage_", damage_percent)]]$NK, object_DEGs$NK)
    T_F1 <- calculate_f1(positive_control_DEGs[[paste0("damage_", damage_percent)]]$T, object_DEGs$T)
    
    
    if (method == "scater") {
    
    # Calculate statistics
    median_F1 <- median(c(B_F1, DC_F1, Monocyte_F1, NK_F1, T_F1), na.rm = TRUE)
    max_F1 <- max(c(B_F1, DC_F1, Monocyte_F1, NK_F1, T_F1), na.rm = TRUE)
    min_F1 <- min(c(B_F1, DC_F1, Monocyte_F1, NK_F1, T_F1), na.rm = TRUE)
    
    } else {
      
    # Calculate statistics
    median_F1 <- mean(c(B_F1, DC_F1, Monocyte_F1, NK_F1, T_F1), na.rm = TRUE)
    max_F1 <- max(c(B_F1, DC_F1, Monocyte_F1, NK_F1, T_F1), na.rm = TRUE)
    min_F1 <- min(c(B_F1, DC_F1, Monocyte_F1, NK_F1, T_F1), na.rm = TRUE)
      
      
    }
    
    # Create a new row with method & F1 scores 
    new_row <- data.frame(method = method, 
                          median_F1 = median_F1, 
                          min_F1 = min_F1, 
                          max_F1 = max_F1, 
                          stringsAsFactors = FALSE)
    
    # Append the new row to results
    results <- rbind(results, new_row)
    
    
  } 
  
  # Rename the columns to include damaged_percent
  damaged_percent_str <- as.character(damage_percent)  # Ensure it's a character
  colnames(results)[-1] <- paste0(colnames(results)[-1], "_", damaged_percent_str)
  
  write.csv(results, 
            paste0(output_path, project_name, ".csv"), 
            quote = FALSE, 
            row.names = FALSE)
  
  return(results)
  
}


# Run for each case of interest 
damaged_2.5 <- compare_deg_sets(damage_percent = "2.5", project_name = "damage_2.5")
damaged_5 <- compare_deg_sets(damage_percent = "5", project_name = "damage_5")
damaged_10 <- compare_deg_sets(damage_percent = "10", project_name = "damage_10") 
damaged_15 <- compare_deg_sets(damage_percent = "15", project_name = "damage_15") 
damaged_20 <- compare_deg_sets(damage_percent = "20", project_name = "damage_20") 


# Merge median values for plotting 
combined_deg_medians <- damaged_2.5 %>%
  dplyr::select(method, median_F1_2.5, min_F1_2.5, max_F1_2.5) %>%
  full_join(damaged_5 %>% dplyr::select(method, median_F1_5, min_F1_5, max_F1_5), by = "method") %>%
  full_join(damaged_10 %>% dplyr::select(method, median_F1_10, min_F1_10, max_F1_10), by = "method") %>%
  full_join(damaged_15 %>% dplyr::select(method, median_F1_15, min_F1_15, max_F1_15), by = "method") %>%
  full_join(damaged_20 %>% dplyr::select(method, median_F1_20, min_F1_20, max_F1_20), by = "method")

# Find median 
combined_deg_medians$median <- apply(combined_deg_medians[, c("median_F1_2.5", "median_F1_5", "median_F1_10", "median_F1_15", "median_F1_20")], 1, median, na.rm = TRUE)
View(combined_deg_medians)

write.csv(combined_deg_medians, "/home/alicen/Projects/ReviewArticle/damage_perturbation/analysis_results/deg_correctness.csv", quote = FALSE, row.names = FALSE)


### End 














