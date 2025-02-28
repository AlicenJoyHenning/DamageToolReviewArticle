# SCRIPT CONTEXT 
#
# We wanted to explore whether/ how well the damaged strategies mitigate the downstream effects of contaminating damaged cells. 
# To do so, we created simulated data and perturbed subsets of the data to resemble damaged cells. As a positive control, 
# we used the unperturbed datasets to define true sets of highly variable genes as well as true sets of significant DEGs. 
# Then, applying the damaged strategies on the perturbed datasets, we did the same for resulting damaged filtered data.
# Finally, the differences between the true and filtered sets of HVG and DEGs were quantified, with strategies having 
# the largest difference from the truth being less suited for mitigating damaged contamination.
# 
#
# This script contains the code to use the scDesign2 model generated using reference data and applying it to create simulated data : 
# 1. Define scDesign2 inputs 
# 2. Run scResign2 to create simulated datasets for the control and stimulated reference data
#     - Includes generating the counts 
#     - Renaming output to create Seurat object with relevant meta data columns 
#     - Calculating QC metrics so the Seurat object can enter straight into damaged strategy testing 
#     - Outputting images of each simulated dataset to see clustering and marker expression 
# 


#-------------------------------------------------------------------------------
# PREPARATIONS
#-------------------------------------------------------------------------------

# Install scDesign2
#library(devtools)
#devtools::install_github("JSB-UCLA/scDesign2")
#install.packages("ggpubr")

# Load all packages 
packages <- c("scDesign2", "cluster", "copula", "dplyr", 
              "reshape2", "ggdendro", "gridExtra", "Seurat", "ggpubr", 
              "gtable", "cowplot", "ggplot2")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    library(pkg)
  }
}

#-------------------------------------------------------------------------------
# Retrieve reference data (again)
#-------------------------------------------------------------------------------

# Sample sheet (guide)
sample_sheet <- read.csv("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/sample_sheet.csv")
sample_sheet$Origin_Condition <- paste0(sample_sheet$Origin, sample_sheet$Condition)
sample_sheet$Name <- paste0(sample_sheet$Origin_Condition, "_R", sample_sheet$Replicate, "_DL", sample_sheet$Damage_level, "_CT", sample_sheet$Celltype_number, "_CN", sample_sheet$Cell_number)
sample_sheet <- sample_sheet %>% distinct(Name, .keep_all = TRUE)

# Produce one plot per replicate to check quality
sample_sheet$Plot <- ifelse(
  grepl("CN2000", sample_sheet$Name) & 
    grepl("R1", sample_sheet$Name), 
  "TRUE", 
  "FALSE"
)

# Core scRNAseq reference datasets of origin & number of cell types included 
# PH : PBMC High responders 
# PL : PBMC Low responders 
# CC : Cell line COVID19 
# CH : Cell line healthy controls

# Matrices
PH_CT6_matrix <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/unfiltered_genes/models/PBMC_high_CT6_reference_matrix.rds")
PH_CT4_matrix <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/unfiltered_genes/models/PBMC_high_CT4_reference_matrix.rds")
PH_CT3_matrix <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/unfiltered_genes/models/PBMC_high_CT3_reference_matrix.rds")
PL_CT6_matrix <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/unfiltered_genes/models/PBMC_low_CT6_reference_matrix.rds")
PL_CT4_matrix <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/unfiltered_genes/models/PBMC_low_CT4_reference_matrix.rds")
PL_CT3_matrix <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/unfiltered_genes/models/PBMC_low_CT3_reference_matrix.rds")
CC_CT1_matrix <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/unfiltered_genes/models/Cellline_covid_CT1_reference_matrix.rds")
CC_CT2_matrix <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/unfiltered_genes/models/Cellline_covid_CT2_reference_matrix.rds")
CH_CT1_matrix <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/unfiltered_genes/models/Cellline_healthy_CT1_reference_matrix.rds")
CH_CT2_matrix <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/unfiltered_genes/models/Cellline_healthy_CT2_reference_matrix.rds")


# Checking gene & cell numbers 
dim(PH_CT6_matrix) #  32738 genes 3000 cells 
dim(PL_CT6_matrix) # 32738 genes 3000 cells
dim(CC_CT2_matrix) # 33538 genes 3000 cells 
dim(CH_CT1_matrix) # 33538 genes 3000 cells 

# Read in models (Fit scDesign2 models.R)
PH_CT6_model <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/unfiltered_genes/models/PBMC_high_CT6_model.rds")
PH_CT4_model <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/unfiltered_genes/models/PBMC_high_CT4_model.rds")
PH_CT3_model <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/unfiltered_genes/models/PBMC_high_CT3_model.rds")
PL_CT6_model <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/unfiltered_genes/models/PBMC_low_CT6_model.rds")
PL_CT4_model <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/unfiltered_genes/models/PBMC_low_CT4_model.rds")
PL_CT3_model <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/unfiltered_genes/models/PBMC_low_CT3_model.rds")
CC_CT1_model <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/unfiltered_genes/models/cellline_covid_CT1_model.rds")
CC_CT2_model <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/unfiltered_genes/models/cellline_covid_CT2_model.rds")
CH_CT1_model <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/unfiltered_genes/models/cellline_healthy_CT1_model.rds")
CH_CT2_model <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/unfiltered_genes/models/cellline_healthy_CT2_model.rds")


#-------------------------------------------------------------------------------
# Run scDesign2
#-------------------------------------------------------------------------------


# Simulate datasets using scDesign2----
# - Generate control matrices
# - Adds damage (simulate_damage)
# - Calculates QC metrics 
# - Plot output (if desired)

# Function to introduce damage perturbation
simulate_damage <- function(count_matrix, 
                            damage_proportion, 
                            organism = "Hsap",
                            lambda = 7
) {
  
  # Calculate the target number of damaged cells given proportion 
  total_cells <- ncol(count_matrix)
  damaged_cell_number <- round(total_cells * damage_proportion)
  
  # Damage cell selections must be distributed across cell types evenly
  cell_types <- as.factor(sub("_.*", "", colnames(count_matrix))) 
  cell_type_counts <- table(cell_types)
  damage_per_type <- round(cell_type_counts * (damage_proportion))
  
  # Adjust to ensure the total matches the target
  total_damaged <- sum(damage_per_type)
  
  # Adjust if necessary
  if (total_damaged != damaged_cell_number) {
    diff <- damaged_cell_number - total_damaged
    adj_indices <- sample(seq_along(damage_per_type), abs(diff), replace = TRUE)
    damage_per_type[adj_indices] <- damage_per_type[adj_indices] + sign(diff)
  }
  
  # Select damaged cells
  damaged_cell_selections <- unlist(lapply(names(damage_per_type), function(ct) {
    cells_of_type <- which(cell_types == ct)
    sample(cells_of_type, size = damage_per_type[ct], replace = FALSE)
  }))
  
  # Make the column names unique (no longer need the celltypes in isolation)
  colnames(count_matrix) <- paste0(colnames(count_matrix), "_", seq_len(dim(count_matrix)[2])) # Keep cell types but ensure uniqueness
  
  # Storage of damage levels for all cell barcodes for plotting later
  damage_label <- data.frame(barcode = colnames(count_matrix)[damaged_cell_selections], status = rep("damaged", length(damaged_cell_selections)))
  undamaged_cell_number_cells <- setdiff(seq_len(total_cells), damaged_cell_selections)
  undamaged_cell_number <- data.frame(barcode = colnames(count_matrix)[undamaged_cell_number_cells], status = rep("control", length(undamaged_cell_number_cells)))
  damage_label <- rbind(damage_label, undamaged_cell_number)
  
  # Retrieve genes corresponding to the organism of interest
  if (organism == "Hsap") {
    mito_pattern <- "^MT-"
    ribo_pattern <- "^(RPS|RPL)"
    nuclear <- c("FIRRE", "NEAT1","XIST", "MALAT1", "MIAT", "MEG3", "KCNQ1OT1", "HOXA11-AS", "FTX") 
  }
  
  # Isolate gene set indices (consistent across cells, not subsetting the matrix)
  mito_idx <- grep(mito_pattern, rownames(count_matrix), ignore.case = FALSE)
  nucl_idx <- which(rownames(count_matrix) %in% nuclear)
  mito_idx <- c(mito_idx, nucl_idx)
  ribo_idx <- grep(ribo_pattern, rownames(count_matrix), ignore.case = FALSE)
  other_idx <- setdiff(seq_len(nrow(count_matrix)), union(mito_idx, ribo_idx))
  
  # Initialize for storing modified counts
  new_matrix <- count_matrix
  
  # Define the number of Monte Carlo samples
  n_samples <- 100000 # tested 1000 10000 100000 & found no helpful shift in shape from increasing
  
  # Loop over the damaged cells and apply the reduction to non-mito genes.
  for (i in seq_along(damaged_cell_selections)) {
    #i = 100  # JUST FOR TESTING
    
    # Index in the count matrix for the cell
    cell <- damaged_cell_selections[i]
    
    # Compute initial gene sums in the cell
    M_i <- sum(count_matrix[mito_idx, cell])  # Mitochondrial
    R_i <- sum(count_matrix[ribo_idx, cell])  # Ribosomal
    O_i <- sum(count_matrix[other_idx, cell]) # Other (non-mito & non-ribo)
    T_i <- R_i + O_i  # Total non-mito counts (everything that must be reduced)
    
    # Monte Carlo sampling to estimate r
    r_samples <- runif(n_samples, 0.01, 0.7)  # 0.01 - 0.7 to 0.1 - 0.7 to 0.05 to 0.5 test if less intense
    
    # Compute A, B, and the absolute error for each sample
    A_values <- round(r_samples * R_i) / (M_i + r_samples * T_i)
    B_values <- M_i / (M_i + round(r_samples * T_i))
    exp_decay_values <- exp(-lambda * A_values)
    errors <- abs(B_values - exp_decay_values)
    
    # Select the best r (minimizing the error)
    best_r <- r_samples[which.min(errors)]
    
    # Apply reduction to non-mito genes (ribo and other genes) in this cell.
    non_mito <- c(ribo_idx, other_idx)
    new_matrix[non_mito, cell] <- round(best_r * count_matrix[non_mito, cell])
    
  }
  
  # isolated indices for correct ordering of damage status 
  matched_indices <- match(colnames(count_matrix), damage_label$barcode)
  
  # QC statistics for all cells
  qc_summary <- data.frame(
    Cell = colnames(count_matrix),
    Damaged_Status = damage_label$status[matched_indices],
    Original_Features = colSums(count_matrix != 0),
    New_Features = colSums(new_matrix != 0),
    Original_MitoProp = colSums(count_matrix[mito_idx, , drop = FALSE]) / colSums(count_matrix),
    New_MitoProp = colSums(new_matrix[mito_idx, , drop = FALSE]) / colSums(new_matrix),
    Original_RiboProp = colSums(count_matrix[ribo_idx, , drop = FALSE]) / colSums(count_matrix),
    New_RiboProp = colSums(new_matrix[ribo_idx, , drop = FALSE]) / colSums(new_matrix)
  )
  
  # Return a list containing the new count matrix and the QC summary.
  return(list(matrix = new_matrix, qc_summary = qc_summary))
  
}


simulate_matrices <- function(
    origin, # "PH" "PL" "CH" "CC"
    replicate, # 1 2 3 
    damage_level, # 2, 5, 10, 20, 50
    celltypes,  # 1, 2, 3, 6
    cellnumber, # 200, 500, 1000, 2000, 5000
    generate_plot,
    output_dir
){
  
  # Define inputs for scDesign
  project_name <- paste0(origin, "_R", replicate, "_DL", damage_level, "_CT", celltypes, "_CN", cellnumber)
  origin_celltype <- paste0(origin, "_CT", celltypes)
  
  # Dataset of origin related to counts & model
  if (origin_celltype == "PH_CT6"){
    count_matrix <- PH_CT6_matrix
    model <- PH_CT6_model
  }
  
  if (origin_celltype == "PH_CT4"){
    count_matrix <- PH_CT4_matrix
    model <- PH_CT4_model
  }
  
  if (origin_celltype == "PH_CT3"){
    count_matrix <- PH_CT3_matrix
    model <- PH_CT3_model
  }
  
  if (origin_celltype == "PL_CT6"){
    count_matrix <- PL_CT6_matrix
    model <- PL_CT6_model
  }
  
  if (origin_celltype == "PL_CT4"){
    count_matrix <- PL_CT4_matrix
    model <- PL_CT4_model
  }
  
  if (origin_celltype == "PL_CT3"){
    count_matrix <- PL_CT3_matrix
    model <- PL_CT3_model
  }
  
  if (origin_celltype == "CC_CT1"){
    count_matrix <- CC_CT1_matrix
    model <- CC_CT1_model
  }
  
  if (origin_celltype == "CC_CT2"){
    count_matrix <- CC_CT2_matrix
    model <- CC_CT2_model
  }
  
  if (origin_celltype == "CH_CT1"){
    count_matrix <- CH_CT1_matrix
    model <- CH_CT1_model
  }
  
  if (origin_celltype == "CH_CT2"){
    count_matrix <- CH_CT2_matrix
    model <- CH_CT2_model
  }
  
  
  # Replicate reproducibility 
  if (replicate == 1){
    set.seed(7)
  }
  
  if (replicate == 2){
    set.seed(77)
  }
  
  if (replicate == 3){
    set.seed(777)
  }
  
  # Cell type selection 
  if (celltypes == 6){
    cell_type_selection <- c("B", "DC", "Monocyte", "NK", "T", "pDC") 
  }
  
  if (celltypes == 4){
    cell_type_selection <- c("B", "Monocyte", "T" , "NK") 
  }
  
  if (celltypes == 3){
    cell_type_selection <- c("B", "Monocyte", "T") 
  }
  
  if (celltypes == 2){
    cell_type_selection <- c("CD8", "CD4") 
  }
  
  if (celltypes == 1){
    cell_type_selection <- c("CD4") 
  }
  
  # Proportion in which cell types exist 
  cell_type_proportion <- table(colnames(count_matrix))[cell_type_selection]
  
  # scDesign2 to create count matrices 
  sim_matrix <- simulate_count_scDesign2(model_params = model, 
                                         n_cell_new = cellnumber, 
                                         sim_method = 'copula',
                                         cell_type_prop = cell_type_proportion)
  
  
  # Transfer gene names & assign cell identifiers 
  rownames(sim_matrix) <- rownames(count_matrix)
  celltypes <- colnames(sim_matrix)
  
  # Apply damage perturbation 
  perturbed <- simulate_damage(
    count_matrix = sim_matrix, 
    damage_proportion = (damage_level / 100)
  )
  
  # Extract damage level & add to cell annotations
  output <- perturbed$qc_summary
  output <- output[, c("Cell", "Damaged_Status")]
  perturbed_celltypes <- data.frame(celltypes)
  perturbed_celltypes$Cell <- output$Cell
  perturbed_celltypes$Annotation <- output$Damaged_Status[match(output$Cell, perturbed_celltypes$Cell)]
  perturbed_celltypes$Perturbed_annotation <- ifelse(perturbed_celltypes$Annotation == "damaged", "damaged", perturbed_celltypes$celltypes)
  
  # Create Seurat object
  seurat <- CreateSeuratObject(counts = perturbed$matrix, assay = "RNA", project = project_name)
  seurat$celltype <-  perturbed_celltypes$Perturbed_annotation
  # damage status 
  
  # Calculate standard quality control metrics ----- 
  
  # Retrieve the corresponding annotations for the organism of interest 
  seurat <- NormalizeData(seurat)
  
  # Calculate the feature percentages and normalised expressions 
  seurat$mt.percent <- PercentageFeatureSet(
    object   = seurat,
    features = rownames(seurat)[grepl("^MT-", rownames(seurat))]
  )
  
  seurat$rb.percent <- PercentageFeatureSet(
    object   = seurat,
    features = rownames(seurat@assays$RNA)[grepl("^RPS|^RPL", rownames(seurat@assays$RNA))],
    assay    = "RNA"
  ) 

  
  # Calculate psuedo-nuclear fraction score for DropletQC 
  seurat$malat1 <- FetchData(seurat, vars = "MALAT1", layer = "data") # data is default 
  seurat$nf_malat1 <- FetchData(seurat, vars = "MALAT1", layer = "counts")
  seurat$nf_malat1 <- (seurat$malat1 - min(seurat$malat1)) / (max(seurat$malat1) - min(seurat$malat1))
  
  # Save final object (resembles output of 1.1 Preprocess for section B)
  saveRDS(seurat, file.path(output_dir, paste0(project_name, "_seurat.rds")))
  
  # Viewing the data ----
  if (generate_plot) {
    
    test <- FindVariableFeatures(seurat) %>% # Already done the normalisation step
      ScaleData() %>% 
      RunPCA() %>%
      FindNeighbors(dims = 1:30) %>% 
      FindClusters() %>% 
      RunUMAP(dims = 1:30)
    
    # Note that these markers applicable for immune/PBMCs and T cell line, not generalizable 
    colours <- c("T" = "#A6BEAE",
                 "Monocyte" = "#88A1CD",
                 "DC" = "#A799C9",
                 "B" = "#BDC5EE",
                 "NK" = "#9DBD73",
                 "pDC" = "#E7948C",
                 "CD4" = "#A6BEAE",
                 "CD8" =  "#BDC5EE",
                 "damaged" = "red"
    )
    
    clusters <- DimPlot(test, group.by = "celltype")  + 
      NoAxes() + ggtitle(project_name) + scale_color_manual(values = colours) + 
      theme(panel.border = element_rect(colour = "black"))
    
    # MArkers 
    if (origin %in% c("PH", "PL")){
      
      # PBMC CELLTYPES
      T_cells <- FeaturePlot(test, "CD3E", cols = c("#E5E1E1", "#88A1CD")) + NoAxes() + NoLegend() + theme(panel.border = element_rect(colour = "black"))
      B_cells <- FeaturePlot(test, "CD79A", cols = c("#E1E1E1", "#88A1CD")) + NoAxes() + NoLegend() + theme(panel.border = element_rect(colour = "black"))
      NK_cells <- FeaturePlot(test, "NKG7", cols = c("#E1E1E1", "#88A1CD")) + NoAxes() + NoLegend() + theme(panel.border = element_rect(colour = "black"))
      Mon_cells <- FeaturePlot(test, "CD14", cols = c("#E1E1E1", "#88A1CD")) + NoAxes() + NoLegend() + theme(panel.border = element_rect(colour = "black"))
      
      plot <- clusters + ((T_cells | B_cells) / (NK_cells | Mon_cells))
      
    } else {
      
      # T CELLTYPES
      CD3E <- FeaturePlot(test, "CD3E", cols = c("#E1E1E1", "#88A1CD")) + NoAxes() + NoLegend() + theme(panel.border = element_rect(colour = "black"))
      CD2 <- FeaturePlot(test, "CD2", cols = c("#E1E1E1", "#88A1CD")) + NoAxes() + NoLegend() + theme(panel.border = element_rect(colour = "black"))
      GNLY <- FeaturePlot(test, "GNLY", cols = c("#E1E1E1", "#88A1CD")) + NoAxes() + NoLegend() + theme(panel.border = element_rect(colour = "black"))
      GZMB <- FeaturePlot(test, "GZMB", cols = c("#E1E1E1", "#88A1CD")) + NoAxes() + NoLegend() + theme(panel.border = element_rect(colour = "black"))
      
      
      plot <- clusters + ((CD3E | CD2) / (GNLY | GZMB)) 
      
    }
    
    ggsave(filename = file.path(output_dir, paste0(project_name, "_markers.png")), 
           plot, 
           width = 8, 
           height = 4, 
           units = "in",
           dpi = 300)
  }
  
  return(seurat)
  
}



# Iterate through sample sheet extracting parameters from entries
simulated_matrices <- list()
for (sample in seq_len(nrow(sample_sheet))) {
  
  store <- simulate_matrices(
    origin = paste0(sample_sheet$Origin[sample], sample_sheet$Condition[sample]), 
    replicate = sample_sheet$Replicate[sample], 
    damage_level = sample_sheet$Damage_level[sample], 
    celltypes = sample_sheet$Celltype_number[sample], 
    cellnumber = sample_sheet$Cell_number[sample], 
    generate_plot = sample_sheet$Plot[sample], 
    output_dir = "/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/output_data/control_simulated_matrices/"
  )
  
  simulated_matrices[[sample]] <- store
}

### End 
