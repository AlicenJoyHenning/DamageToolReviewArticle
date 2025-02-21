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
#     - Calculating QC metrics so the Seurat object can enter straight into damaged stratergy testing 
#     - Outputting images of each simulated dataset to see clustering and marker expression 
# 


#-------------------------------------------------------------------------------
# PREPARATIONS
#-------------------------------------------------------------------------------

# Install scDesign2
#library(devtools)
#devtools::install_github("JSB-UCLA/scDesign2")
#install.packages("copula")

# Load all packages 
packages <- c("AnnotationHub", "scDesign2", "cluster", "copula", "dplyr", 
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
sample_sheet <- read.csv("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/unperturbed_data/sample_sheet.csv")
sample_sheet$Origin_Condition <- paste0(sample_sheet$Origin, sample_sheet$Condition)
sample_sheet$Name <- paste0(sample_sheet$Origin_Condition, "_R", sample_sheet$Replicate, "_SD", sample_sheet$Sequencing_depth, "_CT", sample_sheet$Celltype_number, "_CN", sample_sheet$Cell_number)
sample_sheet_subset <- sample_sheet %>% distinct(Name, .keep_all = TRUE)


# Produce one plot per replicate to check quality
sample_sheet_subset$Plot <- ifelse(
    grepl("CN2000", sample_sheet_subset$Name) & 
    grepl("R1", sample_sheet_subset$Name), 
  "TRUE", 
  "FALSE"
)

# Testing on smaller subsets
# sample_sheet_subset <- subset(sample_sheet_subset, Plot == "TRUE")

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
dim(CH_CT2_matrix) # 33538 genes 3000 cells 

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

simulate_matrices <- function(
    origin, # "PH" "PL" "CH" "CC"
    replicate, # 1 2 3 
    celltypes,  # 1, 2, 3, 6
    cellnumber, # 200, 500, 1000, 2000, 5000
    generate_plot,
    output_dir
){
  
  # Define inputs for scDesign
  project_name <- paste0(origin, "_R", replicate, "_CT", celltypes, "_CN", cellnumber)
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
  
  
  # Transfer gene names & extract cell barcodes 
  # dim(sim_matrix)[1] == dim(sim_matrix)[1] # Must be TRUE
  rownames(sim_matrix) <- rownames(count_matrix)
  celltypes <- colnames(sim_matrix)
  colnames(sim_matrix) <- paste0(colnames(sim_matrix), "_", seq_len(dim(sim_matrix)[2])) # Keep cell types but ensure uniqueness
  write.csv(sim_matrix, file.path(output_dir, paste0(project_name, "_matrix.csv")))
  
  # Testing the mito percentages 
  mito_genes <- grepl("^MT-", rownames(sim_matrix))
  percent_mito <- ((colSums(sim_matrix[mito_genes, , drop = FALSE])) / (colSums(sim_matrix))) * 100
  max(percent_mito)
  
  # Create Seurat object
  seurat <- CreateSeuratObject(counts = sim_matrix, assay = "RNA", project = project_name)
  seurat$celltype <- celltypes
  
  # Calculate standard quality control metrics ----- 
  
  # Retrieve the corresponding annotations for the organism of interest 
  # Define human genes 
  malat1 <- c("MALAT1")
  seurat$malat1 <- FetchData(seurat, vars = malat1)
  
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
  
  seurat$malat1.percent <- PercentageFeatureSet(
    object   = seurat,
    features = malat1,
    assay    = "RNA"
  ) 
  
  # Calculate psuedo-nuclear fraction score for DropletQC 
  seurat$nf_malat1 <- FetchData(seurat, vars = "MALAT1", layer = "counts")
  
  # min-max normalization: scales values so the min value maps to 0 and max maps to 1 (make the expression values more like nf scores)
  seurat$nf_malat1 <- (seurat$nf_malat1 - min(seurat$nf_malat1)) / (max(seurat$nf_malat1) - min(seurat$nf_malat1))
  
  # Save final object (resembles output of 1.1 Preprocess for section B)
  saveRDS(seurat, file.path(output_dir, paste0(project_name, "_seurat.rds")))
  
  # Viewing the data ----
  if (generate_plot) {
    
    test <- NormalizeData(seurat) %>%
      FindVariableFeatures() %>%
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
for (sample in seq_len(nrow(sample_sheet_subset))) {
  
  store <- simulate_matrices(
    origin = paste0(sample_sheet_subset$Origin[sample], sample_sheet_subset$Condition[sample]), 
    replicate = sample_sheet_subset$Replicate[sample], 
    celltypes = sample_sheet_subset$Celltype_number[sample], 
    cellnumber = sample_sheet_subset$Cell_number[sample], 
    generate_plot = sample_sheet_subset$Plot[sample], 
    output_dir = "/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/unperturbed_data/control_simulated_matrices/"
  )
  
  simulated_matrices[[sample]] <- store
}

### End 
