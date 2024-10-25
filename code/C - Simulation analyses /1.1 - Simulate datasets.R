# SCRIPT CONTEXT 
#
# We wished to explore whether/ how well the damaged strategies mitigate the downstream effects 
# of contaminating damaged cells. To do so, we created simulated data (this script) and perturbed defined 
# subsets of the data to resemble damaged cells. As a positive control, we used the unperturbed datasets 
# to define true sets of highly variable genes as well as true sets of significant DEGs. Then, applying 
# the damaged strategies on the perturbed datasets, we did the same for resulting damaged filtered data.
# 
#
# This script contains the code to create simulated data using scDesign 2, a scRNA-seq simulator that 
# preserves gene names (required for all strategies), using from an annotated reference dataset from SeuratData
#   -  IFNB containing treated (IFN-STIMULATED) and control (IFN-CONTROL) cells


#-------------------------------------------------------------------------------
# PREPARATIONS
#-------------------------------------------------------------------------------

message("Loading libraries...")

# Install scDesign2
#library(devtools)
#devtools::install_github("JSB-UCLA/scDesign2")
#install.packages("copula")


# Load all packages 
packages <- c("scDesign2", "copula", "plyr", "reshape2", "gridExtra", 
              "Seurat", "SeuratData", "ggpubr", "cowplot", "ggplot2")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    library(pkg)
  }
}


#-------------------------------------------------------------------------------
# Retrieve & prepare reference data
#-------------------------------------------------------------------------------

message("Retrieve data...")

# SeuratData publicly available data ----
available_data <- AvailableData()
available_data <- subset(available_data, species == "human")

# IFNB-Stimulated & Control PBMCs
InstallData("ifnb")
data("ifnb")
ifnb <- UpdateSeuratObject(ifnb)

# Isolate stimulated and control cells (used to simulate datasets in isolation)
control <- subset(ifnb, stim == "CTRL")    # 6548 cells 
stimulated <- subset(ifnb, stim == "STIM") # 7451 cells 

# Extract counts 
control_matrix <-   as.matrix(control@assays$RNA@counts)
stimulated_matrix <-   as.matrix(stimulated@assays$RNA@counts)

# Recollect celltypes 
ifnb$celltypes <- "na"
ifnb$celltypes <- ifelse(ifnb$seurat_annotations %in% c("CD14 Mono",  "CD16 Mono"), "Monocyte", ifnb$celltypes)
ifnb$celltypes <- ifelse(ifnb$seurat_annotations %in% c("CD4 Memory T", "T activated" , "CD4 Naive T", "CD8 T"), "T", ifnb$celltypes)
ifnb$celltypes <- ifelse(ifnb$seurat_annotations %in% c("B Activated" , "B"), "B", ifnb$celltypes)
ifnb$celltypes <- ifelse(ifnb$seurat_annotations %in% c("pDC","DC"), "DC", ifnb$celltypes)
ifnb$celltypes <- ifelse(ifnb$seurat_annotations == "NK", "NK", ifnb$celltypes)


# Extract celltypes 
celltypes <- as.data.frame(ifnb@meta.data[, c("celltypes")]) # name column annotations
celltypes$cells <- rownames(ifnb@meta.data)
colnames(celltypes)[1] <- "celltype"


# Replace the colnames of the matrices with the corresponding cell types
celltype_map <- setNames(celltypes$celltype, celltypes$cells)
colnames(control_matrix) <- celltype_map[colnames(control_matrix)]
colnames(stimulated_matrix) <- celltype_map[colnames(stimulated_matrix)]

#-------------------------------------------------------------------------------
# Run scDesign2
#-------------------------------------------------------------------------------

# scDesign2 ----

# Cell types to include 
cell_type_selection <- c("Monocyte", "DC", "T", "B", "NK")   

# Proportion which these celltype exist 
cell_type_proportion_control <- table(colnames(control_matrix))[cell_type_selection]
cell_type_proportion_stimulated <- table(colnames(stimulated_matrix))[cell_type_selection]

# Read in models
control_model <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/control_model.rds")
simulated_control <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/stimulated_model.rds")


# Simulate datasets ----

# Controls
target_cells <- runif(1, min = 4000, max = 6000) # Randomly generates target cell number
control_sim_1 <- simulate_count_scDesign2(control_model, 
                                          n_cell_new = target_cells, 
                                          sim_method = 'copula',
                                          cell_type_prop = cell_type_proportion_control)
saveRDS(control_sim_1, "/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/control_sim_1.rds")


target_cells <- runif(1, min = 4000, max = 6000)
control_sim_2 <- simulate_count_scDesign2(control_model, 
                                          n_cell_new = target_cells, 
                                          sim_method = 'copula',
                                          cell_type_prop = cell_type_proportion_control)
saveRDS(control_sim_2, "/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/control_sim_2.rds")


target_cells <- runif(1, min = 4000, max = 6000)
control_sim_3 <- simulate_count_scDesign2(control_model, 
                                          n_cell_new = target_cells, 
                                          sim_method = 'copula',
                                          cell_type_prop = cell_type_proportion_control)
saveRDS(control_sim_3, "/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/control_sim_3.rds")

# Stimulated 
target_cells <- runif(1, min = 4000, max = 6000) # Randomly generates target cell number
stimulated_sim_1 <- simulate_count_scDesign2(stimulated_model, 
                                          n_cell_new = target_cells, 
                                          sim_method = 'copula',
                                          cell_type_prop = cell_type_proportion_stimulated)
saveRDS(stimulated_sim_1, "/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/stimulated_sim_1.rds")

target_cells <- runif(1, min = 4000, max = 6000)
stimulated_sim_2 <- simulate_count_scDesign2(stimulated_model, 
                                          n_cell_new = target_cells, 
                                          sim_method = 'copula',
                                          cell_type_prop = cell_type_proportion_stimulated)
saveRDS(stimulated_sim_2, "/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/stiulated_sim_2.rds")


target_cells <- runif(1, min = 4000, max = 6000)
stimulated_sim_3 <- simulate_count_scDesign2(stimulated_model, 
                                          n_cell_new = target_cells, 
                                          sim_method = 'copula',
                                          cell_type_prop = cell_type_proportion_stimulated)
saveRDS(control_sim_3, "/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/control_sim_3.rds")



### End 


                                              
                                              
                                              
                                              











