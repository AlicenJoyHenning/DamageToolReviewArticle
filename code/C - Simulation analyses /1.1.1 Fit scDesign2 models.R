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
# This script contains the code to create the model for generating simulated data using scDesign 2, 
# a scRNA-seq simulator that preserves gene names (required for all strategies) : 
# 1. Obtain annotated reference dataset 
#   -  SeuratData: IFNB containing treated (IFN-STIMULATED) and control (IFN-CONTROL) cells
# 2. Prepare the dataset for modelling 
# 3. Run the modelling scDesign2 function and save the models 
# 
# NOTE: For scDesign2 fit_model_scDesign2() function, it is advised to run in the background or in the terminal as it takes time (20 ~ 45 minutes in our experience)
# $ Rscript ./'1.1.1 Fit scDesign2 models.R' 


#-------------------------------------------------------------------------------
# PREPARATIONS
#-------------------------------------------------------------------------------

## Install scDesign2
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
set.seed(7)

# SeuratData publicly available data ----
available_data <- AvailableData()
available_data <- subset(available_data, species == "human")

# IFNB-Stimulated & Control PBMCs
InstallData("ifnb")
data("ifnb")
ifnb <- UpdateSeuratObject(ifnb)

# Recollect celltypes  (don't need fine annotations, only coarse)
ifnb$celltypes <- "na"
ifnb$celltypes <- ifelse(ifnb$seurat_annotations %in% c("CD14 Mono",  "CD16 Mono"), "Monocyte", ifnb$celltypes)
ifnb$celltypes <- ifelse(ifnb$seurat_annotations %in% c("CD4 Memory T", "T activated" , "CD4 Naive T", "CD8 T"), "T", ifnb$celltypes)
ifnb$celltypes <- ifelse(ifnb$seurat_annotations %in% c("B Activated" , "B"), "B", ifnb$celltypes)
ifnb$celltypes <- ifelse(ifnb$seurat_annotations %in% c("pDC","DC"), "DC", ifnb$celltypes)
ifnb$celltypes <- ifelse(ifnb$seurat_annotations == "NK", "NK", ifnb$celltypes)


# Isolate stimulated and control cells (used to simulate datasets in isolation- mimic reality)
control <- subset(ifnb, stim == "CTRL")    # 6548 cells 
stimulated <- subset(ifnb, stim == "STIM") # 7451 cells 


saveRDS(control, "/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/control_reference.rds")
saveRDS(stimulated, "/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/stimulated_reference.rds")

# Randomly down sample data for manageable model estimating times 
control_selections <- sample(rownames(control@meta.data), 2000)
control_subset <- subset(control, cells = control_selections)
stimulated_selections <- sample(rownames(stimulated@meta.data), 2000)
stimulated_subset <- subset(stimulated, cells = stimulated_selections)


# Extract counts 
control_matrix <-   as.matrix(control_subset@assays$RNA@counts)
stimulated_matrix <-   as.matrix(stimulated_subset@assays$RNA@counts)

# Extract celltypes 
celltypes <- as.data.frame(ifnb@meta.data[, c("celltypes")]) # name column annotations
celltypes$cells <- rownames(ifnb@meta.data)
colnames(celltypes)[1] <- "celltype"


# Replace the colnames of the matrices with the corresponding cell types
celltype_map <- setNames(celltypes$celltype, celltypes$cells)
colnames(control_matrix) <- celltype_map[colnames(control_matrix)]
colnames(stimulated_matrix) <- celltype_map[colnames(stimulated_matrix)]

saveRDS(control_matrix, "/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/control_reference_matrix.rds")
saveRDS(stimulated_matrix, "/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/stimulated_reference_matrix.rds")


#-------------------------------------------------------------------------------
# Run scDesign2
#-------------------------------------------------------------------------------

# scDesign2 -----

# Defining parameters for model fitting 

# Cell types to include 
cell_type_selection <- c("Monocyte", "DC", "T", "B", "NK")  


# Creating model with parameters
message("Control modelling...")

control_model <- fit_model_scDesign2(data = control_matrix, 
                                     cell_type_sel = cell_type_selection,
                                     sim_method = 'copula', 
                                     ncores = length(cell_type_selection)
)

message("Saving control model...")
saveRDS(control_model, "/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/control_model.rds")

# Stimulated model 
message("Stimulated modelling...")

stimulated_model <- fit_model_scDesign2(data = stimulated_matrix, 
                                     cell_type_sel = cell_type_selection,
                                     sim_method = 'copula', 
                                     ncores = length(cell_type_selection)
)

message("Saving stimulated model...")
saveRDS(stimulated_model, "/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/stimulated_model.rds")


### End
