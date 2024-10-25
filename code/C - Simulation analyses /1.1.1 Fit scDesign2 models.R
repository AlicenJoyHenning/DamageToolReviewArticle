# SCRIPT CONTEXT 
#
# Running the model fitting for the control and stimulated data 


#-------------------------------------------------------------------------------
# PREPARATIONS
#-------------------------------------------------------------------------------

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
set.seed(7)

# SeuratData publicly available data ----
available_data <- AvailableData()
available_data <- subset(available_data, species == "human")

# IFNB-Stimulated & Control PBMCs
InstallData("ifnb")
data("ifnb")
ifnb <- UpdateSeuratObject(ifnb)

# Recollect celltypes 
ifnb$celltypes <- "na"
ifnb$celltypes <- ifelse(ifnb$seurat_annotations %in% c("CD14 Mono",  "CD16 Mono"), "Monocyte", ifnb$celltypes)
ifnb$celltypes <- ifelse(ifnb$seurat_annotations %in% c("CD4 Memory T", "T activated" , "CD4 Naive T", "CD8 T"), "T", ifnb$celltypes)
ifnb$celltypes <- ifelse(ifnb$seurat_annotations %in% c("B Activated" , "B"), "B", ifnb$celltypes)
ifnb$celltypes <- ifelse(ifnb$seurat_annotations %in% c("pDC","DC"), "DC", ifnb$celltypes)
ifnb$celltypes <- ifelse(ifnb$seurat_annotations == "NK", "NK", ifnb$celltypes)


# Isolate stimulated and control cells (used to simulate datasets in isolation)
control <- subset(ifnb, stim == "CTRL")    # 6548 cells 
stimulated <- subset(ifnb, stim == "STIM") # 7451 cells 

# Randomly downsample data for manageable model estimating times 
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
