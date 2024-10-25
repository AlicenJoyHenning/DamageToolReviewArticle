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

# Install scDesign2
library(devtools)
devtools::install_github("JSB-UCLA/scDesign2")
install.packages("copula")


# Load all packages 
packages <- c("scDesign2", "copula", "plyr", "reshape2", "gridExtra", 
              "Seurat", "SeuratData", "ggpubr", "cowplot", "ggplot2")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    library(pkg)
  }
}


#-------------------------------------------------------------------------------
# Retrieve reference data
#-------------------------------------------------------------------------------

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


# scDesign2 requires matrix with celltype information as column labels
# data_mat <- readRDS(system.file("extdata", "mouse_sie_10x.rds", package = "scDesign2")) # ( EXAMPLE)

# Extract counts 
control_matrix <-   as.matrix(control@assays$RNA@counts)

# Extract celltypes 
control_celltypes <- as.data.frame(control@meta.data[, c("seurat_annotations")]) # name column annotations
control_celltypes$cells <- rownames(control@meta.data)
colnames(control_celltypes)[1] <- "celltype"

# Replace the colnames with the corresponding cell types
celltype_map <- setNames(control_celltypes$celltype, control_celltypes$cells)
colnames(control_matrix) <- celltype_map[colnames(control_matrix)]

# Verify 
table(colnames(control_matrix))







