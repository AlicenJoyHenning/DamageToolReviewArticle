# SCRIPT CONTEXT 
#
# Identify a conserved signature of damaged cells to develop a function 
# that when applied to input data will perturb them to resemble damaged cells.  
# The idea is that the only truth we have is to add a known amount of damage to cases, 
# perform HVG and DEG identification compared to the case without perturbations,
# then see how well tool filter the perturbed data to mitigate the downstream 
# confounding effects of the perturbation. 
#
# This script explores the true damaged populations and thier control counterparts, 
# isolates model 
#   Case 1: 0 - 5 % damage
#   Case 2: 5 - 10 % damage
#   Case 3: 10 - 15 % damage
#   Case 4: 15 - 20 % damage
#   Case 5: 20 - 30 % damage


#-------------------------------------------------------------------------------
# PREPARATIONS
#-------------------------------------------------------------------------------

packages <- c("AnnotationHub","biomaRt", "cowplot", "devtools", "dplyr", "ggplot2", "glmGamPoi", "Matrix",
              "png", "Seurat", "SeuratData", "miQC", "SingleCellExperiment", "loupeR")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    library(pkg)
  }
}


#-------------------------------------------------------------------------------
# Load datasets  
#-------------------------------------------------------------------------------

