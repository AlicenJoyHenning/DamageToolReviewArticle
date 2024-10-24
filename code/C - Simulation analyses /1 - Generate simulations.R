# SCRIPT CONTEXT 
#
# Simulate datasets with set number of cell types and known set of DEGs.
# Importantly, be able to add simulated damaged cells to the datasets: 
#   Case 1: 0 - 5 % damage
#   Case 2: 5 - 10 % damage
#   Case 3: 10 - 15 % damage
#   Case 4: 15 - 20 % damage
#   Case 5: 20 - 30 % damage
#
# Simulation tools: 
#   Splatter: https://bioconductor.org/packages/devel/bioc/vignettes/splatter/inst/doc/splatter.html#Session_information
#         - don't know about DEGs or how to add damage, also sce not seurat 
# 


#-------------------------------------------------------------------------------
# PREPARATIONS
#-------------------------------------------------------------------------------

library(devtools)
install_github("Vivianstats/scDesign")
library(scDesign)

#-------------------------------------------------------------------------------
# TRIALS 
#-------------------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("IRanges")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biobase")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("S4Vectors")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SummarizedExperiment")

library(devtools)
devtools::install_github("YosefLab/SymSim")
library(SymSim)

??SymSim

test_sim <- SimulateTrueCounts(ncells_total = 500, ngenes = 2000, randseed =  777)

meta_data <- test_sim$cell_meta
counts <- test_sim$counts
rownames(counts)

