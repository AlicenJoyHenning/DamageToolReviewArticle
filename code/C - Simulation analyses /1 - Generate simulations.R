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

# Getting standard PBMC dataset to build the simulated data
library(BiocManager)
BiocManager::install("scRNAseq")
library(scRNAseq) # only contains 3 tiny datasets 

# Looking for PBMC dataset ~ 6000 cells 
?scRNAseq

library(scater)
help(package = "scater")


# design_data() creates a new scRNA-seq matrix using an existing one 
realcount1 = readRDS(system.file("extdata", "astrocytes.rds", package = "scDesign"))
simcount1 = design_data(realcount = realcount1, 
                        S = 1e7, 
                        ncell = 1000, 
                        ngroup = 1, 
                        ncores = 1)


simcount1[1:3, 1:3]

realcount1

