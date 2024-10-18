# SCRIPT CONTEXT 
#
# To keep tool comparisions as similar as possible, the processed counts are generated 
# first before input into tool testing. This is NB for tools not run in R, ddqc, so 
# standard input count matrices can be generated as a csv file. 
# 
# Note: 
#1. After running this script, the user will need to run the ddqc tool in python.
#     (ii - Run ddqc.ipynb)


#-------------------------------------------------------------------------------
# PREPARATIONS
#-------------------------------------------------------------------------------

# Load libraries -------

packages <- c("cowplot", "devtools", "dplyr", "ggplot2", "glmGamPoi", "Matrix", "robustbase",
              "png", "Seurat", "SoupX", "DoubletFinder", "tidyr", "DoubletFinder", "biomaRt",
              "limiric", "miQC", "SingleCellExperiment",
              "scuttle", "presto", "valiDrops", "DropletQC")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    library(pkg)
  }
}


#-------------------------------------------------------------------------------
# FUNCTION DEFINED 
#-------------------------------------------------------------------------------

# Function for processing ------

 # Create output to run in python
  # ddqc_counts <- seurat@assays$RNA$counts
  # ddqc <- CreateSeuratObject(counts = ddqc_counts, assay = "RNA", 
min.cells = 1)
  # mtx <- suppressWarnings(as.matrix(ddqc@assays$RNA$counts))
  # write.csv(mtx, file = paste0(output_path, "/ddqc/input/", 
project_name, "_ddqc_matrix_data.csv"))

  
  # cat("\u2714  Input for ddqc prepared \n")
