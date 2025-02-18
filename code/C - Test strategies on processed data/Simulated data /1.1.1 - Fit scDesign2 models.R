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
#   -  scRNAseq: PBMC samples from individuals with extensive response to influenza vaccine (HIGH) and 
#      a moderate to non-existent response (LOW) 
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
              "Seurat", "SeuratData", "ggpubr", "cowplot", "ggplot2", "scRNAseq")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    library(pkg)
  }
}


#-------------------------------------------------------------------------------
# Retrieve & prepare reference data
#-------------------------------------------------------------------------------

# scRNAseq package prvoides publicly available data ----
set.seed(7) # random down sampling of cells 

# PBMC from studying looking at vaccine responsiveness in lupus patients where
# two groups exist, those who with high & low antibody responses. 

# Retrieve data from the scRNAseq package 
pbmc <- fetchDataset("kotliarov-pbmc-2020", "2024-04-18")
metadata <- colData(pbmc)
metadata$phenotype_sample <- paste0(metadata$adjmfc.time, "_", metadata$sample)
counts <- pbmc@assays@data$counts
rownames(counts) <- rownames(pbmc)
colnames(counts) <- colnames(pbmc)
counts  <- as(counts, "sparseMatrix")

# Divide into samples according to response (one individual selected per response)

# HIGH RESPONDERS 
high_sample <- subset(metadata, phenotype_sample == "d0 high_207_d0")
high_sample <- rownames(high_sample)
high_counts  <- counts[, high_sample]

# Seeing that they were not included, annotate the cell types 
pbmc_high_seurat <- CreateSeuratObject(high_counts, assay = "RNA")
pbmc_high_seurat <- NormalizeData(pbmc_high_seurat) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:30)

# clusters <- DimPlot(pbmc_high_seurat)
# markers <- DotPlot(pbmc_high_seurat, features = c("MS4A1", "CD3E", "NKG7", "CD14", "ITGAX", "CLEC4C"))
# clusters | markers # View output

pbmc_high_seurat$celltypes <- ifelse(pbmc_high_seurat$seurat_clusters %in% c(0, 1, 4, 5), "T", "-")
pbmc_high_seurat$celltypes <- ifelse(pbmc_high_seurat$seurat_clusters == 2, "NK", pbmc_high_seurat$celltypes)
pbmc_high_seurat$celltypes <- ifelse(pbmc_high_seurat$seurat_clusters == 3, "Monocyte", pbmc_high_seurat$celltypes)
pbmc_high_seurat$celltypes <- ifelse(pbmc_high_seurat$seurat_clusters == 6, "B", pbmc_high_seurat$celltypes)
pbmc_high_seurat$celltypes <- ifelse(pbmc_high_seurat$seurat_clusters %in% c(8, 9), "DC", pbmc_high_seurat$celltypes)
pbmc_high_seurat$celltypes <- ifelse(pbmc_high_seurat$seurat_clusters == 10, "pDC", pbmc_high_seurat$celltypes)

# Extract celltypes 
celltypes <- as.data.frame(pbmc_high_seurat@meta.data[, c("celltypes")]) # name column annotations
celltypes$cells <- rownames(pbmc_high_seurat@meta.data)
colnames(celltypes)[1] <- "celltype"

# Replace the colnames of the matrix with the corresponding cell types
celltype_map <- setNames(celltypes$celltype, celltypes$cells)
colnames(high_counts) <- celltype_map[colnames(high_counts)]


# LOW RESPONDERS 
low_sample <- subset(metadata, phenotype_sample == "d0 low_277_d0")
low_sample <- rownames(low_sample)
low_counts  <- counts[, low_sample]

# Again, seeing that they were not included, annotate the cell types 
pbmc_low_seurat <- CreateSeuratObject(low_counts, assay = "RNA")
pbmc_low_seurat <- NormalizeData(pbmc_low_seurat) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:30)

# clusters <- DimPlot(pbmc_low_seurat)
# markers <- DotPlot(pbmc_low_seurat, features = c("MS4A1", "CD3E", "NKG7", "CD14", "ITGAX", "CLEC4C"))
# clusters | markers # View output

pbmc_low_seurat$celltypes <- ifelse(pbmc_low_seurat$seurat_clusters %in% c(0, 2, 5), "T", "-")
pbmc_low_seurat$celltypes <- ifelse(pbmc_low_seurat$seurat_clusters == 8, "NK", pbmc_low_seurat$celltypes)
pbmc_low_seurat$celltypes <- ifelse(pbmc_low_seurat$seurat_clusters == 1, "Monocyte", pbmc_low_seurat$celltypes)
pbmc_low_seurat$celltypes <- ifelse(pbmc_low_seurat$seurat_clusters %in% c(4, 6), "B", pbmc_low_seurat$celltypes)
pbmc_low_seurat$celltypes <- ifelse(pbmc_low_seurat$seurat_clusters %in% c(7, 9), "DC", pbmc_low_seurat$celltypes)
# pbmc_low_seurat$celltypes <- ifelse(pbmc_low_seurat$seurat_clusters == 6, "pDC", pbmc_low_seurat$celltypes) # too few 

# Extract celltypes 
celltypes <- as.data.frame(pbmc_low_seurat@meta.data[, c("celltypes")]) # name column annotations
celltypes$cells <- rownames(pbmc_low_seurat@meta.data)
colnames(celltypes)[1] <- "celltype"

# Replace the colnames of the matrix with the corresponding cell types
celltype_map <- setNames(celltypes$celltype, celltypes$cells)
colnames(low_counts) <- celltype_map[colnames(low_counts)]


# Randomly down sample data for manageable model estimating times 
low_selections <- sample(colnames(low_counts), 2000)
low_responders_subset <- low_counts[, low_selections] # 32738 features 2000 cells 
high_selections <- sample(colnames(high_counts), 2000)
high_responders_subset <- high_counts[, high_selections] # 32738  features 2000 cells 


saveRDS(pbmc_low_seurat, "/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/PBMC_low_seurat.rds")
saveRDS(pbmc_high_seurat, "/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/PBMC_high_seurat.rds")
saveRDS(low_responders_subset, "/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/PBMC_low_reference_matrix.rds")
saveRDS(high_responders_subset, "/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/PBMC_high_reference_matrix.rds")


#-------------------------------------------------------------------------------
# Run scDesign2
#-------------------------------------------------------------------------------

# scDesign2 -----

# Defining parameters for model fitting 

# Cell types to include 
cell_type_selection <- c("Monocyte", "DC", "T", "B", "NK")  


# Creating model with parameters for the high responders 
message("High responder modelling...")

pbmc_high_model <- fit_model_scDesign2(data = high_responders_subset, 
                                       cell_type_sel = cell_type_selection,
                                       sim_method = 'copula', 
                                       ncores = length(cell_type_selection)
)

message("Saving model...")
saveRDS(pbmc_high_model, "/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/PBMC_high_model.rds")


# Creating model with parameters for the low responders
message("Low responder modelling...")

pbmc_low_model <- fit_model_scDesign2(data = low_responders_subset, 
                                      cell_type_sel = cell_type_selection,
                                      sim_method = 'copula', 
                                      ncores = length(cell_type_selection)
)

message("Saving model...")
saveRDS(pbmc_low_model, "/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/PBMC_low_model.rds")

### End
