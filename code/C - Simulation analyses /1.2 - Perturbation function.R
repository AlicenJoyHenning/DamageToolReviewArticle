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

# Read in simulated data ----
control_sim_1 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/simulated_references/control_sim_1_seurat.rds")
control_sim_2 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/simulated_references/control_sim_2_seurat.rds")
control_sim_3 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/simulated_references/control_sim_3_seurat.rds")
stimulated_sim_1 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/simulated_references/stimulated_sim_1_seurat.rds")
stimulated_sim_2 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/simulated_references/stimulated_sim_2_seurat.rds")
stimulated_sim_3 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/simulated_references/stimulated_sim_3_seurat.rds")


#seurat <- control_sim_1

#-------------------------------------------------------------------------------
# Preparations
#-------------------------------------------------------------------------------

# Gene annotations for mitochondrial and nuclear localised transcripts -----

# Connect to AnnotationHub (ah) and get the query search for the organisms of interest
ah <- AnnotationHub() 
Hsap_ahDb <- query(ah, pattern = c("Homo sapien", "EnsDb"), ignore.case = TRUE) 

# # Extract gene-level information from database of most up-to-date version
Hsap_versions <- mcols(Hsap_ahDb) 
Hsap_latest_version <- tail(rownames(Hsap_versions), n = 1)
Hsap_edb <- ah[[Hsap_latest_version]]
human_annotations <- genes(Hsap_edb, return.type = "data.frame")   

# Extract mitochondrial gene subsets
mt_genes <- human_annotations %>%
  dplyr::filter(grepl("MT-", gene_name) | grepl("mitochondrial", description, ignore.case = TRUE)) %>%
  pull(gene_name)

mt_genes <- append("MALAT1", mt_genes)

# Remove empty entries
mt_genes <- mt_genes[mt_genes != ""]
mito_genes <- mt_genes

# Nuclear localised genes 
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
nuclear_genes <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                       filters = 'go', 
                       values = c('GO:0048471', 'GO:001660'), 
                       mart = ensembl)




#-------------------------------------------------------------------------------
# Damage Perturbation Function 
#-------------------------------------------------------------------------------

# Function to perturb input count matrix to mirror a conventional damaged state by cellular burst. 
# To do so, we: 
# 
# 1. Reduce the expression of all cytoplasmic genes already expressed in the cell, mimicking partial loss of RNA transcripts (decreases UMI count) 
# 2. Remove the expression of a random subset cytoplasmic genes (make expression 0), mimicking complete loss of RNA transcripts (decreases feature count)

# Function to perturb a single cell's gene expression profile
introducePerturb_helper <- function(seurat, mito_genes, percent_damage) {
  
  # Extract cells ----
  
  # Extract percentage of cells at random, these will be altered and then reintroduced 
  damaged_cell_number <- round(((percent_damage/100) * (dim(seurat)[2])), 0)
  damaged_cell_selections <- sample(rownames(seurat@meta.data), damaged_cell_number)
  damaged <- subset(seurat, cells = damaged_cell_selections)
  damaged_matrix <- damaged@assays$RNA$counts
  
  # Remove these cells from the seurat object
  all_cells <- WhichCells(seurat)
  remaining_cells <- setdiff(all_cells, damaged_cell_selections)
  seurat <- subset(seurat, cells = remaining_cells)
  
  
  
  # Perturb cells -----
  
  # Randomly generate a target mitochondrial percentage between 60 and 100 for the cell
  target_mito_pct <- runif(1, min = 60, max = 100)
  
  # Separate mitochondrial and non-mitochondrial gene counts
  all_genes <- rownames(damaged_matrix)
  non_mito_genes <- all_genes[!(all_genes %in% mito_genes)]
  mito_counts <- damaged_matrix[mito_genes]
  non_mito_counts <- damaged_matrix[non_mito_genes]
  
  # Calculate current mitochondrial percentage
  current_mito_pct <- sum(mito_counts) / sum(damaged_matrix) * 100
  
  # If the current percentage is already higher than or equal to the target, return the cell unchanged
  if (current_mito_pct >= target_mito_pct) {
    return(cell_counts)
  }
  
  # Else, reduce non-mitochondrial gene expression to achieve the target mito percentage
  total_count <- sum(cell_counts)
  target_mito_sum <- total_count * target_mito_pct / 100
  mito_sum <- sum(mito_counts)
  reduction_factor <- (total_count - target_mito_sum) / sum(non_mito_counts)
  
  # Perturb the non-mitochondrial gene expression by randomly zeroing out some genes or reducing others
  perturbed_non_mito <- non_mito_counts * reduction_factor
  perturbed_non_mito <- ifelse(runif(length(perturbed_non_mito)) > 0.5, perturbed_non_mito, 0)
  
  # Recombine perturbed non-mitochondrial with unaltered mitochondrial counts
  perturbed_cell_counts <- cell_counts
  perturbed_cell_counts[!mito_genes] <- perturbed_non_mito
  
  return(perturbed_cell_counts)
}

# Function to perturb an entire count matrix
introducePerturb <- function(count_matrix, mito_genes) {
  
  # Convert dgCMatrix to a normal matrix 
  if (inherits(count_matrix, "dgCMatrix")) {
    count_matrix <- as.matrix(count_matrix)
  }
  
  # Apply the helper function to each column (cell) in the matrix
  perturbed_matrix <- apply(count_matrix, 2, introducePerturb_helper, mito_genes = mito_genes)
  
  # Convert the perturbed matrix back to dgCMatrix (better storage)
  perturbed_matrix <- as(perturbed_matrix, "dgCMatrix")
  
  return(perturbed_matrix)
  
}


#-------------------------------------------------------------------------------
# Retrieve datasets 
#-------------------------------------------------------------------------------

# SeuratData publicly available data ----
available_data <- AvailableData()
available_data <- subset(available_data, species == "human")

# IFNB-Stimulated & Control PBMCs
InstallData("ifnb")
data("ifnb")
ifnb <- UpdateSeuratObject(ifnb)

# Isolate control and treated cases from the original data 
control <- subset(ifnb, stim == "CTRL")    # 6548 cells 
stimulated <- subset(ifnb, stim == "STIM") # 7451 cells 


# Create technical replicates by random sub-sampling the data 
set.seed(111)
control_a_cells <- sample(rownames(control@meta.data), 5000)
control_a <- subset(control, cells = control_a_cells)
stimulated_a_cells <- sample(rownames(stimulated@meta.data), 5000)
stimulated_a <- subset(stimulated, cells = stimulated_a_cells) 

set.seed(112)
control_b_cells <- sample(rownames(control@meta.data), 5000)
control_b <- subset(control, cells = control_b_cells)
stimulated_b_cells <- sample(rownames(stimulated@meta.data), 5000)
stimulated_b <- subset(stimulated, cells = stimulated_b_cells) 

set.seed(113)
control_c_cells <- sample(rownames(control@meta.data), 5000)
control_c <- subset(control, cells = control_c_cells)
stimulated_c_cells <- sample(rownames(stimulated@meta.data), 5000)
stimulated_c <- subset(stimulated, cells = stimulated_c_cells) 



set.seed(123)  
# Example usage
# Create a logical vector indicating which genes in the count matrix are mitochondrial

mito_genes <- rownames(hPBMC_counts) %in% mt_genes
perturbed_count_matrix <- introducePerturb_matrix(count_matrix = hPBMC_counts, mito_genes = mito_genes)

test <- CreateSeuratObject(counts = perturbed_count_matrix, assay = "RNA", project = "test")

# Calculate the feature percentages and normalised expressions 
test$mt.percent <- PercentageFeatureSet(
  object   = test,
  features = intersect(mt_genes, rownames(test@assays$RNA)),
  assay    = "RNA"
) 

test$original.mt.percent <- hPBMC$mt.percent




