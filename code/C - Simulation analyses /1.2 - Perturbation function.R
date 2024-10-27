# SCRIPT CONTEXT 
#
# Using a traditional concept of damaged cells, develop a function 
# that, when applied to input data, will create perturbations to resemble damaged cells.  
# The idea is that the only truth we have is to add a known amount of damage to cases, 
# perform HVG and DEG identification compared to the case without perturbations,
# then see how well tool filter the perturbed data to mitigate the downstream 
# confounding effects of the perturbation. 
#
# This script creates the following cases: 
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
annotations <- genes(Hsap_edb, return.type = "data.frame")   

# Extract mitochondrial gene subsets
mt_genes <- annotations %>%
  dplyr::filter(grepl("MT-", gene_name) | grepl("mitochondrial", description, ignore.case = TRUE)) %>%
  pull(gene_name)

# Remove empty entries
mt_genes <- mt_genes[mt_genes != ""]
mito_genes <- mt_genes

# Introduce consistently-expressed and nuclear localized genes 
mito_genes <- append(c("MALAT1", "NEAT1", "XIST", "TSIX", "TUG1", "MEG3", "RN7SL1", "RPPH1"), mito_genes)


#-------------------------------------------------------------------------------
# Damage Perturbation Function 
#-------------------------------------------------------------------------------

# Function to perturb input count matrix to mirror a conventional damaged state by cellular burst. 
# To do so, we: 
# 
# 1. Reduce the expression of all cytoplasmic genes already expressed in the cell, mimicking partial loss of RNA transcripts (decreases UMI count) 
# 2. Remove the expression of a random subset cytoplasmic genes (make expression 0), mimicking complete loss of RNA transcripts (decreases feature count)


perturb_single_cell <- function(count_matrix, mito_genes, zero_out_constant = 0.0) {
  
  # Convert count matrix to matrix for manipulation
  count_matrix <- as.matrix(count_matrix)
  
  # Identify mitochondrial and non-mitochondrial gene indices
  mito_idx <- rownames(count_matrix) %in% mito_genes
  non_mito_idx <- !mito_idx
  
  # Precompute mitochondrial and non-mitochondrial gene sums for each cell
  mito_counts <- colSums(count_matrix[mito_idx, , drop=FALSE])
  non_mito_counts <- colSums(count_matrix[non_mito_idx, , drop=FALSE])
  total_counts <- colSums(count_matrix)
  
  # Randomly assign target mitochondrial percentages for each cell
  target_mito_pct <- runif(ncol(count_matrix), min = 20, max = 98)
  
  for (i in seq_len(ncol(count_matrix))) {
    # Step 1: Calculate target non-mitochondrial gene count to zero
    current_non_mito_genes <- sum(count_matrix[non_mito_idx, i] > 0)
    
    # Define a scaling curve to control zero-out proportion
    max_zero_out_proportion <- 0.32  # testing to see if higher helps smooth it out # 32
    zero_out_proportion <- max_zero_out_proportion * (1 - exp(-0.02 * target_mito_pct[i]))  # Subtle scaling with cap
    
    # Add a small constant if specified
    zero_out_proportion <- min(zero_out_proportion + zero_out_constant, max_zero_out_proportion)
    
    # Calculate the number of genes to zero out, ensuring it doesn't exceed available genes
    genes_to_zero <- min(floor(current_non_mito_genes * zero_out_proportion), current_non_mito_genes)
    
    # Randomly zero out the calculated number of non-mito genes
    if (genes_to_zero > 0) {
      non_mito_gene_indices <- which(non_mito_idx & count_matrix[, i] > 0)
      genes_to_zero_indices <- sample(non_mito_gene_indices, size = genes_to_zero)
      count_matrix[genes_to_zero_indices, i] <- 0
    }
    
    # Step 2: Calculate new mito and non-mito counts after zeroing
    mito_sum <- sum(count_matrix[mito_idx, i])
    non_mito_sum <- sum(count_matrix[non_mito_idx, i])
    total_counts[i] <- mito_sum + non_mito_sum  # Update total counts
    
    # Step 3: Apply scaling to remaining non-mito genes to reach target mito percentage
    target_mito_counts <- total_counts[i] * (target_mito_pct[i] / 100)
    reduction_factor <- (total_counts[i] - target_mito_counts) / max(non_mito_sum, 1)  # No restriction on reduction factor
    
    # Apply the calculated reduction factor directly
    count_matrix[non_mito_idx, i] <- count_matrix[non_mito_idx, i] * reduction_factor
    
    # Optional: display iteration details
    cat("Processed cell", i, "- Target mito %", round(target_mito_pct[i], 2), 
        "- Remaining non-mito genes:", sum(count_matrix[non_mito_idx, i] > 0), "\n")
  }
  
  return(count_matrix)
}

# Function to perturb an entire count matrix
introducePerturb <- function(seurat, 
                             mito_genes, 
                             percent_damage, 
                             project_name, 
                             output_dir = "/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2") {
  
  # Extract percentage of cells at random, these will be altered and then reintroduced 
  damaged_cell_number <- round(((percent_damage / 100) * (dim(seurat)[2])), 0)
  damaged_cell_selections <- sample(colnames(seurat), damaged_cell_number)
  to_be_perturbed <- subset(seurat, cells = damaged_cell_selections)
  to_be_perturbed_matrix <- to_be_perturbed@assays$RNA$counts
  
  # Remove these cells from the seurat object
  all_cells <- colnames(seurat)
  remaining_cells <- setdiff(all_cells, damaged_cell_selections)
  unperturbed <- subset(seurat, cells = remaining_cells)
  unperturbed_matrix <- unperturbed@assays$RNA$counts
  
  # Perturb the selected cells
  perturbed_matrix <- apply(to_be_perturbed_matrix, 2, perturb_single_cell, mito_genes = mito_genes)
  rownames(perturbed_matrix) <- rownames(to_be_perturbed_matrix)
  
  # Reorder the rows of perturbed_matrix to match the order of unperturbed_matrix
  perturbed_matrix <- perturbed_matrix[match(rownames(unperturbed_matrix), rownames(perturbed_matrix)), ]
  
  # Combine perturbed and unperturbed cells
  combined_matrix <- cbind(unperturbed_matrix, perturbed_matrix)
  combined_matrix <- as(combined_matrix, "dgCMatrix")
  combined_seurat <- CreateSeuratObject(counts = combined_matrix, assay = "RNA", project = paste0(seurat@project.name, "_perturbed"))
  
  # Edit the Seurat object to contain meaningful information 
  combined_seurat$orig.ident <- ifelse(colnames(combined_seurat) %in% damaged_cell_selections, "damaged", "cell")
  
  # Calculate standard quality control metrics ----- 
  
  # Isolate ribosomal genes (RPS and RPL)
  rb_genes <- annotations %>%
    dplyr::filter(grepl("^RPS|^RPL", gene_name)) %>%
    pull(gene_name)
  
  
  # Calculate the feature percentages and normalised expressions 
  combined_seurat$mt.percent <- PercentageFeatureSet(
    object   = combined_seurat,
    features = intersect(mito_genes, rownames(combined_seurat@assays$RNA)),
    assay    = "RNA"
  ) 
  
  combined_seurat$rb.percent <- PercentageFeatureSet(
    object   = combined_seurat,
    features = intersect(rb_genes, rownames(combined_seurat@assays$RNA)),
    assay    = "RNA"
  ) 
  
  combined_seurat$malat1.percent <- PercentageFeatureSet(
    object   = combined_seurat,
    features = "MALAT1",
    assay    = "RNA"
  ) 
  
  # Calculate psuedo-nuclear fraction score for DropletQC (identical to previous in section B, testing strategies 1.1)
  combined_seurat$nf_malat1 <- FetchData(combined_seurat, vars = "MALAT1", layer = "counts")
  combined_seurat$nf_malat1 <- (combined_seurat$nf_malat1 - min(combined_seurat$nf_malat1)) / (max(combined_seurat$nf_malat1) - min(combined_seurat$nf_malat1))
  
  # Transfer cell annotations 
  combined_seurat$celltype <- seurat$celltype[match(rownames(combined_seurat@meta.data), rownames(seurat@meta.data))]
  combined_seurat$celltype <- ifelse(combined_seurat$orig.ident == "damaged", "damaged", combined_seurat$celltype)
  
  # Save outputs 
  saveRDS(combined_seurat, file.path(output_dir, paste0(project_name, ".rds")))
  write.csv(combined_matrix, file.path(output_dir, paste0(project_name, "_matrix.csv")))
  
  
  # Visualise the output
  test <- NormalizeData(combined_seurat) %>%
    FindVariableFeatures() %>%
    ScaleData() %>% 
    RunPCA() %>%
    FindNeighbors(dims = 1:30) %>% 
    FindClusters() %>% 
    RunUMAP(dims = 1:30)
  
  # Note that these markers applicable for immune/PBMCs
  colours <- c("T" = "#A6BEAE",
               "Monocyte" = "#88A1CD",
               "DC" = "#A799C9",
               "B" = "#BDC5EE",
               "NK" = "#9DBD73",
               "damaged" = "red"
               )
  
  plot <- DimPlot(test, group.by = "celltype")  + 
    NoAxes() + ggtitle(project_name) + scale_color_manual(values = colours) + 
    theme(panel.border = element_rect(colour = "black"))

  ggsave(filename = file.path(output_dir, paste0(project_name, "_perturbed.png")), 
         plot, 
         width = 5, 
         height = 4, 
         units = "in",
         dpi = 300)
  
  
  return(combined_seurat)
  
}


#-------------------------------------------------------------------------------
# Generate perturbed data 
#-------------------------------------------------------------------------------

# 2.5 % damage added 
introducePerturb(seurat = control_sim_1, percent_damage = 2.5, mito_genes = mito_genes, project_name = "control_sim_1_2.5")
introducePerturb(seurat = control_sim_2, percent_damage = 2.5, mito_genes = mito_genes, project_name = "control_sim_2_2.5")
introducePerturb(seurat = control_sim_3, percent_damage = 2.5, mito_genes = mito_genes, project_name = "control_sim_3_2.5")
introducePerturb(seurat = stimulated_sim_1, percent_damage = 2.5, mito_genes = mito_genes, project_name = "stimulated_sim_1_2.5")
introducePerturb(seurat = stimulated_sim_2, percent_damage = 2.5, mito_genes = mito_genes, project_name = "stimulated_sim_2_2.5")
introducePerturb(seurat = stimulated_sim_3, percent_damage = 2.5, mito_genes = mito_genes, project_name = "stimulated_sim_3_2.5")


# 5 % damage added 
introducePerturb(seurat = control_sim_1, percent_damage = 5, mito_genes = mito_genes, project_name = "control_sim_1_5")
introducePerturb(seurat = control_sim_2, percent_damage = 5, mito_genes = mito_genes, project_name = "control_sim_2_5")
introducePerturb(seurat = control_sim_3, percent_damage = 5, mito_genes = mito_genes, project_name = "control_sim_3_5")
introducePerturb(seurat = stimulated_sim_1, percent_damage = 5, mito_genes = mito_genes, project_name = "stimulated_sim_1_5")
introducePerturb(seurat = stimulated_sim_2, percent_damage = 5, mito_genes = mito_genes, project_name = "stimulated_sim_2_5")
introducePerturb(seurat = stimulated_sim_3, percent_damage = 5, mito_genes = mito_genes, project_name = "stimulated_sim_3_5")

# 10 % damage added 
introducePerturb(seurat = control_sim_1, percent_damage = 10, mito_genes = mito_genes, project_name = "control_sim_1_10")
introducePerturb(seurat = control_sim_2, percent_damage = 10, mito_genes = mito_genes, project_name = "control_sim_2_10")
introducePerturb(seurat = control_sim_3, percent_damage = 10, mito_genes = mito_genes, project_name = "control_sim_3_10")
introducePerturb(seurat = stimulated_sim_1, percent_damage = 10, mito_genes = mito_genes, project_name = "stimulated_sim_1_10")
introducePerturb(seurat = stimulated_sim_2, percent_damage = 10, mito_genes = mito_genes, project_name = "stimulated_sim_2_10")
introducePerturb(seurat = stimulated_sim_3, percent_damage = 10, mito_genes = mito_genes, project_name = "stimulated_sim_3_10")

# 15 % damage added 
introducePerturb(seurat = control_sim_1, percent_damage = 15, mito_genes = mito_genes, project_name = "control_sim_1_15")
introducePerturb(seurat = control_sim_2, percent_damage = 15, mito_genes = mito_genes, project_name = "control_sim_2_15")
introducePerturb(seurat = control_sim_3, percent_damage = 15, mito_genes = mito_genes, project_name = "control_sim_3_15")
introducePerturb(seurat = stimulated_sim_1, percent_damage = 15, mito_genes = mito_genes, project_name = "stimulated_sim_1_15")
introducePerturb(seurat = stimulated_sim_2, percent_damage = 15, mito_genes = mito_genes, project_name = "stimulated_sim_2_15")
introducePerturb(seurat = stimulated_sim_3, percent_damage = 15, mito_genes = mito_genes, project_name = "stimulated_sim_3_15")

# 20 % damage added 
introducePerturb(seurat = control_sim_1, percent_damage = 20, mito_genes = mito_genes, project_name = "control_sim_1_20")
introducePerturb(seurat = control_sim_2, percent_damage = 20, mito_genes = mito_genes, project_name = "control_sim_2_20")
introducePerturb(seurat = control_sim_3, percent_damage = 20, mito_genes = mito_genes, project_name = "control_sim_3_20")
introducePerturb(seurat = stimulated_sim_1, percent_damage = 20, mito_genes = mito_genes, project_name = "stimulated_sim_1_20")
introducePerturb(seurat = stimulated_sim_2, percent_damage = 20, mito_genes = mito_genes, project_name = "stimulated_sim_2_20")
introducePerturb(seurat = stimulated_sim_3, percent_damage = 20, mito_genes = mito_genes, project_name = "stimulated_sim_3_20")







