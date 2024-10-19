# SCRIPT CONTEXT 
#
# Ground truth cases are different to non ground truth cases. Each ground truth sample 
# includes two separately sequenced samples, one for treated & sorted dead cells,
# and the other for untreated, sorted live cells. 
#
# For damaged cell detection, the count matrices of the damaged and control 
# samples are merged and integrated to resemble 'one psuedosample'. This simply
# concatenates the count matrices without losing their original labels. The 
# integration serves only for downstream clustering and visualisation, the 
# original (processed) count matrices for each sample are used for tool testing. 
#
# Further groundtruth processing done: 
# 1. Converting Ensembl gene symbols (ENSG0000....) to standard HGNC symbols (CD79A) where needed
# 2. Merging and integrating damaged and control samples 
#
# Note: After running this script, the user can proceed to tool testing using 
# the processed matrices. First, separately in python. Then, the remaining in R.



#-------------------------------------------------------------------------------
# PREPARATIONS
#-------------------------------------------------------------------------------

# Load libraries -------

packages <- c("cowplot", "devtools", "dplyr",
              "ensembldb", "ggplot2", "glmGamPoi", "Matrix", "png", 
              "presto", "scuttle", "Seurat", "tidyr")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    library(pkg)
  }
}


# Load data -------


HEK293_control <- readRDS("/home/alicen/Projects/ReviewArticle/R_objects/preprocess_output/HEK293_control_processed.rds")
HEK293_apoptotic <- readRDS("/home/alicen/Projects/ReviewArticle/R_objects/preprocess_output/HEK293_apoptotic_processed.rds") 
HEK293_proapoptotic <- readRDS("/home/alicen/Projects/ReviewArticle/R_objects/preprocess_output/HEK293_proapoptotic_processed.rds") 
GM18507_control <- readRDS("/home/alicen/Projects/ReviewArticle/R_objects/preprocess_output/GM18507_control_processed.rds")
GM18507_dying <- readRDS("/home/alicen/Projects/ReviewArticle/R_objects/preprocess_output/GM18507_dying_processed.rds")
GM18507_dead <- readRDS("/home/alicen/Projects/ReviewArticle/R_objects/preprocess_output/GM18507_dead_processed.rds")
PDX_control <- readRDS("/home/alicen/Projects/ReviewArticle/R_objects/preprocess_output/PDX_control_processed.rds")
PDX_dead <- readRDS("/home/alicen/Projects/ReviewArticle/R_objects/preprocess_output/PDX_dead_processed.rds")


#-------------------------------------------------------------------------------
# MERGER FUNCTION DEFINED 
#-------------------------------------------------------------------------------

# Merging and scale nf_malat1 scores for ground truth samples ------

# Function to merge
mergeandscale <- function(input_list, 
                          sample_IDs,
                          project_name,
                          output_dir = "/home/alicen/Projects/ReviewArticle/R_objects/preprocess_groundtruth/"
){
  
  # 1. Merge the samples ----
  
  # Check if input_list and sample_IDs have the same length
  if (length(input_list) != length(sample_IDs)) {
    stop("The length of input_list and sample_IDs must be the same.")
  }
  
  # Initialize the merged object with the first element
  first_obj <- input_list[[1]]
  
  # Loop through the rest of the elements to create their own list 
  remaining_objs <- list()
  
  for (i in 2:length(input_list)) {
    remaining_objs[[(i-1)]] <- input_list[[i]]
  }
  
  # Convert the list to the required format
  remaining_objs <- do.call(c, remaining_objs)
  
  # Run the Seurat merge function
  seurat <- merge(
    x = first_obj,
    y = remaining_objs,
    add.cell.ids = sample_IDs,  
    project = project_name 
  )
  
  # re-join layers after integration
  seurat[["RNA"]] <- JoinLayers(seurat[["RNA"]])
  
  # 2. Scale the nf malat1 scores -----
  seurat$nf_malat1 <- (seurat$nf_malat1 - min(seurat$nf_malat1)) / 
    (max(seurat$nf_malat1) - min(seurat$nf_malat1))
  
  
  # 3. Save 
  
  # Extract count matrix
  counts <- seurat@assays$RNA$counts
  matrix <- CreateSeuratObject(counts = counts, assay = "RNA")
  matrix <- suppressWarnings(as.matrix(matrix@assays$RNA$counts))
  write.csv(matrix, file = paste0("/home/alicen/Projects/ReviewArticle/python/input/", project_name, "_matrix_data.csv"))
  message("\u2714  Input for python prepared")
  
  saveRDS(seurat, paste0(output_dir, project_name, ".rds"))
  return(seurat)
  
}


# Run function 
HEK293_apo_control <- mergeandscale(list(HEK293_control, HEK293_apoptotic), c("control", "apoptotic"), "HEK293_apo_control")
HEK293_pro_control <- mergeandscale(list(HEK293_control, HEK293_proapoptotic), c("control", "proapoptotic"), "HEK293_pro_control")
GM18507_dying_control <- mergeandscale(list(GM18507_control, GM18507_dying), c("control", "dying"), "GM18507_dying_control")
GM18507_dead_control <- mergeandscale(list(GM18507_control, GM18507_dead), c("control", "dead"), "GM18507_dead_control")
PDX_dead_control <- mergeandscale(list(PDX_control, PDX_dead), c("control", "dead"), "PDX_dead_control")


### End 
