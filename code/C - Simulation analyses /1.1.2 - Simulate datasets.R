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
# This script contains the code to use the scDesign2 model generated using reference data and applying it to create simulated data : 
# 1. Define scDesign2 inputs 
# 2. Run scResign2 to create simulated datasets for the control and stimulated reference data
#     - Includes generating the counts 
#     - Renaming output to create Seurat object with relevant meta data columns 
#     - Calculating QC metrics so the Seurat object can enter straight into damaged stratergy testing 
#     - Outputting images of each simulated dataset to see clustering and marker expression 
# 


#-------------------------------------------------------------------------------
# PREPARATIONS
#-------------------------------------------------------------------------------

# Install scDesign2
#library(devtools)
#devtools::install_github("JSB-UCLA/scDesign2")
#install.packages("copula")

# Load all packages 
packages <- c("AnnotationHub", "scDesign2", "copula", "dplyr", "reshape2", "gridExtra", 
              "Seurat", "ggpubr", "gtable", "cowplot", "ggplot2")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    library(pkg)
  }
}


#-------------------------------------------------------------------------------
# Retrieve reference data (again)
#-------------------------------------------------------------------------------

# IFNB SeuratData reference samples 
control <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/simulated_references/control_reference.rds")
stimulated <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/simulated_references/stimulated_reference.rds")
control_matrix <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/simulated_references/control_reference_matrix.rds")
stimulated_matrix <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/simulated_references/stimulated_reference_matrix.rds")

# Read in models (1.1.1 - Fit scDesign2 models.R)
control_model <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/simulated_references/control_model.rds")
stimulated_model <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/simulated_references/stimulated_model.rds")


# Human gene annotations 

# Connect to AnnotationHub (ah) and get the query search for the organisms of interest
ah <- AnnotationHub() 
Hsap_ahDb <- query(ah, pattern = c("Homo sapien", "EnsDb"), ignore.case = TRUE) 

# # Extract gene-level information from database of most up-to-date version
Hsap_versions <- mcols(Hsap_ahDb) 
Hsap_latest_version <- tail(rownames(Hsap_versions), n = 1)
Hsap_edb <- ah[[Hsap_latest_version]]
annotations <- genes(Hsap_edb, return.type = "data.frame") 



#-------------------------------------------------------------------------------
# Run scDesign2
#-------------------------------------------------------------------------------

# Define inputs for scDesign2 ----

# Cell types to include 
cell_type_selection <- c("Monocyte", "DC", "T", "B", "NK")   # coarse PBMC annotations from reference dataset 

# Proportion which these celltype exist 
cell_type_proportion_control <- table(colnames(control_matrix))[cell_type_selection]
cell_type_proportion_stimulated <- table(colnames(stimulated_matrix))[cell_type_selection]


# Simulate datasets using scDesign2----

generate_simulations <- function(model_params,    # model fit created using scDesign2 fitting function 
                                 cell_type_prop,  # instead of mutinomial proportions, mimic proportions from the original data 
                                 project_name, 
                                 replicates = 3, 
                                 output_dir = "/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2") {
  
  # Helper function to generate and save a simulation dataset
  generate_and_save_simulation <- function(sim_number) {
    
    # Enforce (charisma) uniqueness (nerve and talent)
    set.seed(Sys.time() + sim_number)
    
    # Randomize cell number of simulation within reasonable range (4000 to 6000 cells)
    target_cells <- as.integer(runif(1, min = 4000, max = 6000))
    
    # scDesign2 to create count matrices 
    sim_matrix <- simulate_count_scDesign2(model_params = model_params, 
                                           n_cell_new = target_cells, 
                                           sim_method = 'copula',
                                           cell_type_prop = cell_type_prop)
    
    # Easier recording purposes 
    print(dim(sim_matrix))
    
    # Create matrix.csv and Seurat objects (rds) for downstream tool testing 
    # Store current colnames 
    celltypes <- colnames(sim_matrix)
    
    # Transfer gene names 
    if (project_name == "control"){
      
      # Gene names
      rownames(sim_matrix) <- rownames(control@assays$RNA@counts)
      
      # Extract cell barcodes from Seurat object
      cell_barcodes <- rownames(control@meta.data)
      
    } else {
      
      # Gene names 
      rownames(sim_matrix) <- rownames(stimulated@assays$RNA@counts)
      
      # Extract cell barcodes from Seurat object
      cell_barcodes <- rownames(stimulated@meta.data)
      
    } 
    
    # Transfer barcodes to matrix and save
    colnames(sim_matrix) <- cell_barcodes[1:ncol(sim_matrix)]
    write.csv(sim_matrix, file.path(output_dir, paste0(project_name, "_sim_", sim_number, "_matrix.csv")))
    
    # Create Seurat object
    seurat <- CreateSeuratObject(counts = sim_matrix, assay = "RNA", project = paste0(project_name, "_", sim_number))
    seurat$celltype <- celltypes
    
    # Calculate standard quality control metrics ----- 
    
    # Retrieve the corresponding annotations for the organism of interest 
    # Define human genes 
    malat1 <- c("MALAT1")
      
    # Extract appropriate gene subsets
    mt_genes <- annotations %>%
      dplyr::filter(grepl("MT-", gene_name)) %>% 
      pull(gene_name)
      
    # Isolate ribosomal genes (RPS and RPL)
    rb_genes <- annotations %>%
     dplyr::filter(grepl("^RPS|^RPL", gene_name)) %>%
     pull(gene_name)
      
    # combine mt and rb genes
    mt_rb_genes <- unique(c(mt_genes, rb_genes))
      
    
    # Calculate the feature percentages and normalised expressions 
    seurat$mt.percent <- PercentageFeatureSet(
      object   = seurat,
      features = intersect(mt_genes, rownames(seurat@assays$RNA)),
      assay    = "RNA"
    ) 
    
    seurat$rb.percent <- PercentageFeatureSet(
      object   = seurat,
      features = intersect(rb_genes, rownames(seurat@assays$RNA)),
      assay    = "RNA"
    ) 
    
    seurat$malat1.percent <- PercentageFeatureSet(
      object   = seurat,
      features = malat1,
      assay    = "RNA"
    ) 
    
    seurat$malat1 <- FetchData(seurat, vars = malat1)
    
    # Calculate psuedo-nuclear fraction score for DropletQC 
    seurat$nf_malat1 <- FetchData(seurat, vars = "MALAT1", layer = "counts")
    
    # min-max normalization: scales values so the min value maps to 0 and max maps to 1 (make the expression values more like nf scores)
    seurat$nf_malat1 <- (seurat$nf_malat1 - min(seurat$nf_malat1)) / (max(seurat$nf_malat1) - min(seurat$nf_malat1))
    
    
    # Save final object (resembles output of 1.1 Preprocess for section B)
    saveRDS(seurat, file.path(output_dir, paste0(project_name, "_sim_", sim_number, "_seurat.rds")))
    
    
    # Viewing the data ----
    # (my favourite part)
    test <- NormalizeData(seurat) %>%
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
    
    clusters <- DimPlot(test, group.by = "celltype")  + 
      NoAxes() + ggtitle(project_name) + scale_color_manual(values = colours) + 
      theme(panel.border = element_rect(colour = "black"))
    
    # MArkers 
    T_cells <- FeaturePlot(test, "CD3E") + NoAxes() + NoLegend() + theme(panel.border = element_rect(colour = "black"))
    B_cells <- FeaturePlot(test, "MS4A1") + NoAxes() + NoLegend() + theme(panel.border = element_rect(colour = "black"))
    NK_cells <- FeaturePlot(test, "NKG7") + NoAxes() + NoLegend() + theme(panel.border = element_rect(colour = "black"))
    Mon_cells <- FeaturePlot(test, "CD14") + NoAxes() + NoLegend() + theme(panel.border = element_rect(colour = "black"))
    
    plot <- clusters + ((T_cells | B_cells) / (NK_cells | Mon_cells))
    
    ggsave(filename = file.path(output_dir, paste0(project_name, "_sim_", sim_number, "_markers.png")), 
           plot, 
           width = 10, 
           height = 4, 
           units = "in",
           dpi = 300)
  }
  
  # Generate and save the simulations
  for (i in 1:replicates) {
    generate_and_save_simulation(i)
  }
}

# Generate simulated datasets§
generate_simulations(model_params = control_model, 
                     cell_type_prop = cell_type_proportion_control,
                     project_name = "control")

generate_simulations(model_params = stimulated_model, 
                     cell_type_prop = cell_type_proportion_stimulated,
                     project_name = "stimulated")


### End 


                                              
                                              
                                              
                                              











