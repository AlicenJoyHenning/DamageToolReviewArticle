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

# Sample sheet (guide)
sample_sheet <- read.csv("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/unperturbed_data/sample_sheet.csv")
sample_sheet$Origin_Condition <- paste0(sample_sheet$Origin, sample_sheet$Condition)
sample_sheet$Name <- paste0(sample_sheet$Origin_Condition, "_R", sample_sheet$Replicate, "_SD", sample_sheet$Sequencing_depth, "_CT", sample_sheet$Celltype_number, "_CN", sample_sheet$Cell_number)
sample_sheet$Plot <- ifelse()

# Core scRNAseq references
high_pbmc <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/filter_genes_0.2/models/PBMC_high_seurat.rds")
high_pbmc_matrix <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/filter_genes_0.2/models/PBMC_high_reference_matrix.rds")
low_pbmc <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/filter_genes_0.2/models/PBMC_low_seurat.rds")
low_pbmc_matrix <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/filter_genes_0.2/models/PBMC_low_reference_matrix.rds")
covid <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/filter_genes_0.2/models/Cellline_covid_seurat.rds")
covid_matrix <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/filter_genes_0.2/models/Cellline_covid_reference_matrix.rds")
healthy <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/filter_genes_0.2/models/Cellline_healthy_seurat.rds")
healthy_matrix <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/filter_genes_0.2/models/Cellline_healthy_reference_matrix.rds")


# Read in models (Fit scDesign2 models.R)
low_pbmc_model <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/filter_genes_0.2/models/PBMC_low_model.rds")
high_pbmc_model <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/filter_genes_0.2/models/PBMC_high_model.rds")
covid_cellline_model <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/filter_genes_0.2/models/cellline_covid_model.rds")
healthy_cellline_model <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/filter_genes_0.2/models/cellline_healthy_model.rds")


#-------------------------------------------------------------------------------
# Run scDesign2
#-------------------------------------------------------------------------------


# Simulate datasets using scDesign2----

simulate_matrices <- function(
    origin, # "PH" "PL" "CH" "CC"
    replicate, # 1 2 3 
    target_sequencing_depth, # 200, 500, 1000, 2000, 5000
    target_celltypes,  # 1, 2, 3, 6
    target_cellnumber, # 200, 500, 1000, 2000, 5000
){
  # Define inputs for scDesign
  
  # Dataset of origin related to counts & model
  if (origin == "PH"){
    count_matrix <- high_pbmc_matrix
    model <- high_pbmc_model
  }
  
  if (origin == "PL"){
    count_matrix <- low_pbmc_matrix
    model <- low_pbmc_model
  }
  
  if (origin == "CC"){
    count_matrix <- covid_matrix
    model <- covid_cellline_model
  }
  
  if (origin == "CH"){
    count_matrix <- healthy_matrix
    model <- healthy_cellline_model
  }
  
  
  # Cell type selection 
  Cell_Type_6 <- c("B", "DC", "Monocyte", "NK", "T", "pDC") 
  Cell_Type_3 <- c("B", "Monocyte", "T") 
  Cell_Type_2 <- c("CD8", "CD4") 
  Cell_Type_1 <- c("CD4") 
  
  
  
  # Proportion which these cell types exist 
  cell_type_proportion <- table(colnames(high_pbmc_matrix))[cell_type_selection]
  
  
}




    # scDesign2 to create count matrices 
    sim_matrix <- simulate_count_scDesign2(model_params = model_params, 
                                           n_cell_old = 3000, # true for all 
                                           n_cell_new = target_cells, 
                                           total_count_new = target_library_size, 
                                           sim_method = 'copula',
                                           cell_type_prop = cell_type_proportion)
  
    
    # Create matrix.csv and Seurat objects (rds) for downstream tool testing 
    # Store current colnames 
    celltypes <- colnames(sim_matrix)
    
    # Ensure clean cell barcodes
    colnames(count_matrix) <- names(colnames(count_matrix))
    
    # Transfer gene names & extract cell barcodes 
    rownames(sim_matrix) <- rownames(input_matrix)
    
    # Transfer barcodes to matrix and save
    colnames(sim_matrix) <- sample(cell_barcodes, ncol(sim_matrix))
    write.csv(sim_matrix, file.path(output_dir, paste0(project_name, "_sim_", sim_number, "_matrix.csv")))
    
    # Create Seurat object
    #new_matrix <- as.matrix(sim_matrix)
    #rownames(new_matrix) <- make.unique(rownames(sim_matrix))
    #colnames(new_matrix) <- make.unique(colnames(sim_matrix))
    seurat <- CreateSeuratObject(counts = sim_matrix, assay = "RNA", project = paste0(project_name, "_", sim_number))
    seurat$celltype <- celltypes
    
    # Calculate standard quality control metrics ----- 
    
    # Retrieve the corresponding annotations for the organism of interest 
    # Define human genes 
    malat1 <- c("MALAT1")
    seurat$malat1 <- FetchData(seurat, vars = malat1)
    
    # Calculate the feature percentages and normalised expressions 
    seurat$mt.percent <- PercentageFeatureSet(
      object   = seurat,
      features = rownames(seurat@assays$RNA)[grepl("^MT-", rownames(seurat@assays$RNA))],
      assay    = "RNA"
    ) 
    
    seurat$rb.percent <- PercentageFeatureSet(
      object   = seurat,
      features = rownames(seurat@assays$RNA)[grepl("^RPS|^RPL", rownames(seurat@assays$RNA))],
      assay    = "RNA"
    ) 
    
    seurat$malat1.percent <- PercentageFeatureSet(
      object   = seurat,
      features = malat1,
      assay    = "RNA"
    ) 

    # Calculate psuedo-nuclear fraction score for DropletQC 
    seurat$nf_malat1 <- FetchData(seurat, vars = "MALAT1", layer = "counts")
    
    # min-max normalization: scales values so the min value maps to 0 and max maps to 1 (make the expression values more like nf scores)
    seurat$nf_malat1 <- (seurat$nf_malat1 - min(seurat$nf_malat1)) / (max(seurat$nf_malat1) - min(seurat$nf_malat1))
    
    # Save final object (resembles output of 1.1 Preprocess for section B)
    saveRDS(seurat, file.path(output_dir, paste0(project_name, "_sim_", sim_number, "_seurat.rds")))
    
    
    # Viewing the data ----
    
    if (generate_plot) {
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
    
  }
  
  # Generate and save the simulations
  for (i in 1:replicates) {
    generate_and_save_simulation(i)
  }
  
}

# Generate simulated datasets
generate_simulations(model_params = low_pbmc_model, 
                     cell_type_prop = cell_type_proportion_low,
                     project_name = "low")

generate_simulations(model_params = high_pbmc_model, 
                     cell_type_prop = cell_type_proportion_high,
                     project_name = "high")


simulated <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/high_sim_1_seurat.rds")

library_size <- colSums(simulated)

### End 




                                              
                                              
                                              
                                              











