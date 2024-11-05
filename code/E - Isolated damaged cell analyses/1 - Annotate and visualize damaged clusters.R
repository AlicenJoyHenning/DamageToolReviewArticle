# SCRIPT CONTEXT 
#
# In view of the low precision and recall values, we decided to investigate the damaged cells from 
# the ground truth datasets in isolation. Ultimately, we wanted to understand trends in cells that are 
# not being detected. In exploration, we found that reducing the samples by their mitochondrial and 
# ribosomal genes resulted in two distinct populations forming across datasets. The following 
# script repeats this process as well as determines the proportion of each population detected as damaged 
# by each method. 
#
# 1. View the damaged populations 
#    - Reduce each ground truth case according to mitochondrial and ribosomal genes 
#    - Visualise reductions in plots
#    - Extract barcodes from each damaged population
# 
# 2. Examine tools recognition of each damaged population
#    Use damage labels from the ground truth tool benchmarks to calculate what 
#    propotion of each damaged population is identified by the tools


#-------------------------------------------------------------------------------
# PREPARATIONS
#-------------------------------------------------------------------------------

# Load libraries -------

library(Seurat)
library(ggplot2)
library(svglite)
library(AnnotationHub)


# Load data -------

# Damaged cell samples in isolation 
apoptotic_isolated     <- readRDS("/home/alicen/Projects/limiric/damage_left_behind_analysis/groundtruth/R_objects/apoptotic_isolated.rds")
pro_apoptotic_isolated <- readRDS("/home/alicen/Projects/limiric/damage_left_behind_analysis/groundtruth/R_objects/pro_apoptotic_isolated.rds")
dead_SA928_isolated    <- readRDS("/home/alicen/Projects/limiric/damage_left_behind_analysis/groundtruth/R_objects/reduced_dead_SA928_isolated.rds")
dying_SA928_isolated   <- readRDS("/home/alicen/Projects/limiric/damage_left_behind_analysis/groundtruth/R_objects/reduced_dying_SA928_isolated.rds")
dead_SA604_isolated    <- readRDS("/home/alicen/Projects/limiric/damage_left_behind_analysis/groundtruth/R_objects/reduced_dead_SA604_isolated.rds")


#-------------------------------------------------------------------------------
# OBTAIN LATEST ANNOTATIONS  
#-------------------------------------------------------------------------------

# Gene annotations for  QC metrics -----

# Extract gene-level information from database of most up-to-date version
ah <- AnnotationHub() 
Hsap_ahDb <- query(ah, pattern = c("Homo sapien", "EnsDb"), ignore.case = TRUE)  
Hsap_versions <- mcols(Hsap_ahDb) 
Hsap_latest_version <- tail(rownames(Hsap_versions), n = 1)
Hsap_edb <- ah[[Hsap_latest_version]]
annotations <- genes(Hsap_edb, return.type = "data.frame")   


#-------------------------------------------------------------------------------
# EDITING OBJECTS 
#-------------------------------------------------------------------------------

# Generate reductions -------

# Function to correct orig.ident if necessary & add dimensionality reductions for plotting 
edit_objects <- function(seurat, project_name, meta_data_edit = FALSE, reduce = TRUE){
  
  # Convert meta data if necessary 
  if (meta_data_edit){
    
    # Use labels. stored in cell barcodes   
    seurat$orig.ident <- rownames(seurat@meta.data)
    seurat$orig.ident <- sub("__.*", "", seurat$orig.ident)
  }
  
  # Dimensionality reduction with mtrb genes 
  
  # Get gene annotations for mitochondrial (MT) & ribosomal (RB) genes
  mt_genes <- annotations %>%
    dplyr::filter(grepl("MT-", gene_name)) %>% 
    #dplyr::filter(grepl("protein_coding", gene_biotype)) %>% 
    pull(gene_name)
  
  # Isolate ribosomal genes (RPS and RPL)
  rb_genes <- annotations %>%
    dplyr::filter(grepl("^RPS|^RPL", gene_name)) %>%
    #dplyr::filter(grepl("protein_coding", gene_biotype)) %>% 
    pull(gene_name)
  
  # combine mt and rb genes
  mt_rb_genes <- unique(c(mt_genes, rb_genes))
  
  if (reduce){
    
  # Reduce based on mt & rb genes only (this occurs in a separate Seurat object (limiric))
  reduced <- subset(seurat, 
                    features = intersect(mt_rb_genes, rownames(seurat@assays$RNA)))
  
  } else {reduced <- seurat}
  
  reduced <- NormalizeData(reduced) %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:30) %>%
    FindClusters(resolution = 1) %>%
    RunUMAP(dims = 1:30) %>% 
    RunTSNE(dims = 1:30)
  
  # Ensure gene expression is from original, not reduced, object
  reduced@assays$RNA <- seurat@assays$RNA
  
  # Calculate the percentage expression of mitochondrial genes 
  reduced$mt.percent <- PercentageFeatureSet(
    object   = reduced,
    features = intersect(mt_genes, rownames(reduced@assays$RNA)),
    assay    = "RNA"
  )
  
  # Define ribosomal expression
  reduced$rb.percent <- PercentageFeatureSet(
    object   = reduced,
    features = intersect(rb_genes, rownames(reduced@assays$RNA)),
    assay    = "RNA"
  )
  
  # Plot 
  plot <- FeaturePlot(reduced,
              features = "mt.percent",
              reduction = "umap", 
              label = TRUE, 
              label.size = 7) + 
    NoAxes() + NoLegend() +
    xlab("UMAP 1") + ylab("UMAP 2") + ggtitle(project_name) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, vjust = 1),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

  
  print(plot)
  
  return(reduced)
  
}

# Run the function & create labels for the two populations based on graph output
# A : high mito 
# B : low mito 

# Apoptotic
apoptotic_reduced <- edit_objects(apoptotic_isolated, "apoptotic_isolated")
apoptotic_reduced$damaged_population <- ifelse(apoptotic_reduced$seurat_clusters %in% c(1, 2), "B", "A")

# Pro-apoptotic
pro_apoptotic_reduced <- edit_objects(pro_apoptotic_isolated, "pro_apoptotic_isolated")
pro_apoptotic_reduced$damaged_population <- ifelse(pro_apoptotic_reduced$seurat_clusters %in% c(0, 2, 3), "B", "A")

# GM18607 dead
dead_SA928_reduced <- edit_objects(dead_SA928_isolated, "dead_SA928_isolated", meta_data_edit = TRUE)
dead_SA928_reduced$damaged_population <- ifelse(dead_SA928_reduced$seurat_clusters %in% c(2, 3, 5), "B", "A")

# GM18607 dying
dying_SA928_reduced <- edit_objects(dying_SA928_isolated, "dying_SA928_isolated", meta_data_edit = TRUE)
dying_SA928_reduced$damaged_population <- ifelse(dying_SA928_reduced$seurat_clusters %in% c(1), "A", "B")


# PDX dead 
dead_SA604_reduced  <- edit_objects(dead_SA604_isolated, "dead_SA604_isolated", meta_data_edit = TRUE)
dead_SA604_reduced$damaged_population <- ifelse(dead_SA604_reduced $seurat_clusters %in% c(0,1), "B", "A")


#-------------------------------------------------------------------------------
# DAMAGED ISOLATION 
#-------------------------------------------------------------------------------

# Plotting damaged in isolation, annotating clusters for plotting  -------

# Plotting preparations 
damage_colours <- c("A" = "lightgrey", 
                    "B" = "#8DC5BD") # "#CF9EBB") 

# Function to specify plotting criteria
MyDimPlot <- function(seurat, group = "damaged_population", colours = damage_colours){
  
  feature_plot <- FeaturePlot(seurat,
                          reduction = "umap",
                          features = "mt.percent",
                          pt.size = 1,
                          cols = c("#E6E6E6", "#001E5C")
                          ) +
    NoAxes() + 
    xlab("UMAP 1") + ylab("UMAP 2") +
    theme(plot.title = element_blank(),
          plot.subtitle = element_text(hjust = 0.5, vjust = 1),
          panel.border = element_rect(colour = "black", fill=NA, linewidth = 1),
          legend.position = "left")
  
  
  cluster_plot <- DimPlot(seurat,
          reduction = "umap",
          group.by = group,
          pt.size = 1,
         # label = TRUE,
          label.box = TRUE,
          label.size = 5,
          label.color = "black",
         # repel = TRUE
          ) +
    scale_color_manual(values =  colours) +
    scale_fill_manual(values =  colours) + 
    NoAxes() + NoLegend() +
    xlab("UMAP 1") + ylab("UMAP 2") +
    theme(plot.title = element_blank(),
          plot.subtitle = element_text(hjust = 0.5, vjust = 1),
          panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))
  
  combined <- feature_plot | cluster_plot
  
  # View and save 
  print(combined)
  return(combined)
  
}

MyDimPlot(apoptotic_reduced)
MyDimPlot(pro_apoptotic_reduced)
MyDimPlot(dead_SA928_reduced)
MyDimPlot(dying_SA928_reduced)
MyDimPlot(dead_SA604_reduced)

#-------------------------------------------------------------------------------
# TOOL DETECTION
#-------------------------------------------------------------------------------

# Proportion of damaged cells detected in each population  -------

# Load tool output 
apoptotic_tool_output <- read.csv("/home/alicen/Projects/ReviewArticle/benchmark_results/HEK293_apoptotic.csv")
pro_apoptotic_tool_output <- read.csv("/home/alicen/Projects/ReviewArticle/benchmark_results/HEK293_proapoptotic.csv")
GM18507_dead_tool_output <- read.csv("/home/alicen/Projects/ReviewArticle/benchmark_results/GM18507_dead.csv")
GM18507_dying_tool_output <- read.csv("/home/alicen/Projects/ReviewArticle/benchmark_results/GM18507_dying.csv")
PDX_tool_output <- read.csv("/home/alicen/Projects/ReviewArticle/benchmark_results/PDX_dead.csv")


# Function to count the proportion of each population detected by the tools 
find_damaged_population_detected <- function(seurat, tool_output, special = FALSE) {
  
  if (special){
    
    # Extract barcodes of the damaged population of interest
    population_A <- subset(seurat, damaged_population == "A") 
    population_A$barcode <-  rownames(population_A@meta.data)
    population_A <- as.character(sub(".*_(.*)", "\\1", population_A$barcode))
    
    population_B <- subset(seurat, damaged_population == "B")
    population_B$barcode <-  rownames(population_B@meta.data)
    population_B <- as.character(sub(".*_(.*)", "\\1", population_B$barcode))
    
  } else {
    
    # Extract barcodes of the damaged population of interest
    population_A <- subset(seurat, damaged_population == "A") 
    population_A <- rownames(population_A@meta.data)
    population_B <- subset(seurat, damaged_population == "B") 
    population_B <- rownames(population_B@meta.data)
    
  }
    
  # Add to the tool comparison output 
  tool_output$barcode <- sub(".*_(.*)", "\\1", tool_output$X)
    
  tool_output$damaged_population <- ifelse(tool_output$barcode %in% population_A, "A", "-")
  tool_output$damaged_population <- ifelse(tool_output$barcode %in% population_B, "B", tool_output$damaged_population)
  
  
  methods <- c("ddqc", "DropletQC", "ensembleKQC", "miQC", "SampleQC", "scater", "valiDrops",      
               "manual_all", "manual_mito_ribo","manual_mito","manual_malat1","manual_mito_isolated")
  
  detected_proportions <- data.frame(Method = character(),
                                     Population_A = numeric(), 
                                     Population_B = numeric(), 
                                     stringsAsFactors = FALSE)
                                   
  
  for (method in methods) {
  
    tool_output$detected_A <- ifelse(tool_output[[method]] == "damaged" & tool_output$damaged_population == "A", "detected", "missed")
    tool_output$detected_B <- ifelse(tool_output[[method]] == "damaged" & tool_output$damaged_population == "B", "detected", "missed")
    

    proportion_A <- (sum(tool_output$detected_A == "detected")/length(population_A)) * 100
    proportion_B <- (sum(tool_output$detected_B == "detected")/length(population_B)) * 100
    
    # Create a new row 
    new_row <- data.frame(Method = method,
                          Population_A = proportion_A, 
                          Population_B = proportion_B)
    
    # Append the new row to results
    detected_proportions <- rbind(detected_proportions, new_row)
    
}
  
  return(detected_proportions)

}


apoptotic_tool_proportions <- find_damaged_population_detected(apoptotic_reduced, apoptotic_tool_output)
pro_apoptotic_tool_proportions <- find_damaged_population_detected(pro_apoptotic_reduced, pro_apoptotic_tool_output)
dead_SA928_tool_proportions <- find_damaged_population_detected(dead_SA928_reduced, GM18507_dead_tool_output, special = TRUE)
dying_SA928_tool_proportions <- find_damaged_population_detected(dying_SA928_reduced, GM18507_dying_tool_output, special = TRUE)
dead_SA604_tool_proportions <- find_damaged_population_detected(dead_SA604_reduced, PDX_tool_output, special = TRUE)

#-------------------------------------------------------------------------------
# PLOT
#-------------------------------------------------------------------------------


plot_propotions <- function(data, 
                            project_name,
                            output = "/home/alicen/Projects/ReviewArticle/isolated_damaged/" ){

  strategy_order <- c(
    "ddqc", "DropletQC", "ensembleKQC", "miQC", "SampleQC", "scater", "valiDrops",
    "manual_all", "manual_mito_isolated", "manual_mito", "manual_mito_ribo", "manual_malat1"
  )
  
  data$Method <- factor(data$Method, levels = strategy_order)
  data <- data[order(data$Method), ]


  # Reshape the data to long format
  data_long <- data %>%
    pivot_longer(cols = c(Population_A, Population_B), 
                 names_to = "Population", 
                 values_to = "Value") 

  # Create the plot with facets
  combined_plot <- ggplot(data_long, 
                          aes(y = Method, x = Value, fill = Population)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("Population_A" = "#D3D3D3", "Population_B" = "#8DC5BD")) +
    labs(title = "",
         x = "",
         y = "") +
    theme_classic() +
    coord_cartesian(xlim = c(0, 100)) + 
    theme(legend.position = "none",
          axis.text.x = element_text(size = 16, vjust = -0.5),
          axis.ticks = element_blank(),
          axis.title.x = element_text(size = 16, vjust = -1.5),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 1.5),
          panel.spacing = unit(2, "lines"),
          strip.text = element_blank(),  
          strip.background = element_blank() 
    ) +
    facet_wrap(~ Population, scales = "free_x")

  ggsave(filename = paste0(output, project_name, "_proportions.png"),
         plot = combined_plot,
         width = 11, height = 5.5, units = "in")

}


plot_propotions(apoptotic_tool_proportions, "apoptotic")
plot_propotions(pro_apoptotic_tool_proportions, "pro_apoptotic")
plot_propotions(dead_SA928_tool_proportions, "dead_SA928")
plot_propotions(dying_SA928_tool_proportions, "dying_SA928")
plot_propotions(dead_SA604_tool_proportions, "dead_SA604")


