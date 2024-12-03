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

packages <- c("AnnotationHub", "biomaRt", "dplyr", "ggplot2", "Seurat")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    library(pkg)
  }
}


# Load data -------

# Damaged cell samples in isolation 
apoptotic_isolated     <- readRDS("./D_Summarise_Results /data/raw_seurat/HEK293_apoptotic_processed.rds")
pro_apoptotic_isolated <- readRDS("./D_Summarise_Results /data/raw_seurat/HEK293_proapoptotic_processed.rds")
dead_SA928_isolated    <- readRDS("./D_Summarise_Results /data/raw_seurat/GM18507_dead.rds")
dying_SA928_isolated    <- readRDS("./D_Summarise_Results /data/raw_seurat/GM18507_dying.rds")
dead_SA604_isolated    <- readRDS("./D_Summarise_Results /data/raw_seurat/PDX_dead.rds")



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
edit_objects <- function(seurat, 
                         project_name, 
                         meta_data_edit = FALSE, 
                         reduce = TRUE, 
                         method = "UMAP"){
  
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
  
  if (method == "UMAP"){
  
  reduced <- NormalizeData(reduced) %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:30) %>%
    FindClusters(resolution = 1.2) %>% # High for finer differentiation 
    RunUMAP(dims = 1:30)  
  
  }
  
  if (method == "TSNE"){
    
    reduced <- NormalizeData(reduced) %>%
      FindVariableFeatures() %>%
      ScaleData() %>%
      RunPCA() %>%
      FindNeighbors(dims = 1:30) %>%
      FindClusters(resolution = 1.2) %>% # High for finer differentiation 
      RunTSNE(dims = 1:30)  
    
  }
  
  
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
dim(apoptotic_isolated)
apoptotic_reduced <- edit_objects(apoptotic_isolated, "apoptotic_isolated", method = "UMAP", reduce = FALSE)
apoptotic_reduced$damaged_population <- ifelse(apoptotic_reduced$seurat_clusters %in% c(1, 2, 5, 6), "B", "A")

# Pro-apoptotic
pro_apoptotic_reduced <- edit_objects(pro_apoptotic_isolated, "pro_apoptotic_isolated",  method = "UMAP", reduce = FALSE)
pro_apoptotic_reduced$damaged_population <- ifelse(pro_apoptotic_reduced$seurat_clusters %in% c(2, 5, 10), "B", "A")

# GM18607 dead
dead_SA928_reduced <- edit_objects(dead_SA928_isolated, "dead_SA928_isolated", method = "UMAP", reduce = FALSE, meta_data_edit = TRUE)
dead_SA928_reduced$damaged_population <- ifelse(dead_SA928_reduced$seurat_clusters %in% c(1, 4, 6, 8), "B", "A")

# GM18607 dying
dying_SA928_reduced <- edit_objects(dying_SA928_isolated, "dying_SA928_isolated", method = "UMAP", reduce = FALSE, meta_data_edit = TRUE)
dying_SA928_reduced$damaged_population <- ifelse(dying_SA928_reduced$seurat_clusters %in% c(2, 5), "A", "B")


# PDX dead 
dead_SA604_reduced  <- edit_objects(dead_SA604_isolated, "dead_SA604_isolated", method = "UMAP", reduce = FALSE, meta_data_edit = TRUE)
dead_SA604_reduced$damaged_population <- ifelse(dead_SA604_reduced $seurat_clusters %in% c(0, 2, 5), "B", "A")


#-------------------------------------------------------------------------------
# DAMAGED ISOLATION 
#-------------------------------------------------------------------------------

# Plotting damaged in isolation, annotating clusters for plotting  -------


# Function to specify plotting criteria
MyDimPlot <- function(seurat, 
                      group = "damaged_population", 
                      colours = damage_colours, 
                      project_name,
                      pt.size = 1,
                      output = "./D_Summarise_Results /img/isolated_damaged/" ){
  
  # Mitochondrial percentage
  mt_feature_plot <- FeaturePlot(seurat,
                          reduction = "umap",
                          features = "mt.percent",
                          pt.size = pt.size,
                          cols = c("#E6E6E6", "#001E5C")
                          ) +
    NoAxes() + 
    xlab("UMAP 1") + ylab("UMAP 2") +
    theme(plot.title = element_blank(),
          plot.subtitle = element_text(hjust = 0.5, vjust = 1),
          panel.border = element_rect(colour = "black", fill=NA, linewidth = 1),
          legend.position = "left")
  
  
  cluster_colours <- c("A" = "#E6E6E6",
                       "B" = "#001E5C")
  
  mt_violin <- VlnPlot(
    seurat, 
    pt.size = 0, 
    features = c("mt.percent"), 
    group.by = "damaged_population"
  ) + 
    scale_fill_manual(values = cluster_colours) +
    theme(
      axis.text = element_blank(), 
      axis.text.x = element_blank(),  
      axis.ticks = element_blank(),
      axis.title = element_blank(),     
      plot.title = element_blank(),         
      panel.border = element_rect(color = "black", fill = NA) 
    )
  
  
  
  # Ribosomal percentage
  rb_feature_plot <- FeaturePlot(seurat,
                              reduction = "umap",
                              features = "rb.percent",
                              pt.size = pt.size,
                              cols = c("#E7E7E7", "#B07D9A")
  ) +
    NoAxes() + 
    xlab("UMAP 1") + ylab("UMAP 2") +
    theme(plot.title = element_blank(),
          plot.subtitle = element_text(hjust = 0.5, vjust = 1),
          panel.border = element_rect(colour = "black", fill=NA, linewidth = 1),
          legend.position = "left")
  
  rb_violin <- VlnPlot(
    seurat, 
    pt.size = 0, 
    features = c("rb.percent"), 
    group.by = "damaged_population"
  ) + 
    scale_fill_manual(values = cluster_colours) +
    theme(
      axis.text = element_blank(), 
      axis.text.x = element_blank(),  
      axis.ticks = element_blank(),
      axis.title = element_blank(),     
      plot.title = element_blank(),         
      panel.border = element_rect(color = "black", fill = NA) 
    )
  
  
  # Library size 
  library_feature_plot <- FeaturePlot(seurat,
                                 reduction = "umap",
                                 features = "nFeature_RNA",
                                 pt.size = pt.size,
                                 cols = c("#E6E6E6", "#7F7F7F")
  ) +
    NoAxes() + 
    xlab("UMAP 1") + ylab("UMAP 2") +
    theme(plot.title = element_blank(),
          plot.subtitle = element_text(hjust = 0.5, vjust = 1),
          panel.border = element_rect(colour = "black", fill=NA, linewidth = 1),
          legend.position = "left")
  
  
  # Log transform library size
  seurat$library <- log10(seurat$nFeature_RNA)
  
  library_violin <- VlnPlot(
    seurat, 
    pt.size = 0, 
    features = c("library"), 
    group.by = "damaged_population"
  ) + 
    scale_fill_manual(values = cluster_colours) +
    theme(
      axis.text = element_blank(), 
      axis.text.x = element_blank(),  
      axis.ticks = element_blank(),
      axis.title = element_blank(),     
      plot.title = element_blank(),         
      panel.border = element_rect(color = "black", fill = NA) 
    )
  
  # MALAT1
  seurat$malat1.percent <- PercentageFeatureSet(
    object   = seurat,
    features = "MALAT1",
    assay    = "RNA"
  ) 
  
  
  # malat1_feature_plot <- FeaturePlot(seurat,
  #                                     reduction = "umap",
  #                                     features = "MALAT1",
  #                                     pt.size = pt.size,
  #                                     cols = c("#E6E6E6", "#B07D9A")
  # ) +
  #   NoAxes() + 
  #   xlab("UMAP 1") + ylab("UMAP 2") +
  #   theme(plot.title = element_blank(),
  #         plot.subtitle = element_text(hjust = 0.5, vjust = 1),
  #         panel.border = element_rect(colour = "black", fill=NA, linewidth = 1),
  #         legend.position = "left")
  
  
  
  malat1_violin <- VlnPlot(
    seurat, 
    pt.size = 0, 
    features = c("malat1.percent"), 
    group.by = "damaged_population"
  ) + 
    scale_fill_manual(values = cluster_colours) +
    theme(
      axis.text = element_blank(), 
      axis.text.x = element_blank(),  
      axis.ticks = element_blank(),
      axis.title = element_blank(),     
      plot.title = element_blank(),         
      panel.border = element_rect(color = "black", fill = NA) 
    )
  
  
  
  # Labels 
  # Plotting preparations 
  damage_colours <- c("A" = "#E6E6E6",
                     "B" = "#001E5C")
  
  cluster_plot <- DimPlot(seurat,
          reduction = "umap",
          group.by = group,
          pt.size = pt.size,
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
  
  combined <- mt_feature_plot | library_feature_plot  | rb_feature_plot | cluster_plot 
  
  # View and save 
  print(combined)
  ggsave(filename = paste0(output, project_name, ".png"),
         plot = combined,
         width = 10, height = 3.6, units = "in")
  
  
  return(list(mt_feature = mt_feature_plot, 
              mt_violin = mt_violin, 
              rb_feature = rb_feature_plot, 
              rb_violin = rb_violin,
             # malat1_feature = malat1_feature_plot, 
              malat1_violin = malat1_violin,
              library_feature = library_feature_plot, 
              library_violin = library_violin,
              cluster_plot = cluster_plot
              ))
}



apoptotic_plots <- MyDimPlot(apoptotic_reduced, project_name = "apoptotic")
pro_apoptotic_plots <- MyDimPlot(pro_apoptotic_reduced, project_name = "pro_apoptotic")
SA928_plots <- MyDimPlot(seurat = dead_SA928_reduced, project_name = "dead_SA928", pt.size = 1.7)
SA928_dying_plots <- MyDimPlot(dying_SA928_reduced, project_name = "dying_SA928", pt.size = 1.7)
SA604_plots <- MyDimPlot(dead_SA604_reduced, project_name = "dead_SA604", pt.size = 1.7)

# Dimplots 
mt_plots <- apoptotic_plots$mt_feature | (pro_apoptotic_plots$mt_feature + NoLegend()) | (SA928_plots$mt_feature + NoLegend()) | (SA928_dying_plots$mt_feature + NoLegend())| (SA604_plots$mt_feature + NoLegend())
rb_plots <- apoptotic_plots$rb_feature | (pro_apoptotic_plots$rb_feature + NoLegend()) | (SA928_plots$rb_feature + NoLegend()) | (SA928_dying_plots$rb_feature + NoLegend())| (SA604_plots$rb_feature + NoLegend())
library_plots <- apoptotic_plots$library_feature | (pro_apoptotic_plots$library_feature + NoLegend()) | (SA928_plots$library_feature + NoLegend()) | (SA928_dying_plots$library_feature + NoLegend()) | (SA604_plots$library_feature + NoLegend())
cluster_plots <- apoptotic_plots$cluster_plot | (pro_apoptotic_plots$cluster_plot + NoLegend()) | (SA928_plots$cluster_plot + NoLegend()) | (SA928_dying_plots$cluster_plot + NoLegend()) | (SA604_plots$cluster_plot + NoLegend())


ggsave(filename = "./D_Summarise_Results /img/isolated_damaged/mt_plots.png",
       plot = mt_plots,
       width = 14, height = 2.8, units = "in")

ggsave(filename = "./D_Summarise_Results /img/isolated_damaged/rb_plots.png",
       plot = rb_plots,
       width = 14, height = 2.8, units = "in")

ggsave(filename = "./D_Summarise_Results /img/isolated_damaged/library_plots.png",
       plot = library_plots,
       width = 14, height = 2.8, units = "in")

ggsave(filename = "./D_Summarise_Results /img/isolated_damaged/cluster_plots.png",
       plot = cluster_plots,
       width = 14, height = 2.8, units = "in")

# Violin plots 
mt_plots <- (apoptotic_plots$mt_violin + NoLegend()) | (pro_apoptotic_plots$mt_violin + NoLegend()) | (SA928_plots$mt_violin + NoLegend()) | (SA928_dying_plots$mt_violin + NoLegend())| (SA604_plots$mt_violin + NoLegend())
rb_plots <- (apoptotic_plots$rb_violin + NoLegend()) | (pro_apoptotic_plots$rb_violin + NoLegend()) | (SA928_plots$rb_violin + NoLegend()) | (SA928_dying_plots$rb_violin + NoLegend())| (SA604_plots$rb_violin + NoLegend())
malat1_plots <- (apoptotic_plots$malat1_violin + NoLegend()) | (pro_apoptotic_plots$malat1_violin + NoLegend()) | (SA928_plots$malat1_violin + NoLegend()) | (SA928_dying_plots$malat1_violin + NoLegend())| (SA604_plots$malat1_violin + NoLegend())
library_plots <- (apoptotic_plots$library_violin + NoLegend()) | (pro_apoptotic_plots$library_violin + NoLegend()) | (SA928_plots$library_violin + NoLegend()) | (SA928_dying_plots$library_violin + NoLegend()) | (SA604_plots$library_violin + NoLegend())
cluster_plots <- (apoptotic_plots$cluster_plot) | (pro_apoptotic_plots$cluster_plot + NoLegend()) | (SA928_plots$cluster_plot + NoLegend()) | (SA928_dying_plots$cluster_plot + NoLegend()) | (SA604_plots$cluster_plot + NoLegend())



ggsave(filename = "./D_Summarise_Results /img/isolated_damaged/mt_violins.png",
       plot = mt_plots,
       width = 14, height = 2.8, units = "in")

ggsave(filename = "./D_Summarise_Results /img/isolated_damaged/rb_violins.png",
       plot = rb_plots,
       width = 14, height = 2.8, units = "in")

ggsave(filename = "./D_Summarise_Results /img/isolated_damaged/library_violins.png",
       plot = library_plots,
       width = 14, height = 2.8, units = "in")

ggsave(filename = "./D_Summarise_Results /img/isolated_damaged/malat1_violins.png",
       plot = malat1_plots,
       width = 14, height = 2.8, units = "in")

ggsave(filename = "./D_Summarise_Results /img/isolated_damaged/cluster_plots.png",
       plot = cluster_plots,
       width = 14, height = 2.8, units = "in")


#-------------------------------------------------------------------------------
# TOOL DETECTION
#-------------------------------------------------------------------------------

# Proportion of damaged cells detected in each population  -------

# Load tool output 
apoptotic_tool_output <- read.csv("./C_Test_Strategies/data/benchmark_output/HEK293_apoptotic.csv")
pro_apoptotic_tool_output <- read.csv("./C_Test_Strategies/data/benchmark_output/HEK293_proapoptotic.csv")
GM18507_dead_tool_output <- read.csv("./C_Test_Strategies/data/benchmark_output/GM18507_dead.csv")
GM18507_dying_tool_output <- read.csv("./C_Test_Strategies/data/benchmark_output/GM18507_dying.csv")
PDX_tool_output <- read.csv("./C_Test_Strategies/data/benchmark_output/PDX_dead.csv")


# Function to count the proportion of each population detected by the tools ( population-specific detected cells / size of population )
find_damaged_population_detected <- function(seurat, 
                                             tool_output, 
                                             special = FALSE   # Zenodo files are 'special' and require additional formatting for compatibility 
                                             ) {
  
  if (special){
    
    # Extract barcodes of the damaged population of interest
    population_A <- subset(seurat, damaged_population == "A") 
    population_A$barcode <-  rownames(population_A@meta.data)
    population_A <- str_extract( population_A$barcode, "(?<=_)[A-Z]+(?=-)")
    
    population_B <- subset(seurat, damaged_population == "B")
    population_B$barcode <-  rownames(population_B@meta.data)
    population_B <- str_extract( population_B$barcode, "(?<=_)[A-Z]+(?=-)")
    
    # Add to the tool comparison output 
    tool_output$barcode <- str_extract(tool_output$X, "(?<=__)[A-Z]+(?=-)")
    tool_output$damaged_population <- ifelse(tool_output$barcode %in% population_A, "A", "-")
    tool_output$damaged_population <- ifelse(tool_output$barcode %in% population_B, "B", tool_output$damaged_population)
    
    
  } else {
    
    # Extract barcodes of the damaged population of interest
    population_A <- subset(seurat, damaged_population == "A") 
    population_A <- rownames(population_A@meta.data)
    population_B <- subset(seurat, damaged_population == "B") 
    population_B <- rownames(population_B@meta.data)
    
    # Add to the tool comparison output 
    tool_output$barcode <- str_extract(tool_output$X, "(?<=_)[A-Z]+(?=_)")
    tool_output$damaged_population <- ifelse( tool_output$barcode %in% population_A, "A", "-")
    tool_output$damaged_population <- ifelse( tool_output$barcode %in% population_B, "B", tool_output$damaged_population)
    
  }
    

  methods <- c("ddqc", "DropletQC", "ensembleKQC", "miQC", "SampleQC",  "scater",  "valiDrops", 
               "manual_fixed_mito", "manual_adaptive_mito", "manual_mito_ribo",  "manual_mito_ribo_library", "manual_library", "manual_malat1", "manual_malat1_mito_ribo")
  
  
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


apoptotic_tool_proportions <- find_damaged_population_detected(seurat = apoptotic_reduced, tool_output = apoptotic_tool_output)
pro_apoptotic_tool_proportions <- find_damaged_population_detected(pro_apoptotic_reduced, pro_apoptotic_tool_output)
dead_SA928_tool_proportions <- find_damaged_population_detected(seurat = dead_SA928_reduced, tool_output = GM18507_dead_tool_output, special = TRUE)
dying_SA928_tool_proportions <- find_damaged_population_detected(dying_SA928_reduced, GM18507_dying_tool_output, special = TRUE)
dead_SA604_tool_proportions <- find_damaged_population_detected(dead_SA604_reduced, PDX_tool_output, special = TRUE)


# Function to count the proportion of each population detected by the tools  ( population-specific detected cells / size of tool-specific detected cells )
find_damaged_popoulation_proportion <- function(seurat, 
                                                tool_output, 
                                                focus = FALSE,     # Whether to include other categories, or not and 'focus' on population A and B
                                                special = FALSE    # Zenodo files are 'special' and require additional formatting for compatibility 
                                                ) {
  
  if (special){
    
    # Extract barcodes of the damaged population of interest
    population_A <- subset(seurat, damaged_population == "A") 
    population_A$barcode <-  rownames(population_A@meta.data)
    population_A <- str_extract( population_A$barcode, "(?<=_)[A-Z]+(?=-)")
    
    population_B <- subset(seurat, damaged_population == "B")
    population_B$barcode <-  rownames(population_B@meta.data)
    population_B <- str_extract( population_B$barcode, "(?<=_)[A-Z]+(?=-)")
    
    # Add to the tool comparison output 
    tool_output$barcode <- str_extract(tool_output$X, "(?<=__)[A-Z]+(?=-)")
    tool_output$damaged_population <- ifelse(tool_output$barcode %in% population_A, "A", "-")
    tool_output$damaged_population <- ifelse(tool_output$barcode %in% population_B, "B", tool_output$damaged_population)
    
    
  } else {
    
    # Extract barcodes of the damaged population of interest
    population_A <- subset(seurat, damaged_population == "A") 
    population_A <- rownames(population_A@meta.data)
    population_B <- subset(seurat, damaged_population == "B") 
    population_B <- rownames(population_B@meta.data)
    
    # Add to the tool comparison output 
    tool_output$barcode <- str_extract(tool_output$X, "(?<=_)[A-Z]+(?=_)")
    tool_output$damaged_population <- ifelse( tool_output$barcode %in% population_A, "A", "-")
    tool_output$damaged_population <- ifelse( tool_output$barcode %in% population_B, "B", tool_output$damaged_population)
    
  }
  
  methods <- c("ddqc", "DropletQC", "ensembleKQC", "miQC", "SampleQC",  "scater",  "valiDrops", 
               "manual_fixed_mito", "manual_adaptive_mito", "manual_mito_ribo",  "manual_mito_ribo_library", "manual_library", "manual_malat1", "manual_malat1_mito_ribo")
  
  
  if (focus){
  
  detected_proportions <- data.frame(Method = character(),
                                     Population_A = numeric(), 
                                     Population_B = numeric(), 
                                     Population_other = numeric(), 
                                     stringsAsFactors = FALSE)
  
  
  for (method in methods) {
    
    tool_output$detected_A <- ifelse(tool_output[[method]] == "damaged" & tool_output$damaged_population == "A", "detected", "missed")
    tool_output$detected_B <- ifelse(tool_output[[method]] == "damaged" & tool_output$damaged_population == "B", "detected", "missed")
    tool_output$detected_other <- ifelse(tool_output[[method]] == "damaged" & (tool_output$damaged_population != "B" & tool_output$damaged_population != "A"), "detected", "missed")
    
    
    total_cell_detection <- sum(tool_output[[method]] == "damaged") 
      
    proportion_A <- (sum(tool_output$detected_A == "detected")/total_cell_detection) * 100
    proportion_B <- (sum(tool_output$detected_B == "detected")/total_cell_detection) * 100
    proportion_other <- (sum(tool_output$detected_other == "detected")/total_cell_detection) * 100
    
    # Create a new row 
    new_row <- data.frame(Method = method,
                          Population_A = proportion_A, 
                          Population_B = proportion_B,
                          Population_other = proportion_other)
    
    # Append the new row to results
    detected_proportions <- rbind(detected_proportions, new_row)
    
  }
  
  } else {
    
    detected_proportions <- data.frame(Method = character(),
                                       Population_A = numeric(), 
                                       Population_B = numeric(), 
                                       stringsAsFactors = FALSE)
    
    
    for (method in methods) {
      
      tool_output$detected_A <- ifelse(tool_output[[method]] == "damaged" & tool_output$damaged_population == "A", "detected", "missed")
      tool_output$detected_B <- ifelse(tool_output[[method]] == "damaged" & tool_output$damaged_population == "B", "detected", "missed")
      tool_output$detected_other <- ifelse(tool_output[[method]] == "damaged" & (tool_output$damaged_population != "B" & tool_output$damaged_population != "A"), "detected", "missed")
      
      
      total_cell_detection <- sum(tool_output$detected_A == "detected") + sum(tool_output$detected_B == "detected")
      
      proportion_A <- (sum(tool_output$detected_A == "detected")/total_cell_detection) * 100
      proportion_B <- (sum(tool_output$detected_B == "detected")/total_cell_detection) * 100
      
      # Create a new row 
      new_row <- data.frame(Method = method,
                            Population_A = proportion_A, 
                            Population_B = proportion_B)
      
      # Append the new row to results
      detected_proportions <- rbind(detected_proportions, new_row)
      
    }
    
  }
  
  
  return(detected_proportions)
  
}


apoptotic_tool_proportions_alt <- find_damaged_popoulation_proportion(seurat = apoptotic_reduced, tool_output = apoptotic_tool_output)
pro_apoptotic_tool_proportions_alt <- find_damaged_popoulation_proportion(pro_apoptotic_reduced, pro_apoptotic_tool_output)
dead_SA928_tool_proportions_alt <- find_damaged_popoulation_proportion(dead_SA928_reduced, GM18507_dead_tool_output, special = TRUE)
dying_SA928_tool_proportions_alt <- find_damaged_popoulation_proportion(dying_SA928_reduced, GM18507_dying_tool_output, special = TRUE)
dead_SA604_tool_proportions_alt <- find_damaged_popoulation_proportion(dead_SA604_reduced, PDX_tool_output, special = TRUE)

#-------------------------------------------------------------------------------
# PLOT
#-------------------------------------------------------------------------------

# Separate bar plots for each proportion ----

plot_propotions <- function(data, 
                            project_name,
                            output = "./D_Summarise_Results /img/isolated_damaged/" ){

  strategy_order <- c("ddqc", "DropletQC", "ensembleKQC", "miQC", "SampleQC",  "scater",  "valiDrops", 
                      "manual_fixed_mito", "manual_adaptive_mito", "manual_mito_ribo",  "manual_mito_ribo_library", "manual_library", "manual_malat1", "manual_malat1_mito_ribo")
  
  
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
    scale_fill_manual(values = c("Population_A" = "#E6E6E6", "Population_B" = "#001E5C")) +
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




# Stacked bar plot for tool-specific proportion -----


plot_proportions_alt <- function(data, 
                                 project_name,
                                 focus = FALSE,  # Whether to include other categories, or not and 'focus' on population A and B
                                 output = "./D_Summarise_Results /img/isolated_damaged/") {
  
  strategy_order <- c("ddqc", "DropletQC", "ensembleKQC", "miQC", "SampleQC",  "scater",  "valiDrops", 
                      "manual_fixed_mito", "manual_adaptive_mito", "manual_mito_ribo",  "manual_mito_ribo_library", "manual_library", "manual_malat1", "manual_malat1_mito_ribo")
  
  
  data$Method <- factor(data$Method, levels = strategy_order)
  data <- data[order(data$Method), ]
  
  if (focus){
    
  # Reshape the data to long format
  data_long <- data %>%
    pivot_longer(cols = c(Population_A, Population_B, Population_other), 
                 names_to = "Population", 
                 values_to = "Value") 
  
  # Create the stacked bar plot
  combined_plot <- ggplot(data_long, 
                          aes(x = Value, y = Method, fill = Population)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("Population_A" = "#E6E6E6", "Population_B" = "#001E5C", "Population_other" = "#808C98")) +
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
    )
  
  ggsave(filename = paste0(output, project_name, "_proportions_stacks.png"),
         plot = combined_plot,
         width = 6, height = 5, units = "in")
  
  } else {
    
    # Reshape the data to long format
    data_long <- data %>%
      pivot_longer(cols = c(Population_A, Population_B), 
                   names_to = "Population", 
                   values_to = "Value") 
    
    # Create the stacked bar plot
    combined_plot <- ggplot(data_long, 
                            aes(x = Value, y = Method, fill = Population)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("Population_A" = "#E6E6E6", "Population_B" = "#001E5C")) + 
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
      )
    
    ggsave(filename = paste0(output, project_name, "_proportions_stacks_focused.png"),
           plot = combined_plot,
           width = 7, height = 8, units = "in")
    
  }

  
}


plot_proportions_alt(apoptotic_tool_proportions_alt, "apoptotic")
plot_proportions_alt(pro_apoptotic_tool_proportions_alt, "pro_apoptotic")
plot_proportions_alt(dead_SA928_tool_proportions_alt, "dead_SA928")
plot_proportions_alt(dying_SA928_tool_proportions_alt, "dying_SA928")
plot_proportions_alt(dead_SA604_tool_proportions_alt, "dead_SA604")


### End

