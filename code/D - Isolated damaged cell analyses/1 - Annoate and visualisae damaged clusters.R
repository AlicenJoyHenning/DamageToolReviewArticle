# SCRIPT CONTEXT 
#
# Benchmarking damaged cell detection tools. Here, viewing damaged cells from 
# ground truth datasets in isolation & integrated with controls. Idea is to 
# create plot grid where you can see: 
#
# 1. The 2 clusters that form in damaged cells when viewed in isolation 
# 2. How these 2 clusters exist in the integrated objects 
# 3. Which damaged cells are detected by the damaged cell detection tools 
# 4. How the expression of certain QC measures differ between labelled populations 

#-------------------------------------------------------------------------------
# PREPARATIONS
#-------------------------------------------------------------------------------

# Load libraries -------

library(Seurat)
library(ggplot2)
library(svglite)

# Load data -------

# Damaged cell samples in isolation 
apoptotic_isolated     <- 
readRDS("/home/alicen/Projects/limiric/damage_left_behind_analysis//groundtruth/R_objects/apoptotic_isolated.rds")
pro_apoptotic_isolated <- 
readRDS("/home/alicen/Projects/limiric/damage_left_behind_analysis//groundtruth/R_objects/pro_apoptotic_isolated.rds")
dead_SA928_isolated    <- 
readRDS("/home/alicen/Projects/limiric/damage_left_behind_analysis/SA928_dead_isolated.rds")
dying_SA928_isolated   <- 
readRDS("/home/alicen/Projects/limiric/damage_left_behind_analysis/SA928_dying_isolated.rds")
dead_SA604_isolated    <- 
readRDS("/home/alicen/Projects/limiric/damage_left_behind_analysis/SA604_dead_isolated.rds")

# Integrated datasets 
apoptotic     <- 
readRDS("/home/alicen/Projects/limiric/damage_left_behind_analysis/groundtruth/R_objects/apoptotic.rds")
pro_apoptotic <- 
readRDS("/home/alicen/Projects/limiric/damage_left_behind_analysis/groundtruth/R_objects/pro_apoptotic.rds")
SA928_dead    <- 
readRDS("/home/alicen/Projects/limiric/damage_left_behind_analysis/groundtruth/R_objects/SA928_dead_live.rds")
SA928_dying   <- 
readRDS("/home/alicen/Projects/limiric/damage_left_behind_analysis/groundtruth/R_objects/SA928_dying_live.rds")


SA604_dead    <- 
readRDS("/home/alicen/Projects/limiric/damage_left_behind_analysis/groundtruth/R_objects/SA604_dead_live.rds")


#-------------------------------------------------------------------------------
# EDITING OBJECTS 
#-------------------------------------------------------------------------------

# Edit and save the objects with reductions for plotting -------

# Function to correct orig.ident if necessary & add dimensionality 
reductions for plotting 
edit_objects <- function(seurat, project_name, meta_data_edit = FALSE, 
reduce = TRUE){
  
  # Convert meta data if necessary 
  if (meta_data_edit){
    # Use labels. stored in cell barcodes   
    seurat$orig.ident <- rownames(seurat@meta.data)
    seurat$orig.ident <- sub("__.*", "", seurat$orig.ident)
  }
  
  
  # Dimensionality reduction with mtrb genes 
  
  # Get gene annotations for mitochondrial (MT) & ribosomal (RB) genes
  data("human_annotations", package = "limiric")
  annotations <- human_annotations
  
  mt_genes <- annotations %>%
    filter(grepl("MT-", gene_name)) %>% 
    pull(gene_name)
  
  # Isolate ribosomal genes (RPS and RPL)
  rb_genes <- annotations %>%
    filter(grepl("^RPS|^RPL", gene_name)) %>%
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
    FindClusters() %>%
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
  
  # Save objects 
  saveRDS(reduced, 
          
paste0("/home/alicen/Projects/limiric/damage_left_behind_analysis/groundtruth/R_objects/reduced_", 
project_name, ".rds"))

  return(reduced)
}

# Run the function 
apoptotic_isolated     <- edit_objects(apoptotic_isolated, 
"apoptotic_isolated")
pro_apoptotic_isolated <- edit_objects(pro_apoptotic_isolated, 
"pro_apoptotic_isolated")
dead_SA928_isolated    <- edit_objects(dead_SA928_isolated, 
"dead_SA928_isolated", meta_data_edit = TRUE)
dying_SA928_isolated   <- edit_objects(dying_SA928_isolated, 
"dying_SA928_isolated", meta_data_edit = TRUE)
dead_SA604_isolated    <- edit_objects(dead_SA604_isolated, 
"dead_SA604_isolated", meta_data_edit = TRUE)

apoptotic     <- edit_objects(apoptotic, "apoptotic")
pro_apoptotic <- edit_objects(pro_apoptotic, "pro_apoptotic")
# SA928_dead, SA928_dying, and SA604_dead already reduced


#-------------------------------------------------------------------------------
# DAMAGED ISOLATION 
#-------------------------------------------------------------------------------

# Plotting damaged in isolation, annotating clusters for plotting  -------


# Plotting preparations 
damage_colours <- c("A" = "#8FCBE2", "B" = "#BE9FC0", "C"  = "#C5C5C5", 
"D" = "#C5C5C5")
detected_damage_colours <- c("cell" = "#C5C5C5","NA" = "#C5C5C5","damaged" 
= "#C00054")


# Function to specify plotting criteria
MyDimPlot <- function(seurat, group = "label", colours = damage_colours){
  
  plot <- DimPlot(seurat,
          reduction = "umap",
          group.by = group,
          pt.size = 0.5,
         # label = TRUE,
          label.box = TRUE,
          label.size = 5,
          label.color = "white",
          repel = TRUE) +
    scale_color_manual(values =  colours) +
    scale_fill_manual(values =  colours) + 
    NoAxes() + NoLegend() +
    xlab("UMAP 1") + ylab("UMAP 2") +
    theme(plot.title = element_blank(),
          plot.subtitle = element_text(hjust = 0.5, vjust = 1),
          panel.border = element_rect(colour = "black", fill=NA, linewidth 
=1))
  
  return(plot)
  
}

# Function to make violin plot 
seurat <- apoptotic_isolated
MyViolin <- fuction(seurat){
  
  # Extract meta data 
  df <- seurat@meta.data

  # Create violin plot
  mt_violin <- ggplot(df, aes( x= label, y = mt.percent, fill = label)) + 
    geom_violin() + 
    theme_classic() +
    scale_fill_manual(values =  c("A" = "#8FCBE2", "B" = "#BE9FC0", "C"  = 
"#EFEFEF", "D" = "#A1A1A1")) + 
    theme(
      legend.position = "none",
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, linewidth =1)
    )
  
  # Create violin plot
  rb_violin <- ggplot(df, aes( x= label, y = rb.percent, fill = label)) + 
    geom_violin() + 
    theme_classic() +
    scale_fill_manual(values =  c("A" = "#8FCBE2", "B" = "#BE9FC0", "C"  = 
"#EFEFEF", "D" = "#EFEFEF")) + 
    theme(
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, linewidth =1)
    )  

  

  violin_plots <- mt_violin + rb_violin
  
}


# APOPTOTIC -----

# Annotate damaged clusters
FeaturePlot(apoptotic_isolated, reduction = "umap", features = 
"mt.percent", label = TRUE)      
apoptotic_isolated$label <- ifelse(apoptotic_isolated$seurat_clusters %in% 
c(0, 3),"B", "A")

# Rotate UMAP coordinates for matching the integrated UMAP
apoptotic_isolated@reductions$umap@cell.embeddings[, 1] <- 
-apoptotic_isolated@reductions$umap@cell.embeddings[, 1]  # Negate UMAP 1
apoptotic_isolated@reductions$umap@cell.embeddings[, 2] <- 
-apoptotic_isolated@reductions$umap@cell.embeddings[, 2]  # Negate UMAP 2
apoptotic_isolated_clusters <- MyDimPlot(apoptotic_isolated)

# Transfer cluster annotations from damaged to integrated object
apoptotic$label <- ifelse(
  sub("^[^_]+_", "", apoptotic$barcode) %in% 
rownames(apoptotic_isolated@meta.data), 
  apoptotic_isolated$label[match(sub("^[^_]+_", "", apoptotic$barcode), 
rownames(apoptotic_isolated@meta.data))], 
  "C"
)

# Translate UMAP coordinates about the y axis to match damaged plot
apoptotic@reductions$umap@cell.embeddings[, 1] <- 
-apoptotic@reductions$umap@cell.embeddings[, 1]  # Negate UMAP 1
FeaturePlot(apoptotic, reduction = "umap", features = "mt.percent", label 
= TRUE)      
apoptotic_clusters <- MyDimPlot(apoptotic)


# Create two control population categories (c, D)
DimPlot(apoptotic, label = TRUE)
apoptotic$label <- ifelse(apoptotic$seurat_clusters %in% c(), "C", 
apoptotic$label)
apoptotic$label <- ifelse((apoptotic$seurat_clusters %in% c(0, 1, 6, 7)) & 
(apoptotic$orig.ident == "HEK293_control"), "D", apoptotic$label)



# Compile labels of damaged cells detected across methods 
apoptotic$damage <- ifelse(
  rowSums(apoptotic@meta.data[, c("limiric", "miQC", "valiDrops", 
"DropletQC", "ddqc", 
                        "manual_all", "manual_mito_ribo", "manual_mito", 
                        "manual_malat1", "manual_mito_isolated")] == 
"damaged") >= 2, 
  "damaged",
  "cell"
)

apoptotic_tool_labels <- apoptotic@meta.data[, c("barcode", "damage")]
apoptotic_isolated$damage <- ifelse(rownames(apoptotic_isolated@meta.data) 
%in% apoptotic_tool_labels$barcode, apoptotic_tool_labels$damage, "cell")
apoptotic_tool_results <- MyDimPlot(apoptotic, group = "damage", 
detected_damage_colours)

# # Visualize on subset 
# apoptotic_A_subset <- subset(apoptotic, orig.ident == "apoptotic" & 
label == "A")
# apoptotic_B_subset <- subset(apoptotic, orig.ident == "apoptotic" & 
label == "B")
# 
# apoptotic_A_subset_tool_results <- MyDimPlot(apoptotic_A_subset, group = 
"damage", detected_damage_colours)
# apoptotic_B_subset_tool_results <- MyDimPlot(apoptotic_B_subset, group = 
"damage", detected_damage_colours)



# All plots together 
combined_plot <- 
  violin_plots |
  apoptotic_isolated_clusters | 
  apoptotic_clusters | 
  apoptotic_tool_results 
combined_plot


# PRO-APOPTOTIC -----

# Annotate damaged clusters
FeaturePlot(pro_apoptotic_isolated, reduction = "umap", features = 
"mt.percent", label = TRUE)      
pro_apoptotic_isolated$label <- 
ifelse(pro_apoptotic_isolated$seurat_clusters %in% c(0, 2, 6),"B", "A")

# Rotate UMAP coordinates for matching the integrated UMAP
# apoptotic_isolated@reductions$umap@cell.embeddings[, 1] <- 
-apoptotic_isolated@reductions$umap@cell.embeddings[, 1]  # Negate UMAP 1
# apoptotic_isolated@reductions$umap@cell.embeddings[, 2] <- 
-apoptotic_isolated@reductions$umap@cell.embeddings[, 2]  # Negate UMAP 2
pro_apoptotic_isolated_clusters <- MyDimPlot(pro_apoptotic_isolated)


# Transfer cluster annotations from damaged to integrated object
pro_apoptotic@meta.data$barcode <- sub("^[^_]+_[^_]+_", "", 
rownames(pro_apoptotic@meta.data))
pro_apoptotic$label <- ifelse(
  sub("^[^_]+_", "", pro_apoptotic$barcode) %in% 
rownames(pro_apoptotic_isolated@meta.data), 
  pro_apoptotic_isolated$label[match(sub("^[^_]+_", "", 
pro_apoptotic$barcode), rownames(pro_apoptotic_isolated@meta.data))], 
  "C"
)

# Plot the integrated objects with colours of cluster labels 
pro_apoptotic_clusters <- MyDimPlot(pro_apoptotic)

# Compile labels of damaged cells detected across methods 
pro_apoptotic$damage <- ifelse(
  rowSums(pro_apoptotic@meta.data[, c("limiric", "miQC", "valiDrops", 
"DropletQC", "ddqc", 
                                     "manual_all", "manual_mito_ribo", 
"manual_mito", 
                                    "manual_malat1", 
"manual_mito_isolated")] == "damaged") >= 3, 
  "damaged","cell"
)

pro_apoptotic_tool_labels <- pro_apoptotic@meta.data[, c("barcode", 
"damage")]
pro_apoptotic_isolated$damage <- 
ifelse(rownames(pro_apoptotic_isolated@meta.data) %in% 
pro_apoptotic_tool_labels$barcode, pro_apoptotic_tool_labels$damage, 
"cell")
pro_apoptotic_tool_results <- MyDimPlot(pro_apoptotic, group = "damage", 
detected_damage_colours)
pro_apoptotic_tool_results


# Isolated view of damaged cell populations 
pro_apoptotic_subset <- subset(pro_apoptotic, orig.ident == 
"pro_apoptotic")
pro_apoptotic_subset_tool_results <- MyDimPlot(pro_apoptotic_subset, group 
= "damage", detected_damage_colours)

# All plots together 
combined_plot <- pro_apoptotic_isolated_clusters | 
  pro_apoptotic_clusters | 
  pro_apoptotic_tool_results  
combined_plot



# SA604 -----

# Annotate damaged clusters
FeaturePlot(dead_SA604_isolated, reduction = "umap", features = 
"mt.percent", label = TRUE)      
dead_SA604_isolated$label <- ifelse(dead_SA604_isolated$seurat_clusters 
%in% c(0, 1),"B", "A")

# Rotate UMAP coordinates for matching the integrated UMAP
# apoptotic_isolated@reductions$umap@cell.embeddings[, 1] <- 
-apoptotic_isolated@reductions$umap@cell.embeddings[, 1]  # Negate UMAP 1
# apoptotic_isolated@reductions$umap@cell.embeddings[, 2] <- 
-apoptotic_isolated@reductions$umap@cell.embeddings[, 2]  # Negate UMAP 2
dead_SA604_isolated_clusters <- MyDimPlot(dead_SA604_isolated)


# Transfer cluster annotations from damaged to integrated object
SA604_dead@meta.data$barcode <- sub("^[^_]+_[^_]+_", "", 
rownames(SA604_dead@meta.data))
SA604_dead$label <- ifelse(
  sub("^[^_]+_", "", SA604_dead$barcode) %in% 
rownames(dead_SA604_isolated@meta.data), 
  dead_SA604_isolated$label[match(sub("^[^_]+_", "", SA604_dead$barcode), 
rownames(dead_SA604_isolated@meta.data))], 
  "C"
)

# Plot the integrated objects with colours of cluster labels 
SA604_dead_clusters <- MyDimPlot(SA604_dead)

# Compile labels of damaged cells detected across methods 
pro_apoptotic$damage <- ifelse(
  rowSums(pro_apoptotic@meta.data[, c("limiric", "miQC", "valiDrops", 
"DropletQC", "ddqc", 
                                      "manual_all", "manual_mito_ribo", 
"manual_mito", 
                                      "manual_malat1", 
"manual_mito_isolated")] == "damaged") >= 3, 
  "damaged","cell"
)

pro_apoptotic_tool_labels <- pro_apoptotic@meta.data[, c("barcode", 
"damage")]
pro_apoptotic_isolated$damage <- 
ifelse(rownames(pro_apoptotic_isolated@meta.data) %in% 
pro_apoptotic_tool_labels$barcode, pro_apoptotic_tool_labels$damage, 
"cell")
pro_apoptotic_tool_results <- MyDimPlot(pro_apoptotic, group = "damage", 
detected_damage_colours)
pro_apoptotic_tool_results


# Isolated view of damaged cell populations 
pro_apoptotic_subset <- subset(pro_apoptotic, orig.ident == 
"pro_apoptotic")
pro_apoptotic_subset_tool_results <- MyDimPlot(pro_apoptotic_subset, group 
= "damage", detected_damage_colours)

# All plots together 
combined_plot <- pro_apoptotic_isolated_clusters | 
  pro_apoptotic_clusters | 
  pro_apoptotic_tool_results  
combined_plot


# SA928 dead -----

# Annotate damaged clusters
FeaturePlot(dead_SA928_isolated , reduction = "umap", features = 
"mt.percent", label = TRUE)      
dead_SA928_isolated$label <- ifelse(dead_SA928_isolated$seurat_clusters 
%in% c(2, 3),"B", "A")

# Rotate UMAP coordinates for matching the integrated UMAP
# apoptotic_isolated@reductions$umap@cell.embeddings[, 1] <- 
-apoptotic_isolated@reductions$umap@cell.embeddings[, 1]  # Negate UMAP 1
# apoptotic_isolated@reductions$umap@cell.embeddings[, 2] <- 
-apoptotic_isolated@reductions$umap@cell.embeddings[, 2]  # Negate UMAP 2
dead_SA928_isolated_clusters <- MyDimPlot(dead_SA928_isolated)


# Transfer cluster annotations from damaged to integrated object
SA928_dead@meta.data$barcode <- sub("^[^_]+_[^_]+_", "", 
rownames(SA928_dead@meta.data))
SA928_dead$label <- ifelse(
  sub("^[^_]+_", "", SA928_dead$barcode) %in% 
rownames(dead_SA928_isolated@meta.data), 
  dead_SA928_isolated$label[match(sub("^[^_]+_", "", SA928_dead$barcode), 
rownames(dead_SA928_isolated@meta.data))], 
  "C"
)

# Plot the integrated objects with colours of cluster labels 
SA928_dead_clusters <- MyDimPlot(SA928_dead)

# Compile labels of damaged cells detected across methods 
pro_apoptotic$damage <- ifelse(
  rowSums(pro_apoptotic@meta.data[, c("limiric", "miQC", "valiDrops", 
"DropletQC", "ddqc", 
                                      "manual_all", "manual_mito_ribo", 
"manual_mito", 
                                      "manual_malat1", 
"manual_mito_isolated")] == "damaged") >= 3, 
  "damaged","cell"
)

pro_apoptotic_tool_labels <- pro_apoptotic@meta.data[, c("barcode", 
"damage")]
pro_apoptotic_isolated$damage <- 
ifelse(rownames(pro_apoptotic_isolated@meta.data) %in% 
pro_apoptotic_tool_labels$barcode, pro_apoptotic_tool_labels$damage, 
"cell")
pro_apoptotic_tool_results <- MyDimPlot(pro_apoptotic, group = "damage", 
detected_damage_colours)
pro_apoptotic_tool_results


# Isolated view of damaged cell populations 
pro_apoptotic_subset <- subset(pro_apoptotic, orig.ident == 
"pro_apoptotic")
pro_apoptotic_subset_tool_results <- MyDimPlot(pro_apoptotic_subset, group 
= "damage", detected_damage_colours)

# All plots together 
combined_plot <- pro_apoptotic_isolated_clusters | 
  pro_apoptotic_clusters | 
  pro_apoptotic_tool_results  
combined_plot
