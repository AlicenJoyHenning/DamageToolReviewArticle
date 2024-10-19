# SCRIPT CONTEXT 
#
# Helper function for the benchmarking pipeline to plot the output of benchmark function.
# Plots the output labels of each tool according to three different visualisations.

#-------------------------------------------------------------------------------
# OBTAIN LATEST ANNOTATIONS  
#-------------------------------------------------------------------------------

# Gene annotations for  QC metrics -----

# Connect to AnnotationHub (ah) and get the query search for the organisms of interest
ah <- AnnotationHub() 
Hsap_ahDb <- query(ah, pattern = c("Homo sapien", "EnsDb"), ignore.case = TRUE) 
Mmus_ahDb <- query(ah, pattern = c("Mus musculus", "EnsDb"), ignore.case = TRUE) 

# # Extract gene-level information from database of most up-to-date version
Hsap_versions <- mcols(Hsap_ahDb) 
Hsap_latest_version <- tail(rownames(Hsap_versions), n = 1)
Hsap_edb <- ah[[Hsap_latest_version]]
human_annotations <- genes(Hsap_edb, return.type = "data.frame")   

Mmus_versions <- mcols(Mmus_ahDb) 
Mmus_latest_version <- tail(rownames(Mmus_versions), n = 1)
Mmus_edb <- ah[[Mmus_latest_version]]
mouse_annotations <- genes(Mmus_edb, return.type = "data.frame")  


#-------------------------------------------------------------------------------
# HELPER FUNCTION DEFINED 
#-------------------------------------------------------------------------------

# Function for plotting using processed output from the benchmarking pipeline

BenchPlot <- function(seurat, 
                      output_file, 
                      organism = "Hsap", 
                      methods = c("ddqc", "DropletQC", "ensembleKQC", "miQC", "valiDrops", 
                                  "manual_all", "manual_mito_ribo", "manual_mito", "manual_malat1", "manual_mito_isolated")
                      
                      ) {
  
  # Prepare data frame for plotting 
  df <- seurat@meta.data
  df$mt.percent <- as.numeric(df$mt.percent)
  df$nFeature_RNA <- as.numeric(df$nFeature_RNA)
  df$nCount_RNA <- as.numeric(df$nCount_RNA)
  df$nf_malat1 <- as.numeric(df$nf_malat1)

  
  if ("nf" %in% colnames(df)){ df$nf <- as.numeric(df$nf) } else { df$nf <- as.numeric(log10(df$nf_malat1)) }
  

  # Define scatter plot themes
  plot_theme <- theme(
    plot.caption = element_text(hjust = 0.5, vjust = -0.5, size = 16),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(hjust = 0.5, vjust = -0.2, size = 16),
    axis.title.y = element_text(hjust = 0.5, size = 16),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "below",
    legend.box.background = element_rect(colour = "black"),
    strip.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.line = element_line(colour = "black")
  )
  
  
  # Function to create individual plots
  create_plot <- function(df, x, y, method, title) {
    ggplot(df, aes_string(x = x, y = y, color = method)) +
      geom_point(size = 0.5) +
      scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF")) +
      labs(title = title) + 
      xlab("") + ylab("") +
      plot_theme
  }
  
  # Create plots for nFeature_RNA vs mito.ratio
  plots1 <- lapply(methods, function(method) {
    create_plot(df, "nFeature_RNA", "mt.percent", method, method)
  })
  
  # Create plots for nf vs nCount_RNA
  plots2 <- lapply(methods, function(method) {
    create_plot(df, "nf_malat1", "nCount_RNA", method, method) 
  })
  
  # Create plots for nf vs nCount_RNA
  plots3 <- lapply(methods, function(method) {
    create_plot(df, "nf", "nCount_RNA", method, method) 
  })
  
  # Create plots for tSNE
  if (organism == "Hsap"){
    
    # Define human genes 
    annotations <- human_annotations 
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
    
  }
  if (organism == "Mmus"){
    
    # Define human genes 
    annotations <- mouse_annotations 
    malat1 <- c("Malat1")
    
    # Extract appropriate gene subsets
    mt_genes <- annotations %>%
      dplyr::filter(grepl("mt-", gene_name)) %>% 
      pull(gene_name)
    
    # Isolate ribosomal genes (RPS and RPL)
    rb_genes <- annotations %>%
      dplyr::filter(grepl("^Rsp|^Rpl", gene_name)) %>%
      pull(gene_name)
    
    # combine mt and rb genes
    mt_rb_genes <- unique(c(mt_genes, rb_genes))
    
  }
  
  # Reduce based on mt & rb genes only (this occurs in a separate Seurat object (limiric))
  seurat_view <- subset(seurat, features = intersect(mt_rb_genes, rownames(seurat@assays$RNA)))
  
  seurat_view <- NormalizeData(seurat_view, verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    FindNeighbors(dims = 1:10, verbose = FALSE) %>%
    FindClusters(resolution = 1, verbose = FALSE) %>%
    RunTSNE(dims = 1:10, verbose = FALSE, check_duplicates = FALSE)
  
  # Transfer labels 
  seurat_view$ddqc <- seurat$ddqc
  seurat_view$DropletQC <- seurat$DropletQC
  seurat_view$ensembleKQC <- seurat$ensembleKQC
  seurat_view$miQC <-  seurat$miQC
  seurat_view$valiDrops <- seurat$valiDrops
  seurat_view$manual_all <- seurat$manual_all
  seurat_view$manual_mito_ribo <- seurat$manual_mito_ribo
  seurat_view$manual_mito <- seurat$manual_mito
  seurat_view$manual_mito_isolated <- seurat$manual_mito_isolated
  seurat_view$manual_malat1 <- seurat$manual_malat1
  
  
  # Create tSNE plots
  plots4 <- lapply(methods, function(method) {
    DimPlot(seurat_view, reduction = 'tsne', group.by = method) + 
      scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF")) +
      NoAxes() + NoLegend() +
      theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
  })
  
  # Combine plots ------
  combined_plot1 <- plot_grid(plotlist = plots1, ncol = 10, rel_heights = c(1, 1))
  combined_plot2 <- plot_grid(plotlist = plots2, ncol = 10, rel_heights = c(1, 1))
  combined_plot3 <- plot_grid(plotlist = plots3, ncol = 10, rel_heights = c(1, 1))
  combined_plot4 <- plot_grid(plotlist = plots4, ncol = 10, rel_heights = c(1, 1))
  
  combined_plot <- plot_grid(
    combined_plot1, combined_plot2, combined_plot3, combined_plot4,
    ncol = 1, 
    rel_heights = c(1, 1, 1, 1)
  )
  
  combined_plot <- ggdraw(combined_plot) + 
    theme(
      plot.background = element_rect(fill = 'white', colour = NA)
    )
  
  # Save the plot
  ggsave(output_file, combined_plot, width = 36, height = 14, units = "in", dpi = 300, limitsize = FALSE)
  
  return(combined_plot)
  
}
