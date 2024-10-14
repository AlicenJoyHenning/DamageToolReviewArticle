# Plot the labels of each tool according to the three different visualisaitons 


BenchPlot <- function(seurat, output_file, organism) {
  
  # Prepare data
  df <- seurat@meta.data
  df$mito.ratio <- as.numeric(df$mito.ratio)
  df$nFeature_RNA <- as.numeric(df$nFeature_RNA)
  df$nf <- as.numeric(df$nf)
  df$nCount_RNA <- as.numeric(df$nCount_RNA)
  
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
  
  # List of methods of interest
  methods <- c("ddqc", "DropletQC", "limiric", "miQC", "valiDrops", "manual_all", "manual_mito_ribo", "manual_mito", "manual_malat1", "manual_mito_isolated")
  
  
  # Function to create individual plots
  create_plot <- function(df, x, y, method, title) {
    ggplot(df, aes_string(x = x, y = y, color = method)) +
      geom_point(size = 0.5) +
      scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF", "empty_droplet" = "darkgrey", "Uncertain" = "grey")) +
      labs(title = title) + 
      xlab("") + ylab("") +
      plot_theme
  }
  
  # Create plots for nFeature_RNA vs mito.ratio
  plots1 <- lapply(methods, function(method) {
    create_plot(df, "nFeature_RNA", "mito.ratio", method, method)
  })
  
  # Create plots for nf vs nCount_RNA
  plots2 <- lapply(methods, function(method) {
    create_plot(df, "nf", "nCount_RNA", method, method) + scale_y_log10()
  })
  
  # Create plots for tSNE
  if (organism == "Hsap") { data("mt_rb_genes", package = "limiric")  }
  
  if (organism == "Mmus") {
    
    # Get mouse annotations
    data("mouse_annotations", package = "limiric")
    
    annotations <- mouse_annotations
    
    # Get gene annotations for mitochondrial genes
    mt_genes <- annotations[grep("mt-", annotations$gene_name, perl = TRUE), ]
    mt_genes <- mt_genes[grepl("mitochondrially encoded", mt_genes$description, perl = TRUE), ]
    mt_genes <- mt_genes %>% pull(gene_name)
    
    # isolate ribosomal genes
    rb_genes <- annotations[grepl("ribosomal", annotations$description, perl = TRUE), ]
    rb_genes <- rb_genes[grepl("protein_coding", rb_genes$gene_biotype, perl = TRUE), ]
    rb_genes <- rb_genes %>% pull(gene_name)
    
    # combine mt and rb genes
    mt_rb_genes <- c(mt_genes, rb_genes)
    mt_rb_genes <- unique(mt_rb_genes)
    
  }
  
  # Reduce based on mt & rb genes only (this occurs in a separate Seurat object (limiric))
  seurat_view <- subset(seurat, 
                        features = intersect(mt_rb_genes, rownames(seurat@assays$RNA)))
  
  seurat_view <- NormalizeData(seurat_view, verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    FindNeighbors(dims = 1:10, verbose = FALSE) %>%
    FindClusters(resolution = 1, verbose = FALSE) %>%
    RunTSNE(dims = 1:10, verbose = FALSE, check_duplicates = FALSE)
  
  # Transfer labels 
  seurat_view$limiric <- seurat$limiric
  seurat_view$DropletQC <- seurat$DropletQC
  seurat_view$Seurat <- seurat$threshold
  seurat_view$miQC <-  seurat$miQC
  seurat_view$valiDrops <- seurat$valiDrops
  seurat_view$ddqc <- seurat$ddqc
  
  # Create tSNE plots
  plots3 <- lapply(methods, function(method) {
    DimPlot(seurat_view, reduction = 'tsne', group.by = method) + 
      scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF", "empty_droplet" = "darkgrey", "Uncertain" = "grey")) +
      NoAxes() + NoLegend() +
      theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
  })
  
  # Combine plots ------
  combined_plot1 <- plot_grid(plotlist = plots1, ncol = 6, rel_heights = c(1, 1))
  combined_plot2 <- plot_grid(plotlist = plots2, ncol = 6, rel_heights = c(1, 1))
  combined_plot3 <- plot_grid(plotlist = plots3, ncol = 6, rel_heights = c(1, 1))
  
  combined_plot <- plot_grid(
    combined_plot1, combined_plot2, combined_plot3,
    ncol = 1, 
    rel_heights = c(1, 1, 1)
  )
  
  combined_plot <- ggdraw(combined_plot) + 
    theme(
      plot.background = element_rect(fill = 'white', colour = NA)
    )
  
  # Save the plot
  ggsave(output_file, combined_plot, width = 24, height = 10, units = "in", dpi = 300)
  
  return(combined_plot)
  
}
