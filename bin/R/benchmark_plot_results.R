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
  
  # Create individual plots
  Seurat <- ggplot(df, aes(x = nFeature_RNA, y = mito.ratio, color = threshold)) +
    geom_point(size = 0.5) +
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF")) +
    labs(title = "Threshold") + 
    xlab("") + ylab("") +
    plot_theme
  
  DropletQC <- ggplot(df, aes(x = nFeature_RNA, y = mito.ratio, color = DropletQC)) +
    geom_point(size = 0.5) +
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF", "empty_droplet" = "darkgrey")) +
    labs(title = "DropletQC") + 
    xlab("") + ylab("") +
    plot_theme
  
  limiric <- ggplot(df, aes(x = nFeature_RNA, y = mito.ratio, color = limiric)) +
    geom_point(size = 0.5) +
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF")) +
    labs(title = "limiric") + 
    xlab("") + ylab("") +
    plot_theme
  
  miQC <- ggplot(df, aes(x = nFeature_RNA, y = mito.ratio, color = miQC)) +
    geom_point(size = 0.5) +
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF")) +
    labs(title = "miqc") + 
    xlab("") + ylab("") +
    plot_theme
  
  valiDrops <- ggplot(df, aes(x = nFeature_RNA, y = mito.ratio, color = valiDrops)) +
    geom_point(size = 0.5) +
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF")) +
    #scale_color_manual(values = c("live" = "#D5D5D5", "dead" = "#9AA9FF")) +
    labs(title = "valiDrops") + 
    xlab("") + ylab("") +
    plot_theme
  
  ddqc <- ggplot(df, aes(x = nFeature_RNA, y = mito.ratio, color = ddqc)) +
    geom_point(size = 0.5) +
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF", "Uncertain" = "grey")) +
    # scale_color_manual(values = c("True" = "#D5D5D5", "False" = "#9AA9FF", "Uncertain" = "grey")) +
    labs(title = "ddqc") + 
    xlab("") + ylab("") +
    plot_theme
  
  # Additional plots
  Seurat2 <- ggplot(df, aes(x = nf, y = nCount_RNA, color = threshold)) +
    scale_y_log10() + geom_point(size = 0.5) +
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF")) +
    labs(title = "Seurat") + 
    xlab("") + ylab("") +
    plot_theme
  
  DropletQC2 <- ggplot(df, aes(x = nf, y = nCount_RNA, color = DropletQC)) +
    scale_y_log10() + geom_point(size = 0.5) +
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF", "empty_droplet" = "darkgrey")) +
    # scale_color_manual(values = c("cell" = "#D5D5D5", "damaged_cell" = "#9AA9FF", "empty_droplet" = "darkgrey")) +
    labs(title = "DropletQC") +
    xlab("") + ylab("") +
    plot_theme
  
  limiric2 <- ggplot(df, aes(x = nf, y = nCount_RNA, color = limiric)) +
    scale_y_log10() + geom_point(size = 0.5) +
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF")) +
    labs(title = "limiric") + 
    xlab("") + ylab("") +
    plot_theme
  
  miQC2 <- ggplot(df, aes(x = nf, y = nCount_RNA, color = miQC)) +
    scale_y_log10() + geom_point(size = 0.5) +
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF")) +
    labs(title = "miqc") + 
    xlab("") + ylab("") +
    plot_theme
  
  valiDrops2 <- ggplot(df, aes(x = nf, y = nCount_RNA, color = valiDrops)) +
    scale_y_log10() + geom_point(size = 0.5) +
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF")) +
   #  scale_color_manual(values = c("live" = "#D5D5D5", "dead" = "#9AA9FF")) +
    labs(title = "valiDrops") + 
    xlab("") + ylab("") +
    plot_theme
  
  ddqc2 <- ggplot(df, aes(x = nf, y = nCount_RNA, color = ddqc)) +
    scale_y_log10() + geom_point(size = 0.5) +
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF", "Uncertain" = "grey")) +
    # scale_color_manual(values = c("True" = "#D5D5D5", "False" = "#9AA9FF", "Uncertain" = "grey")) +
    labs(title = "ddqc") + 
    xlab("") + ylab("") +
    plot_theme
  
  # limiric plots ------
  
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
  
  # View output 
  limiric3 <- DimPlot(seurat_view, 
                      reduction = 'tsne', 
                      group.by = "limiric") + 
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF")) +
    NoAxes() + NoLegend() +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
  
  DropletQC3 <- DimPlot(seurat_view, 
                        reduction = 'tsne', 
                        group.by = "DropletQC") + 
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF", "empty_droplet" = "darkgrey")) +
    # scale_color_manual(values = c("cell" = "#D5D5D5", "damaged_cell" = "#9AA9FF", "empty_droplet" = "darkgrey")) +
    NoAxes() + NoLegend() +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
  
  Seurat3 <- DimPlot(seurat_view, 
                     reduction = 'tsne', 
                     group.by = "Seurat") + 
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF")) +
    NoAxes() + NoLegend() +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
  
  miQC3 <- DimPlot(seurat_view, 
                   reduction = 'tsne', 
                   group.by = "miQC") + 
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF")) +
    NoAxes() + NoLegend() +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
  
  valiDrops3 <- DimPlot(seurat_view, 
                        reduction = 'tsne', 
                        group.by = "valiDrops") + 
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF")) +
    # scale_color_manual(values = c("live" = "#D5D5D5", "dead" = "#9AA9FF")) +
    NoAxes() + NoLegend() +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
  
  ddqc3 <- DimPlot(seurat_view, 
                   reduction = 'tsne', 
                   group.by = "ddqc") + 
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF", "Uncertain" = "grey")) +
    # scale_color_manual(values = c("True" = "#D5D5D5", "False" = "#9AA9FF", "Uncertain" = "grey")) +
    NoAxes() + NoLegend() +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
  
  
  # Combine plots ------
  combined_plot <- plot_grid(
    limiric, DropletQC, Seurat, miQC, valiDrops, ddqc,
    limiric2, DropletQC2, Seurat2, miQC2, valiDrops2, ddqc2,
    limiric3, DropletQC3, Seurat3, miQC3, valiDrops3, ddqc3,
    ncol = 6, 
    rel_heights = c(1, 1)
  )
  
  combined_plot <- ggdraw(combined_plot) + 
    theme(
      plot.background = element_rect(fill = 'white', colour = NA)
    )
  
  # Save the plot
  ggsave(output_file, combined_plot, width = 24, height = 10, units = "in", dpi = 300)
  
  return(combined_plot)
}



# Same but ccount for ground truth ------
GroundBench <- function(seurat, output_file, organism) {
  
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

  
  # Create individual plots
  groundtruth <- ggplot(df, aes(x = nFeature_RNA, y = mito.ratio, color = orig.ident)) +
    geom_point(size = 0.5) +
    scale_color_manual(values = c("live" = "#D5D5D5", "dead" = "#AAD4B5")) +
    labs(title = "Truth") + 
    xlab("") + ylab("") +
    plot_theme
  
  # Create individual plots
  Seurat <- ggplot(df, aes(x = nFeature_RNA, y = mito.ratio, color = threshold)) +
    geom_point(size = 0.5) +
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF")) +
    labs(title = "Threshold") + 
    xlab("") + ylab("") +
    plot_theme
  
  
  # Create individual plots
  Seurat <- ggplot(df, aes(x = nFeature_RNA, y = mito.ratio, color = threshold)) +
    geom_point(size = 0.5) +
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF")) +
    labs(title = "Threshold") + 
    xlab("") + ylab("") +
    plot_theme
  
  DropletQC <- ggplot(df, aes(x = nFeature_RNA, y = mito.ratio, color = DropletQC)) +
    geom_point(size = 0.5) +
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF", "empty_droplet" = "darkgrey")) +
    labs(title = "DropletQC") + 
    xlab("") + ylab("") +
    plot_theme
  
  limiric <- ggplot(df, aes(x = nFeature_RNA, y = mito.ratio, color = limiric)) +
    geom_point(size = 0.5) +
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF")) +
    labs(title = "limiric") + 
    xlab("") + ylab("") +
    plot_theme
  
  miQC <- ggplot(df, aes(x = nFeature_RNA, y = mito.ratio, color = miQC)) +
    geom_point(size = 0.5) +
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF")) +
    labs(title = "miqc") + 
    xlab("") + ylab("") +
    plot_theme
  
  valiDrops <- ggplot(df, aes(x = nFeature_RNA, y = mito.ratio, color = valiDrops)) +
    geom_point(size = 0.5) +
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF")) +
    #scale_color_manual(values = c("live" = "#D5D5D5", "dead" = "#9AA9FF")) +
    labs(title = "valiDrops") + 
    xlab("") + ylab("") +
    plot_theme
  
  ddqc <- ggplot(df, aes(x = nFeature_RNA, y = mito.ratio, color = ddqc)) +
    geom_point(size = 0.5) +
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF", "Uncertain" = "grey")) +
    # scale_color_manual(values = c("True" = "#D5D5D5", "False" = "#9AA9FF", "Uncertain" = "grey")) +
    labs(title = "ddqc") + 
    xlab("") + ylab("") +
    plot_theme
  
  # Nuclear fraction plots ------
  
  # Create individual plots
  groundtruth2 <- ggplot(df, aes(x = nf, y = nCount_RNA, color = orig.ident)) +
    scale_y_log10() + geom_point(size = 0.5) +
    scale_color_manual(values = c("live" = "#D5D5D5", "dead" = "#AAD4B5")) +
    labs(title = "Truth") + 
    xlab("") + ylab("") +
    plot_theme
  
  Seurat2 <- ggplot(df, aes(x = nf, y = nCount_RNA, color = threshold)) +
    scale_y_log10() + geom_point(size = 0.5) +
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF")) +
    labs(title = "Seurat") + 
    xlab("") + ylab("") +
    plot_theme
  
  DropletQC2 <- ggplot(df, aes(x = nf, y = nCount_RNA, color = DropletQC)) +
    scale_y_log10() + geom_point(size = 0.5) +
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF", "empty_droplet" = "darkgrey")) +
    # scale_color_manual(values = c("cell" = "#D5D5D5", "damaged_cell" = "#9AA9FF", "empty_droplet" = "darkgrey")) +
    labs(title = "DropletQC") +
    xlab("") + ylab("") +
    plot_theme
  
  limiric2 <- ggplot(df, aes(x = nf, y = nCount_RNA, color = limiric)) +
    scale_y_log10() + geom_point(size = 0.5) +
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF")) +
    labs(title = "limiric") + 
    xlab("") + ylab("") +
    plot_theme
  
  miQC2 <- ggplot(df, aes(x = nf, y = nCount_RNA, color = miQC)) +
    scale_y_log10() + geom_point(size = 0.5) +
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF")) +
    labs(title = "miqc") + 
    xlab("") + ylab("") +
    plot_theme
  
  valiDrops2 <- ggplot(df, aes(x = nf, y = nCount_RNA, color = valiDrops)) +
    scale_y_log10() + geom_point(size = 0.5) +
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF")) +
    #  scale_color_manual(values = c("live" = "#D5D5D5", "dead" = "#9AA9FF")) +
    labs(title = "valiDrops") + 
    xlab("") + ylab("") +
    plot_theme
  
  ddqc2 <- ggplot(df, aes(x = nf, y = nCount_RNA, color = ddqc)) +
    scale_y_log10() + geom_point(size = 0.5) +
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF", "Uncertain" = "grey")) +
    # scale_color_manual(values = c("True" = "#D5D5D5", "False" = "#9AA9FF", "Uncertain" = "grey")) +
    labs(title = "ddqc") + 
    xlab("") + ylab("") +
    plot_theme
  
  # limiric plots ------
  
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
  seurat_view$orig.ident <- seurat$orig.ident
  seurat_view$limiric <- seurat$limiric
  seurat_view$DropletQC <- seurat$DropletQC
  seurat_view$Seurat <- seurat$threshold
  seurat_view$miQC <-  seurat$miQC
  seurat_view$valiDrops <- seurat$valiDrops
  seurat_view$ddqc <- seurat$ddqc
  
  # View output 
  groundtruth3 <- DimPlot(seurat_view, 
                      reduction = 'tsne', 
                      group.by = "orig.ident") + 
    scale_color_manual(values = c("live" = "#D5D5D5", "dead" = "#AAD4B5")) +
    NoAxes() + NoLegend() +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
  
  limiric3 <- DimPlot(seurat_view, 
                      reduction = 'tsne', 
                      group.by = "limiric") + 
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF")) +
    NoAxes() + NoLegend() +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
  
  DropletQC3 <- DimPlot(seurat_view, 
                        reduction = 'tsne', 
                        group.by = "DropletQC") + 
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF", "empty_droplet" = "darkgrey")) +
    # scale_color_manual(values = c("cell" = "#D5D5D5", "damaged_cell" = "#9AA9FF", "empty_droplet" = "darkgrey")) +
    NoAxes() + NoLegend() +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
  
  Seurat3 <- DimPlot(seurat_view, 
                     reduction = 'tsne', 
                     group.by = "Seurat") + 
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF")) +
    NoAxes() + NoLegend() +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
  
  miQC3 <- DimPlot(seurat_view, 
                   reduction = 'tsne', 
                   group.by = "miQC") + 
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF")) +
    NoAxes() + NoLegend() +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
  
  valiDrops3 <- DimPlot(seurat_view, 
                        reduction = 'tsne', 
                        group.by = "valiDrops") + 
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF")) +
    # scale_color_manual(values = c("live" = "#D5D5D5", "dead" = "#9AA9FF")) +
    NoAxes() + NoLegend() +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
  
  ddqc3 <- DimPlot(seurat_view, 
                   reduction = 'tsne', 
                   group.by = "ddqc") + 
    scale_color_manual(values = c("cell" = "#D5D5D5", "damaged" = "#9AA9FF", "Uncertain" = "grey")) +
    # scale_color_manual(values = c("True" = "#D5D5D5", "False" = "#9AA9FF", "Uncertain" = "grey")) +
    NoAxes() + NoLegend() +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
  
  
  # Combine plots ------
  combined_plot <- plot_grid(
    groundtruth, limiric, DropletQC, Seurat, miQC, valiDrops, ddqc,
    groundtruth2, limiric2, DropletQC2, Seurat2, miQC2, valiDrops2, ddqc2,
    groundtruth3, limiric3, DropletQC3, Seurat3, miQC3, valiDrops3, ddqc3,
    ncol = 7, 
    rel_heights = c(1, 1)
  )
  
  combined_plot <- ggdraw(combined_plot) + 
    theme(
      plot.background = element_rect(fill = 'white', colour = NA)
    )
  
  # Save the plot
  ggsave(output_file, combined_plot, width = 24, height = 10, units = "in", dpi = 300)
  
  return(combined_plot)
}
