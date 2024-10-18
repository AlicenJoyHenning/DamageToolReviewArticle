# Investigating damaged cell differences 

# A: WHAT IS the difference between the 2 dead cell populations in each case? (DGEA & PCA)
# B: When viewed WITH healthy cells, is the "hidden" damaged cell population visible in each case?
# C: When viewing ALL DEAD cell populations together, how consistently do the two types of dead cell populations 
divide?


library(Seurat)
library(biomaRt) # gene name conversion
library(limiric)
library(dplyr)
library(ggplot2)

library(loupeR)
library(hdf5r)
library(scCustomize)
# library(Matrix)


# Read in data ------

# TENX049_SA928_001_sceset
dead_SA928 <- 
readRDS("/home/alicen/Projects/limiric/test_data/tumour_groundtruth/batchfx-for-zenodo/TENX049_SA928_003_sceset_v3_raw.rds")
dead_SA928 <- as.Seurat(dead_SA928) 

dying_SA928 <- 
readRDS("/home/alicen/Projects/limiric/test_data/tumour_groundtruth/batchfx-for-zenodo/TENX049_SA928_002_sceset_v3_raw.rds")
dying_SA928 <- as.Seurat(dying_SA928)
live_SA928 <- 
readRDS("/home/alicen/Projects/limiric/test_data/tumour_groundtruth/batchfx-for-zenodo/TENX049_SA928_001_sceset_v3_raw.rds")
live_SA928 <- as.Seurat(live_SA928)


# SA604 ground truth datasets 
live_SA604_rep1 <- 
readRDS("/home/alicen/Projects/limiric/test_data/tumour_groundtruth/batchfx-for-zenodo/TENX019_SA604X7XB02089_001_sceset_v3_raw.rds")
live_SA604_rep1 <- as.Seurat(live_SA604_rep1)
live_SA604_rep2 <- 
readRDS("/home/alicen/Projects/limiric/test_data/tumour_groundtruth/batchfx-for-zenodo/TENX019_SA604X7XB02089_002_sceset_v3_raw.rds")
live_SA604_rep2 <- as.Seurat(live_SA604_rep2)
live_SA604_rep3 <- 
readRDS("/home/alicen/Projects/limiric/test_data/tumour_groundtruth/batchfx-for-zenodo/TENX019_SA604X7XB02089_003_sceset_v3_raw.rds")
live_SA604_rep3 <- as.Seurat(live_SA604_rep3)
dead_SA604 <- 
readRDS("/home/alicen/Projects/limiric/test_data/tumour_groundtruth/batchfx-for-zenodo/TENX019_SA604X7XB02089_004_sceset_v3_raw.rds")
dead_SA604 <- as.Seurat(dead_SA604)


# HEK293 data 
apoptotic <- Read10X("/home/alicen/Projects/limiric/test_data/ground_truth/apoptotic/filtered/")
apoptotic <- CreateSeuratObject(counts = apoptotic,
                                assay  = "RNA",
                                min.cells = 1)
dim(apoptotic)

pro_apoptotic <- Read10X("/home/alicen/Projects/limiric/test_data/ground_truth/proapoptotic/filtered/")
pro_apoptotic <- CreateSeuratObject(counts = pro_apoptotic,
                                    assay  = "RNA",
                                    min.cells = 1)

healthy <- Read10X("/home/alicen/Projects/limiric/test_data/ground_truth/healthy/filtered/")
healthy <- CreateSeuratObject(counts = healthy,
                              assay  = "RNA",
                              min.cells = 1)

# Correct non-standard gene naming in zenodo cases and remove excess information  ----

# Function to edit meta.data columns

edit_metadata <- function(dataset) {
  
  # Edit the column names -----
  dataset@meta.data <- dataset@meta.data[, c("cell_status", "id", "nCount_originalexp", "nFeature_originalexp")]
  
  # Edit the genes from ensembl to hgnc -----
  
  # Correct gene symbols to recognizable format 
  dataset@assays$RNA <- dataset@assays$originalexp
  countsData <- dataset@assays$RNA$counts
  
  # Create Emsembl to HGNC gene naming map 
  ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
  bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), 
              values = rownames(countsData), 
              mart = ensembl)
  
  # Perform mapping
  hgnc.symbols <- bm$hgnc_symbol[match(rownames(countsData), bm$ensembl_gene_id)] # matching ENSG0... genes to 
gene names like DUX4
  countsData <- as.matrix(countsData)
  rownames(countsData) <- hgnc.symbols
  
  # Remove empty and duplicated rows 
  countsData <- countsData[!is.na(rownames(countsData)), ]
  countsData <- countsData[(rownames(countsData) != ""), ]
  countsData <- countsData[unique(rownames(countsData)), ]
  
  length(rownames(countsData))
  length(unique(rownames(countsData)))
  
  dataset <- CreateSeuratObject(counts = countsData, 
                                project = "groundtruth")
  
  return(dataset)
}

dead_SA928 <- edit_metadata(dead_SA928)
dying_SA928 <- edit_metadata(dying_SA928)
live_SA928 <- edit_metadata(live_SA928)

# live_SA604_rep1 <- edit_metadata(live_SA604_rep1) # different treatment 
live_SA604_rep2 <- edit_metadata(live_SA604_rep2)
# live_SA604_rep3 <- edit_metadata(live_SA604_rep3) # same treatment but very few cells 
dead_SA604 <- edit_metadata(dead_SA604)


# A. Dimensionality reduction for individual damaged cells, interested to see how (whether) dead cells form 
clusters ----

# Get labels 
data("human_annotations", package = "limiric")
annotations <- human_annotations

# Get gene annotations for mitochondrial (MT) & ribosomal (RB) genes
mt_gene_annotations <- annotations[grep("MT-", annotations$gene_name, perl=TRUE),]
mt_gene_annotations <- mt_gene_annotations[grepl("protein_coding", mt_gene_annotations$gene_biotype, perl=TRUE),]
mt_genes <- mt_gene_annotations %>% pull(gene_name)

# isolate ribosomal genes
rps_genes <- annotations[grep("RPS", annotations$gene_name, perl=TRUE),]
rps_genes <- rps_genes[grepl("protein_coding", rps_genes$gene_biotype, perl=TRUE),]
rps_genes <- rps_genes %>% pull(gene_name)
rpl_genes <- annotations[grep("RPL", annotations$gene_name, perl=TRUE),]
rpl_genes <- rpl_genes[grepl("protein_coding", rpl_genes$gene_biotype, perl=TRUE),]
rpl_genes <- rpl_genes %>% pull(gene_name)
rb_genes  <- c(rps_genes, rpl_genes)

# combine mt and rb genes
mt_rb_genes <- c(mt_genes, rb_genes)
mt_rb_genes <- unique(mt_rb_genes)


# Function to reduced dimensions of JUST damaged cell populations 
lets_see <- function(seurat){
  
  # Define mitochondrial expression
  seurat$mt.percent <- PercentageFeatureSet(
    object   = seurat,
    features = intersect(mt_genes, rownames(seurat@assays$RNA)),
    assay    = "RNA"
  )
  
  # Define ribosomal expression
  seurat$rb.percent <- PercentageFeatureSet(
    object   = seurat,
    features = intersect(rb_genes, rownames(seurat@assays$RNA)),
    assay    = "RNA"
  )

  
  
  # Reduce based on mt & rb genes only (this occurs in a separate Seurat object (limiric))
  limiric <- subset(seurat, features = intersect(mt_rb_genes, rownames(seurat@assays$RNA)))
  
  limiric <- NormalizeData(limiric) %>%
            FindVariableFeatures() %>%
            ScaleData() %>%
            RunPCA() %>%
            FindNeighbors(dims = 1:30) %>%
            FindClusters() %>%
            RunUMAP(dims = 1:30)
  
  
  
  # Visualise
  colours_full <- c("#6C79F0","#A1A9F5","#ADDBB6", "#72BC7B", "#70B1D2","#9FCBE1", 
"#C5E7A7","#ABDC7E","#C9AFDF","#9D71B3",
                    "#6C79F0","#A1A9F5","#ADDBB6", "#72BC7B", "#70B1D2","#9FCBE1", 
"#C5E7A7","#ABDC7E","#C9AFDF","#9D71B3")
  
  # For label consistency with PreProcess plots
  cells <- length(Cells(seurat))
  
  plot_clusters <- DimPlot(limiric,
                           pt.size = 1,
                           label = TRUE,
                           label.box = TRUE,
                           label.size = 5,
                           label.color = "white",
                           repel = TRUE) +
    scale_color_manual(values = colours_full) +
    scale_fill_manual(values = colours_full) + 
    NoAxes() + NoLegend() +
    xlab("UMAP 1") + ylab("UMAP 2") +
    theme(plot.title = element_text(hjust = 0.5, size = rel(1), face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, vjust = 1),
          panel.border = element_rect(colour = "black", fill=NA, linewidth =1))
  
  

  
  mt_plot <- FeaturePlot(limiric,
                           pt.size = 1,
                           features = "mt.percent") +
    NoAxes() + NoLegend() +
    xlab("UMAP 1") + ylab("UMAP 2") +
    theme(plot.title = element_text(hjust = 0.5, size = rel(1), face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, vjust = 1),
          panel.border = element_rect(colour = "black", fill=NA, linewidth =1))
  
  
  rb_plot <- FeaturePlot(limiric,
                         pt.size = 1,
                         features = "rb.percent") +
    NoAxes() + NoLegend() +
    xlab("UMAP 1") + ylab("UMAP 2") +
    theme(plot.title = element_text(hjust = 0.5, size = rel(1), face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, vjust = 1),
          panel.border = element_rect(colour = "black", fill=NA, linewidth =1))
  
  
  plot <- plot_clusters + mt_plot + rb_plot
  
  
  # Calculate cell cycle scores 
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  seurat <- NormalizeData(seurat) %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>% 
    CellCycleScoring(s.features = intersect(s.genes, rownames(seurat@assays$RNA)), 
                     g2m.features = intersect(g2m.genes, rownames(seurat@assays$RNA)))
  
  
  # Ensure the output seurat object has mtrb reduced embedding but RNA content of original 
  limiric$Phase <- seurat$Phase
  limiric$mt.percent <- seurat$mt.percent
  limiric$rb.percent <- seurat$rb.percent
  limiric@assays$RNA <- seurat@assays$RNA
  
  
  return(list(seurat = limiric, plot = plot))
         
}

# Run the function
dead_SA928_output <- lets_see(seurat = dead_SA928)
dying_SA928_output <- lets_see(dying_SA928)
dead_SA604_output <- lets_see(dead_SA604)
apoptotic_output <- lets_see(apoptotic)
pro_apoptotic_output <- lets_see(pro_apoptotic)

# Checking clustering 
dead_SA928_output$seurat
dying_SA928_output$plot 
dead_SA604_output$plot 
apoptotic_output$plot 
pro_apoptotic_output$plot 


  
# Save LoupeBrowser objects for viewing (individual can be done straight)
dead_SA928_output$seurat$mitochondrial <- dead_SA928_output$seurat$mt.percent
create_loupe_from_seurat(dead_SA928_output$seurat,
                         output_dir = "/home/alicen/Projects/limiric/groundtruth/Loupe/reduced_by_mtrb_genes",
                         output_name = "SA928_dead_mtrb",
                         force = TRUE
)

create_loupe_from_seurat(dying_SA928$seurat,
                         output_dir = "/home/alicen/Projects/limiric/groundtruth/Loupe",
                         output_name = "SA928_dying",
                         force = TRUE
)

create_loupe_from_seurat(dead_SA604_output$seurat,
                         output_dir = "/home/alicen/Projects/limiric/groundtruth/Loupe/reduced_by_mtrb_genes/",
                         output_name = "SA604_dead",
                         force = TRUE
)


create_loupe_from_seurat(apoptotic$seurat,
                         output_dir = "/home/alicen/Projects/limiric/groundtruth/Loupe",
                         output_name = "HEK293_apoptotic",
                         force = TRUE
)

create_loupe_from_seurat(pro_apoptotic$seurat,
                         output_dir = "/home/alicen/Projects/limiric/groundtruth/Loupe",
                         output_name = "HEK293_proapoptotic",
                         force = TRUE
)



# like percent mt and ribo, also percent apoptotic? 


# B. Merge live and dead cells into groups corresponding to datasets (cell line) of origin -----

# Edit to make it according to mitochondrial & ribosomal genes 

# Function to merge, integrate, and visualise groups 
merge_and_see <- function(input_list = SA928, 
                          mtrb_reduce = FALSE,
                          sample_IDs = c("live", "dying", "dead"),
                          project_name = "SA928",
                          output_dir = "/home/alicen/Projects/limiric/groundtruth/R_objects/"
){

  # Merge the samples (accounting for # of samples)
  
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
  merged_obj <- merge(
      x = first_obj,
      y = remaining_objs,
      add.cell.ids = sample_IDs,  
      project = project_name 
    )
  
  
  if (mtrb_reduce) {
    
    # Storage 
    full_merged_obj <- merged_obj
  
    # Reduce based on mt & rb genes only (this occurs in a separate Seurat object (limiric))
    merged_obj <- subset(merged_obj, features = intersect(mt_rb_genes, rownames(merged_obj@assays$RNA)))
    
    
  }
  
  
  # Attempt to integrate 
  integrated_obj <- NormalizeData(merged_obj) %>%
                    FindVariableFeatures() %>%
                    ScaleData() %>%
                    RunPCA() %>%
                    IntegrateLayers(method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated")

  
  # re-join layers after integration
  integrated_obj[["RNA"]] <- JoinLayers(integrated_obj[["RNA"]])
  
  # Create meta data column for sample name 
  integrated_obj$orig.ident <- rownames(integrated_obj@meta.data)
  integrated_obj$orig.ident <- sub("_.*", "", integrated_obj$orig.ident)
  
  
  
  
  # Dimensionality reduction for visualisation (UMAP)
  seurat <- FindNeighbors(integrated_obj, reduction = "integrated", dims = 1:30) %>%
    FindClusters(reduction = "integrated", dims = 1:30) %>%
    RunUMAP(reduction = "integrated", dims = 1:30)

  
  # Visualise
  colours_full <- c("#6C79F0","#A1A9F5","#ADDBB6", "#72BC7B", "#70B1D2","#9FCBE1", "#C5E7A7","#C9AFDF","#9D71B3")
                    
  # For label consistency with PreProcess plots
  cells <- length(Cells(seurat))
  message("Sample ", project_name, " has ", cells, " cells.")
  
  plot_clusters <- DimPlot(seurat,
                           group.by = "orig.ident",
                           pt.size = 1,
                           label = TRUE,
                           label.box = TRUE,
                           label.size = 5,
                           label.color = "white",
                           repel = TRUE) +
    scale_color_manual(values = colours_full) +
    scale_fill_manual(values = colours_full) + 
    NoAxes() + NoLegend() +
    xlab("UMAP 1") + ylab("UMAP 2") +
    theme(plot.title = element_text(hjust = 0.5, size = rel(1), face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, vjust = 1),
          panel.border = element_rect(colour = "black", fill=NA, linewidth =1))
  
  
  if (mtrb_reduce) {
    
    # Define mitochondrial expression
    seurat$mt.percent <- PercentageFeatureSet(
      object   = full_merged_obj,
      features = intersect(mt_genes, rownames(full_merged_obj@assays$RNA)),
      assay    = "RNA")
    
    
    # Define ribosomal expression
    seurat$rb.percent <- PercentageFeatureSet(
      object   = full_merged_obj,
      features = intersect(rb_genes, rownames(full_merged_obj@assays$RNA)),
      assay    = "RNA")

    seurat@assays$RNA <- integrated_obj@assays$RNA
    
  }
  
  else { 
    
  # Define mitochondrial expression
  seurat$mt.percent <- PercentageFeatureSet(
    object   = seurat,
    features = intersect(mt_genes, rownames(seurat@assays$RNA)),
    assay    = "RNA"
  )
  
  # Define ribosomal expression
  seurat$rb.percent <- PercentageFeatureSet(
    object   = seurat,
    features = intersect(rb_genes, rownames(seurat@assays$RNA)),
    assay    = "RNA"
  )
  
  }
  
  mt_plot <- FeaturePlot(seurat,
                         pt.size = 1,
                         features = "mt.percent") +
    NoAxes() + NoLegend() +
    xlab("UMAP 1") + ylab("UMAP 2") +
    theme(plot.title = element_text(hjust = 0.5, size = rel(1), face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, vjust = 1),
          panel.border = element_rect(colour = "black", fill=NA, linewidth =1))
  
  
  rb_plot <- FeaturePlot(seurat,
                         pt.size = 1,
                         features = "rb.percent") +
    NoAxes() + NoLegend() +
    xlab("UMAP 1") + ylab("UMAP 2") +
    theme(plot.title = element_text(hjust = 0.5, size = rel(1), face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, vjust = 1),
          panel.border = element_rect(colour = "black", fill=NA, linewidth =1))
  
  
  plot <- plot_clusters + mt_plot + rb_plot
  
  saveRDS(seurat, 
          paste0(output_dir, project_name, ".rds"))
  
  return(list(seurat = seurat,
              plot = plot))
  
}


# Cases
SA928_all <- list(live_SA928, dying_SA928, dead_SA928)
SA928_livendead <- list(live_SA928, dead_SA928)
SA928_livendying <- list(live_SA928, dying_SA928)
SA604 <- list(live_SA604_rep2, dead_SA604)
HEK293_all <- list(healthy, pro_apoptotic, apoptotic)
HEK293_livenapoptotic <-  list(healthy, apoptotic)
HEK293_livenproapoptotic <-  list(healthy, pro_apoptotic)



# Run
SA928_all <- merge_and_see(SA928_all, mtrb_reduce = TRUE, c("live", "dying", "dead"), "SA928_all")
SA928_livendead <- merge_and_see(SA928_livendead, mtrb_reduce = FALSE, c("live", "dead"), "SA928_livendead")
SA928_livendying <- merge_and_see(SA928_livendying, mtrb_reduce = FALSE, c("live", "dying"), "SA928_livendying")
SA604 <- merge_and_see(SA604, mtrb_reduce = FALSE, c("live", "dead"), "SA694_livendead")
HEK293_all <- merge_and_see(HEK293_all, mtrb_reduce = TRUE, c("live", "pro_apoptotic", "apoptotic"), "HEK293_all")
HEK293_livenapoptotic <- merge_and_see(HEK293_livenapoptotic, mtrb_reduce = TRUE, c("live", "apoptotic"), 
"HEK293_livenapoptotic")
HEK293_livenproapoptotic <- merge_and_see(HEK293_livenproapoptotic, mtrb_reduce = TRUE, c("live", 
"pro_apoptotic"), "HEK293_livenproapoptotic")


SA928_all$plot
SA928_livendead$plot
SA928_livendying$plot
SA604$plot
HEK293_livenapoptotic$plot
HEK293_livenproapoptotic$plot
  
  

# Create meta data column with the cluster information of the damaged cells found in isolation  

SA928_dead_labels <- data.frame(clusters = dead_SA928$seurat$seurat_clusters)
SA928_dead_labels$cells <- rownames(SA928_dead_labels)
SA928_dead_labels$cells <- paste0("dead_", SA928_dead_labels$cells)
SA928_dead_labels$clusters <- paste(SA928_dead_labels$clusters, "dead", sep = "_")

SA928_dying_labels <- data.frame(clusters = dying_SA928$seurat$seurat_clusters)
SA928_dying_labels$cells <- rownames(SA928_dying_labels)
SA928_dying_labels$cells <- paste0("dying_", SA928_dying_labels$cells)
SA928_dying_labels$clusters <- paste(SA928_dying_labels$clusters, "dying", sep = "_")

SA604_labels <- data.frame(clusters = dead_SA604$seurat$seurat_clusters)
SA604_labels$cells <- rownames(SA604_labels)
SA604_labels$cells <- paste0("dead_", SA604_labels$cells)
SA604_labels$clusters <- paste(SA604_labels$clusters, "dead", sep = "_")

apoptotic_labels <- data.frame(clusters = apoptotic$seurat$seurat_clusters)
apoptotic_labels$cells <- rownames(apoptotic_labels)
apoptotic_labels$cells <- paste0("apoptotic_", apoptotic_labels$cells)
apoptotic_labels$clusters <- paste(apoptotic_labels$clusters, "apoptotic", sep = "_")

pro_apoptotic_labels <- data.frame(clusters = pro_apoptotic$seurat$seurat_clusters)
pro_apoptotic_labels$cells <- rownames(pro_apoptotic_labels)
pro_apoptotic_labels$cells <- paste0("pro_apoptotic_", pro_apoptotic_labels$cells)
pro_apoptotic_labels$clusters <- paste(pro_apoptotic_labels$clusters, "pro_apoptotic", sep = "_")


###

# Create meta data column for labels from isolated damaged cells 
# All cases for SA928
SA928_all$seurat$labels <- SA928_dead_labels$clusters[match(rownames(SA928_all$seurat@meta.data), 
SA928_dead_labels$cells)]
SA928_all$seurat$labels <- ifelse(SA928_all$seurat$orig.ident == "dying", # condition: if cells match
                                  SA928_dying_labels$clusters[match(rownames(SA928_all$seurat@meta.data), 
SA928_dying_labels$cells)], # use these labels else 
                                  SA928_all$seurat$labels)
SA928_all$seurat$labels <- ifelse(SA928_all$seurat$orig.ident == "live", # condition: if cells match
                                  "live",
                                  SA928_all$seurat$labels)

# Dead
SA928_livendead$seurat$labels <- SA928_dead_labels$clusters[match(rownames(SA928_livendead$seurat@meta.data), 
SA928_dead_labels$cells)]
SA928_livendead$seurat$labels <- ifelse(SA928_livendead$seurat$orig.ident == "live", # condition: if cells match
                                  "live",
                                  SA928_livendead$seurat$labels) 


SA928_livendying$seurat$labels <- SA928_dying_labels$clusters[match(rownames(SA928_livendying$seurat@meta.data), 
SA928_dying_labels$cells)]
SA928_livendying$seurat$labels <- ifelse(SA928_livendying$seurat$orig.ident == "live", # condition: if cells match
                                        "live",
                                        SA928_livendying$seurat$labels) 


# SA604 cases 
SA604$seurat$labels <- SA604_labels$clusters[match(rownames(SA604$seurat@meta.data), SA604_labels$cells)]
SA604$seurat$labels <- ifelse(SA604$seurat$orig.ident == "live", # condition: if cells match
                                         "live",
                              SA604$seurat$labels) 

# HEK293 cases 
HEK293_all$seurat$labels <- apoptotic_labels$clusters[match(rownames(HEK293_all$seurat@meta.data), 
apoptotic_labels$cells)]
HEK293_all$seurat$labels <- ifelse(HEK293_all$seurat$orig.ident == "pro", # condition: if cells match
                                   pro_apoptotic_labels$clusters[match(rownames(HEK293_all$seurat@meta.data), 
pro_apoptotic_labels$cells)], # use these labels else 
                                   HEK293_all$seurat$labels)
HEK293_all$seurat$labels <- ifelse(HEK293_all$seurat$orig.ident == "live", # condition: if cells match
                                  "live",
                                  HEK293_all$seurat$labels)

# Apoptotic cases 
HEK293_livenapoptotic$seurat$labels <- 
apoptotic_labels$clusters[match(rownames(HEK293_livenapoptotic$seurat@meta.data), apoptotic_labels$cells)]
HEK293_livenapoptotic$seurat$labels <- ifelse(HEK293_livenapoptotic$seurat$orig.ident == "live", # condition: if 
cells match
                                   "live",
                                   HEK293_livenapoptotic$seurat$labels)


# Pro-apoptotic cases
HEK293_livenproapoptotic$seurat$labels <- 
pro_apoptotic_labels$clusters[match(rownames(HEK293_livenproapoptotic$seurat@meta.data), 
pro_apoptotic_labels$cells)]
HEK293_livenproapoptotic$seurat$labels <- ifelse(HEK293_livenproapoptotic$seurat$orig.ident == "live", # 
condition: if cells match
                                              "live",
                                              HEK293_livenproapoptotic$seurat$labels)


# Create output Loupe Browser objects 
create_loupe_from_seurat(SA928_all$seurat,
                         output_dir = "/home/alicen/Projects/limiric/groundtruth/Loupe/reduced_by_mtrb_genes/",
                         output_name = "integrated_SA928_all_mtrb",
                         force = TRUE
)

create_loupe_from_seurat(SA928_livendead$seurat,
                         output_dir = "/home/alicen/Projects/limiric/groundtruth/Loupe/reduced_by_mtrb_genes/",
                         output_name = "integrated_SA928_livendead_mtrb",
                         force = TRUE
)

create_loupe_from_seurat(SA928_livendying$seurat,
                         output_dir = "/home/alicen/Projects/limiric/groundtruth/Loupe/reduced_by_mtrb_genes/",
                         output_name = "integrated_SA928_livendying_mtrb",
                         force = TRUE
)

create_loupe_from_seurat(SA604$seurat,
                         output_dir = "/home/alicen/Projects/limiric/groundtruth/Loupe/reduced_by_mtrb_genes/",
                         output_name = "integrated_SA604_mtrb",
                         force = TRUE
)

create_loupe_from_seurat(HEK293_all$seurat,
                         output_dir = "/home/alicen/Projects/limiric/groundtruth/Loupe",
                         output_name = "integrated_HEK293_all_mtrb",
                         force = TRUE
)

create_loupe_from_seurat(HEK293_livenapoptotic$seurat,
                         output_dir = "/home/alicen/Projects/limiric/groundtruth/Loupe/reduced_by_mtrb_genes/",
                         output_name = "integrated_HEK293_livenapoptotic_mtrb",
                         force = TRUE
)

create_loupe_from_seurat(HEK293_livenproapoptotic$seurat,
                         output_dir = "/home/alicen/Projects/limiric/groundtruth/Loupe/reduced_by_mtrb_genes/",
                         output_name = "integrated_HEK293_livenproapoptotic_mtrb",
                         force = TRUE
)



# C. Merge all damaged datasets -----

damaged_input <-  list(dead_SA604$seurat, dead_SA928$seurat, dying_SA928$seurat, apoptotic$seurat, 
pro_apoptotic$seurat)

# Run
damaged <- merge_and_see(input_list = damaged_input, 
                         sample_IDs = c("SA604_dead", "SA928_dead", "SA928_dying", "HEK293_apoptotic", 
"HEK293_proapoptotic"), 
                         project_name = "damaged")

# Editing labels fit better into the damaged output
SA604_labels$cells <- paste0("SA604_", SA604_labels$cells)
SA604_labels$clusters <- paste0("SA604_", SA604_labels$clusters)
SA928_dead_labels$cells <- paste0("SA928_", SA928_dead_labels$cells)
SA928_dead_labels$clusters <- paste0("SA928_", SA928_dead_labels$clusters)
SA928_dying_labels$cells <- paste0("SA928_", SA928_dying_labels$cells)
SA928_dying_labels$clusters <- paste0("SA928_", SA928_dying_labels$clusters)
apoptotic_labels$cells <- paste0("HEK293_", apoptotic_labels$cells)
apoptotic_labels$clusters <- paste0("HEK293_", apoptotic_labels$clusters)

pro_apoptotic_labels <- data.frame(clusters = pro_apoptotic$seurat$seurat_clusters)
pro_apoptotic_labels$cells <- rownames(pro_apoptotic_labels)
pro_apoptotic_labels$cells <- paste0("proapoptotic_", pro_apoptotic_labels$cells)
pro_apoptotic_labels$clusters <- paste(pro_apoptotic_labels$clusters, "pro_apoptotic", sep = "_")
pro_apoptotic_labels$cells <- paste0("HEK293_", pro_apoptotic_labels$cells)
pro_apoptotic_labels$clusters <- paste0("HEK293_", pro_apoptotic_labels$clusters)


  
# Merge into one column 
damaged$seurat$labels <- SA604_labels$clusters[match(rownames(damaged$seurat@meta.data), SA604_labels$cells)]


damaged$seurat$labels <- ifelse(damaged$seurat$orig.ident == "SA928" & is.na(damaged$seurat$labels), 
                                SA928_dead_labels$clusters[match(rownames(damaged$seurat@meta.data), 
SA928_dead_labels$cells)], 
                                damaged$seurat$labels)

damaged$seurat$labels <- ifelse(damaged$seurat$orig.ident == "SA928" & is.na(damaged$seurat$labels), 
                                SA928_dying_labels$clusters[match(rownames(damaged$seurat@meta.data), 
SA928_dying_labels$cells)], 
                                damaged$seurat$labels)

damaged$seurat$labels <- ifelse(damaged$seurat$orig.ident == "HEK293" & is.na(damaged$seurat$labels), 
                                apoptotic_labels$clusters[match(rownames(damaged$seurat@meta.data), 
apoptotic_labels$cells)], 
                                damaged$seurat$labels)

damaged$seurat$labels <- ifelse(damaged$seurat$orig.ident == "HEK293" & is.na(damaged$seurat$labels), 
                                pro_apoptotic_labels$clusters[match(rownames(damaged$seurat@meta.data), 
pro_apoptotic_labels$cells)], 
                                damaged$seurat$labels)


# Save as a cloupe

create_loupe_from_seurat(damaged$seurat,
                         output_dir = "/home/alicen/Projects/limiric/groundtruth/Loupe",
                         output_name = "integrated_damaged",
                         force = TRUE
)






