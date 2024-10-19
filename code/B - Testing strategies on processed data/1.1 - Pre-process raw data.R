# SCRIPT CONTEXT 
#
# To keep tool comparisons fair, the same processed count matrix for each sample is 
# used as input into all tool testing. In R, the processed count matrices are stored
# in 'Seurat' objects. For tools not run in R, like `ddqc` and `EnsembleKQC`, count matrices 
# are converted to a csv. In all cases, matrices house identical counts. 
#
#
# Preparations: 
# 1. Gene annotations for the human and mouse genome are obtained trhough bioMart
# 2. Converting Ensembl gene symbols (ENSG0000....) to standard HGNC symbols (CD79A) where needed
#
#
# Processing done: 
# 1. SoupX ambient RNA correction, requiring raw count matrices. 
#    Note: If raw count path not present, this step will be skipped.
# 2. Filtering out features that are only expressed in fewer than 9 cells 
# 3. Calculating nuclear fraction (or alternative) scores
# 4. Calculating QC metrics mt.percent, rb.percent, MALAT1 percent 
# 5. DoubletFinder to filter potential doublets (first filtering done)
# 6. Remove red blood cells if present 
# 7. Extract processed count matrix for non-R tools (ddqc)


#-------------------------------------------------------------------------------
# PREPARATIONS
#-------------------------------------------------------------------------------

# Load libraries -------

packages <- c("AnnotationHub", "biomaRt", "cowplot", "devtools", "dplyr", "DoubletFinder", 
              "ensembldb", "ggplot2", "glmGamPoi", "Matrix", "png", "presto", "scuttle", "Seurat", "SingleCellExperiment", "SoupX", "tidyr")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    library(pkg)
  }
}


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
# GENE SYMBOL FUNCTION DEFINED 
#-------------------------------------------------------------------------------

# Converting to traditional gene symbols ------

# Input processed files (Zenodo) have non-conventional gene names that need to be changed

renamegenes <- function(processed_path, project_name) {
  
  # Read in processed object
  sample <- readRDS(processed_path)
  sample <- as.Seurat(sample)
  
  # Filter and edit the column names -----
  sample@meta.data <- sample@meta.data[, c("cell_status", "id", "nCount_originalexp", "nFeature_originalexp")]
  
  # Edit the genes from ensembl to hgnc -----
  
  # Correct gene symbols to recognizable format 
  sample@assays$RNA <- sample@assays$originalexp
  counts <- sample@assays$RNA$counts
  
  # Create Emsembl to HGNC gene naming map 
  ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
  bm <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), 
              values = rownames(counts), 
              mart = ensembl)
  
  # Perform mapping
  hgnc.symbols <- bm$hgnc_symbol[match(rownames(counts), bm$ensembl_gene_id)] # matching ENSG0... genes to gene names like DUX4
  counts <- as.matrix(counts)
  rownames(counts) <- hgnc.symbols
  
  # Remove empty and duplicated rows 
  counts <- counts[!is.na(rownames(counts)), ]
  counts <- counts[(rownames(counts) != ""), ]
  counts <- counts[unique(rownames(counts)), ]
  
  length(rownames(counts))
  length(unique(rownames(counts)))
  
  seurat <- CreateSeuratObject(counts = counts, 
                               project = project_name)
  
  seurat$orig.ident <- project_name
  
  return(seurat) # not saving because it needs to be processed, this is temporary form 
}

# Run the conversions 
GM18507_control <- renamegenes(project_name = "GM18507_control", processed_path = "/home/alicen/Projects/limiric/test_data/tumour_groundtruth/batchfx-for-zenodo/TENX049_SA928_001_sceset_v3_raw.rds")
GM18507_dying <- renamegenes(project_name = "GM18507_dying", processed_path = "/home/alicen/Projects/limiric/test_data/tumour_groundtruth/batchfx-for-zenodo/TENX049_SA928_002_sceset_v3_raw.rds")
GM18507_dead <- renamegenes(project_name = "GM18507_dead", processed_path = "/home/alicen/Projects/limiric/test_data/tumour_groundtruth/batchfx-for-zenodo/TENX049_SA928_003_sceset_v3_raw.rds")
PDX_control <- renamegenes(project_name = "PDX_control", processed_path = "/home/alicen/Projects/limiric/test_data/tumour_groundtruth/batchfx-for-zenodo/TENX019_SA604X7XB02089_002_sceset_v3_raw.rds")
PDX_dead <- renamegenes(project_name = "PDX_dead", processed_path = "/home/alicen/Projects/limiric/test_data/tumour_groundtruth/batchfx-for-zenodo/TENX019_SA604X7XB02089_004_sceset_v3_raw.rds")



#-------------------------------------------------------------------------------
# PREPROCESSING FUNCTION DEFINED 
#-------------------------------------------------------------------------------


# Function for processing ------

preprocess <- function(
    project_name,           # string with dataset identifier 
    organism = "Hsap",      # Either Hsap or Mmus
    hemo_threshold = 50,    # Adjustable for more stringent filtering of red blood cells 
    SoupX = TRUE,           # TRUE or FALSE 
    Seurat = NULL,          # For the case where no counts of any kind, Seurat object can be used 
    scale_nf_malat1 = TRUE, # For ground truth, no scaling before merging 0 - 1 based on both control and dead cells
    raw_path = NULL,        # path to .gz raw output of STARsolo
    filtered_path = NULL,   # path to .gz filtered output of STARsolo
    velocyto_path = NULL,   # path to .gz velocyto output of STARsolo
    output_path             # where all outputs are saved 
){
  
  message("Begin pre-processing for ", project_name, "...")
  
  # 1. SoupX Correction and 2. Feature filtering -------
  
  # Some cases may want this step to be skipped, making it optional 
  if (SoupX)  {
    
    # Using Seurat read in the matrices from the STARsolo output (must be zipped input files)
    table_of_counts <- suppressWarnings(Read10X(filtered_path))
    table_of_droplets <- suppressWarnings(Read10X(raw_path))
    
    # Create the soup channel (sc)
    sc <- SoupChannel(table_of_droplets, table_of_counts, calcSoupProfile = FALSE)
    
    # Estimate the contamination
    sc <- estimateSoup(sc)
    
    # Use Seurat to cluster the filtered matrix, although not essential it is recommended to get better estimations
    seurat_soup <- suppressWarnings(CreateSeuratObject(table_of_counts, min.cells = 10))
    seurat_soup <- suppressWarnings(SCTransform(seurat_soup, verbose = FALSE) %>%
                                      RunPCA(verbose = FALSE) %>%
                                      RunUMAP(dims = 1:30, verbose = FALSE) %>%
                                      FindNeighbors(dims = 1:30, verbose = FALSE) %>%
                                      FindClusters(verbose = FALSE))
    
    # Adding the cluster embeddings to the SoupX object
    meta.data <- seurat_soup@meta.data
    umap.embedding <- seurat_soup@reductions$umap@cell.embeddings
    sc <- suppressWarnings(setClusters(sc, setNames(meta.data$seurat_clusters, rownames(meta.data))))
    sc <- suppressWarnings(setDR(sc, umap.embedding, c("UMAP_1", "UMAP_2")))
    
    # With defined clusters, run the main SoupX function to calculate the contamination fraction rho where rho E (0, 1) and the closer to 1, the more contaminated
    sc <- autoEstCont(sc, verbose = FALSE, doPlot = FALSE)
    rho_estimate <- sc$metaData$rho[1] # record the rho in PreProcess output
    
    # Silencing the output ...
    # Open a connection to a temporary file for writing
    tmp_conn <- file(tempfile(), open = "wt")
    
    # Redirect standard output and messages to the temporary file
    sink(tmp_conn)
    sink(tmp_conn, type = "message")
    
    # Call the (now silenced) verbose function
    # Output integer matrix of soup-corrected reads (unzipped output) where contaminated reads are removed
    adj.matrix <- suppressWarnings(adjustCounts(sc, roundToInt = T))
    
    # Reset output redirection
    sink(NULL)
    sink(NULL, type = "message")
    
    # Close the connection
    close(tmp_conn)
    
    # Output results in Seurat object for continued workflow
    seurat <- suppressWarnings(CreateSeuratObject(counts = adj.matrix, # SoupX corrected count matrix
                                                  min.cells = 10,  # A feature must be expressed in at least 10 cells to be retained in the count matrix     
                                                  project = project_name))
    
    # Normalize counts 
    seurat <- NormalizeData(seurat)
    
    # Terminal output : 
    cell_number <- length(Cells(seurat))
    message("\u2714  Ambient correction complete for ", cell_number, " cells")
    message("\u2714  SoupX contamination estimate of ", rho_estimate)
    
  }
  if (SoupX == FALSE & is.null(Seurat)) {
    
    # Using Seurat read in the matrices from the STARsolo output for filtered (TOC) and raw (TOD) counts  (must be zipped input files)
    table_of_counts <- suppressWarnings(Read10X(filtered_path))
    
    seurat <- CreateSeuratObject(counts = table_of_counts, 
                             min.cells = 5, # A feature must be expressed in at least 5 cells to be retained in the count matrix
                             project = project_name)
    
    # Estimates for feature filtering 
    # 0 : ~ 60 000 features
    # 1 : ~ 30 000 features
    # 2 : ~ 25 000 features
    # 3 : ~ 25 000 features
    # 4 : ~ 22 000 features 
    # 5 : ~ 20 000 features 
    
    # Normalize counts 
    seurat <- NormalizeData(seurat)
    
    cell_number <- length(Cells(seurat))

    message("\u2714 ", cell_number, " cells detected, ambient correction skipped")

  }
  if (SoupX == FALSE & !is.null(Seurat)){ 
    
    seurat <- Seurat
  
  }
  
  
  # 3. Calculate nuclear fraction -----
  
  if (!is.null(velocyto_path)){
    
    message("Calculating nuclear fraction scores...")
    
    # Define parameters for file paths and gene names
    spliced_file   <- file.path(velocyto_path, "spliced.mtx.gz")
    barcodes_file  <- file.path(velocyto_path, "barcodes.tsv.gz")
    features_file  <- file.path(velocyto_path, "features.tsv.gz")
    unspliced_file <- file.path(velocyto_path, "unspliced.mtx.gz")
    
    # Add data from spliced and unspliced counts to make Seurat objects
    spliced <- ReadMtx(mtx      = spliced_file,
                       cells    = barcodes_file,
                       features = features_file)
    
    spliced <- suppressWarnings(CreateSeuratObject(
      counts    = spliced,
      project   = "spliced"))
    
    unspliced <- ReadMtx(mtx      = unspliced_file,
                         cells    = barcodes_file,
                         features = features_file)
    
    unspliced <- suppressWarnings(CreateSeuratObject(
      counts    = unspliced,
      project   = "unspliced"))
    
    # Calculate nuclear fraction
    ExonSum   <- Matrix::colSums(spliced[['RNA']]$counts)   # summing over all the genes for each cell (1 reading per cell)
    IntronSum <- Matrix::colSums(unspliced[['RNA']]$counts)
    NuclearFraction <- IntronSum / (ExonSum + IntronSum)
    nf <- data.frame(barcode = rownames(unspliced@meta.data), nf = NuclearFraction)

    # Add to seurat (ensure row order (barcode) correct)
    seurat$nf <- nf$nf[match(rownames(seurat@meta.data), nf$barcode)]
    
  }
  
  
  # Inspired by Z. Clarke and G. Bager, 2024; using nuclear lnRNA gene expression MALAT1 and NEAT1 as substitutes for nf scores if they absent (since we still want to test DropletQC) 
  # Note: This is run for all cases to see, when nf is calculated, how how much it differs from this score. 
  
  if (organism == "Hsap"){ seurat$nf_malat1 <- FetchData(seurat, vars = "MALAT1", layer = "counts")}
  if (organism == "Mmus"){ seurat$nf_malat1 <- FetchData(seurat, vars = "Malat1", layer = "counts")}
    
  # min-max normalization: scales values so the min value maps to 0 and max maps to 1 (make the expression values more like nf scores)
  if (scale_nf_malat1) {
    
   seurat$nf_malat1 <- (seurat$nf_malat1 - min(seurat$nf_malat1)) / 
    (max(seurat$nf_malat1) - min(seurat$nf_malat1))
  
  }
  
  # 4. Calculate standard quality control metrics ----- 
  
  message("Calculating quality scores...")

  # Retrieve the corresponding annotations for the organism of interest 
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
  
  # Calculate the feature percentages and normalised expressions 
  seurat$mt.percent <- PercentageFeatureSet(
    object   = seurat,
    features = intersect(mt_genes, rownames(seurat@assays$RNA)),
    assay    = "RNA"
  ) 

  # mt_genes_present <- intersect(mt_genes, rownames(seurat@assays$RNA))
  # mt_gene_expression <- FetchData(seurat, vars = mt_genes_present, layer = counts)
  # seurat$mt <- apply(mt_gene_expression, 1, mean)
  
  seurat$rb.percent <- PercentageFeatureSet(
    object   = seurat,
    features = intersect(rb_genes, rownames(seurat@assays$RNA)),
    assay    = "RNA"
  ) 
  
  # rb_genes_present <- intersect(rb_genes, rownames(seurat@assays$RNA))
  # rb_gene_expression <- FetchData(seurat, vars = rb_genes_present)
  # seurat$rb <- apply(rb_gene_expression, 1, mean)
  
  seurat$malat1.percent <- PercentageFeatureSet(
    object   = seurat,
    features = malat1,
    assay    = "RNA"
  ) 
  
  seurat$malat1 <- FetchData(seurat, vars = malat1)
  
  
  # 5. DoubletFinder filtering -----
  
  message("Doublet Finder running...")

  # Prepare seurat object (required) for DropletQC
  seurat_DF <- NormalizeData(seurat, verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    FindNeighbors(dims = 1:10, verbose = FALSE) %>%
    FindClusters(verbose = FALSE) %>%
    RunUMAP(dims = 1:10, verbose = FALSE)

  # Open a connection to a temporary file for writing (output of DoubletFinder)
  tmp_conn <- file(tempfile(), open = "wt")

  # Redirect standard output and messages to the temporary file
  sink(tmp_conn)
  sink(tmp_conn, type = "message")

  # pK Identification (no ground-truth)
  sweep.res.list <- paramSweep(seurat_DF, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)

  pK <- bcmvn %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
   dplyr::filter(BCmetric == max(BCmetric)) %>%
   dplyr::select(pK)
  pK <- as.numeric(as.character(pK[[1]]))


  # Homotypic Doublet Proportion Estimate
  annotations <- seurat_DF@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.075 * nrow(seurat_DF@meta.data))  # Assuming 7.5% doublet formation rate (can be changed for dataset specificity)
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  nExp_poi.adj <- nExp_poi.adj - (nExp_poi.adj * (1 / 7))

# Run DF itself
  seurat_DF <- doubletFinder(
    seurat_DF,
    PCs = 1:10,
    pN = 0.25,
    pK = pK,
    nExp = nExp_poi.adj,
    reuse.pANN = FALSE,
  )

  # Reset output redirection
  sink(NULL)
  sink(NULL, type = "message")

  # Close the connection
  close(tmp_conn)


  # Meta data for actual Seurat object
  df_col_name <- grep("^DF\\.classifications", colnames(seurat_DF@meta.data), value = TRUE)
  doublet_classifications <- seurat_DF@meta.data[[df_col_name]]
  seurat$DF <- doublet_classifications

  # Remove doublets before running any of the detection tools
  seurat <- subset(seurat, DF == "Singlet")
  seurat$DF <- NULL
 
  
  # 6. Remove red blood cells -----
  
  message("Filtering red blood cells...")
  
  if (organism == "Hsap") { hemo_gene <- c("HBA1", "HBA2", "HBB") } # hemoglobin subunit genes 
  if (organism == "Mmus") { hemo_gene <- c("Hba-a1", "Hba-a2", "Hbb-bt", "Hbb-bs") } # lower case & different 
  
  # Calculate proportion of hemo gene expression & store as meta data
  seurat$RBC <- PercentageFeatureSet(
    object   = seurat,
    features = intersect(hemo_gene, rownames(seurat@assays$RNA)),
    assay    = "RNA"
  )
  
  # Filter and remove column 
  seurat$RBC <- ifelse(seurat$RBC >=  hemo_threshold, "RBC", "non-RBC")

  # Calculate the number of RBCs 
  unfiltered_RBC <- length(Cells(seurat))
  seurat <- subset(seurat, RBC == "non-RBC")
  RBC_number <- unfiltered_RBC - length(Cells(seurat))
  message("\u2714 ", RBC_number, " red blood cells removed")
  
  seurat$RBC <- NULL
  
    
  # 7. Create output to run in python tools ----
  
  # Extract count matrix
  counts <- seurat@assays$RNA$counts
  matrix <- CreateSeuratObject(counts = counts, assay = "RNA")
  matrix <- suppressWarnings(as.matrix(matrix@assays$RNA$counts))
  write.csv(matrix, file = paste0(output_path, "/python/input/", project_name, "_matrix_data.csv"))
  message("\u2714  Input for python prepared")
  
  # Save processed Seurat object 
  saveRDS(seurat, file = paste0(output_path, "/R_objects/", project_name, "_processed.rds"))

  return(seurat)
  
}



# Samples 
# A549 
cellline_A549 <- preprocess(project_name = "A549",
                           organism = "Hsap",
                           SoupX = TRUE,
                           raw_path = "/home/alicen/Projects/limiric/test_data/cellline_A549/raw/",
                           filtered_path = "/home/alicen/Projects/limiric/test_data/cellline_A549/filtered/",
                           velocyto_path = "/home/alicen/Projects/limiric/test_data/cellline_A549/velocyto/",
                           output_path = "/home/alicen/Projects/ReviewArticle"
)

# HCT-116 
cellline_HCT116 <- preprocess(project_name = "HCT116",
                             organism = "Hsap",
                             SoupX = FALSE,
                             raw_path = "/home/alicen/Projects/limiric/test_data/cellline_HCT-116/raw/",
                             filtered_path = "/home/alicen/Projects/limiric/test_data/cellline_HCT-116/filtered/",
                             velocyto_path = "/home/alicen/Projects/limiric/test_data/cellline_HCT-116/velocyto/",
                             output_path = "/home/alicen/Projects/ReviewArticle"
)


# Jurkat 
cellline_jurkat <- preprocess(project_name = "jurkat",
                             organism = "Hsap",
                             SoupX = FALSE,
                             raw_path = "/home/alicen/Projects/limiric/test_data/cellline_jurkat/raw/",
                             filtered_path = "/home/alicen/Projects/limiric/test_data/cellline_jurkat/filtered/",
                             velocyto_path = "/home/alicen/Projects/limiric/test_data/cellline_jurkat/velocyto/",
                             output_path = "/home/alicen/Projects/ReviewArticle"
)



# 3 Diseased tissue 
diseased_liver <- preprocess(project_name = "dLiver",
                            organism = "Hsap",
                            raw_path = "/home/alicen/Projects/limiric/test_data/diseased_liver/raw/",
                            filtered_path = "/home/alicen/Projects/limiric/test_data/diseased_liver/filtered/",
                            velocyto_path = "/home/alicen/Projects/limiric/test_data/diseased_liver/velocyto/",
                            output_path = "/home/alicen/Projects/ReviewArticle"
)


diseased_lung <- preprocess(project_name = "dLung",
                           organism = "Hsap",
                           raw_path = "/home/alicen/Projects/limiric/test_data/diseased_lung/raw/",
                           filtered_path = "/home/alicen/Projects/limiric/test_data/diseased_lung/filtered/",
                           velocyto_path = "/home/alicen/Projects/limiric/test_data/diseased_lung/velocyto/",
                           output_path = "/home/alicen/Projects/ReviewArticle"
)


diseased_PBMC <- preprocess(project_name = "dPBMC",
                           organism = "Hsap",
                           raw_path = "/home/alicen/Projects/limiric/test_data/diseased_PBMC/raw/",
                           filtered_path = "/home/alicen/Projects/limiric/test_data/diseased_PBMC/filtered/",
                           velocyto_path = "/home/alicen/Projects/limiric/test_data/diseased_PBMC/velocyto/",
                           output_path = "/home/alicen/Projects/ReviewArticle"
)


# 3 healthy tissues 
healthy_liver <- preprocess(project_name = "hLiver",
                           organism = "Hsap",
                           raw_path = "/home/alicen/Projects/limiric/test_data/healthy_liver/raw/",
                           filtered_path = "/home/alicen/Projects/limiric/test_data/healthy_liver/filtered/",
                           velocyto_path = "/home/alicen/Projects/limiric/test_data/healthy_liver/velocyto/",
                           output_path = "/home/alicen/Projects/ReviewArticle"
)

healthy_lung <- preprocess(project_name = "hLung",
                          organism = "Hsap",
                          SoupX = FALSE,
                          raw_path = "/home/alicen/Projects/limiric/test_data/healthy_lung/raw/",
                          filtered_path = "/home/alicen/Projects/limiric/test_data/healthy_lung/filtered/",
                          velocyto_path = "/home/alicen/Projects/limiric/test_data/healthy_lung/velocyto/",
                          output_path = "/home/alicen/Projects/ReviewArticle"
)

healthy_PBMC <- preprocess(project_name = "hPBMC",
                          organism = "Hsap",
                          raw_path = "/home/alicen/Projects/limiric/test_data/healthy_PBMC/raw/",
                          filtered_path = "/home/alicen/Projects/limiric/test_data/healthy_PBMC/filtered/",
                          velocyto_path = "/home/alicen/Projects/limiric/test_data/healthy_PBMC/velocyto/",
                          output_path = "/home/alicen/Projects/ReviewArticle"
)


# 3 Mouse samples
mouse_liver <- preprocess(project_name = "mLiver",
                         organism = "Mmus",
                         raw_path = "/home/alicen/Projects/limiric/test_data/mouse_liver/raw/",
                         filtered_path = "/home/alicen/Projects/limiric/test_data/mouse_liver/filtered/",
                         velocyto_path = "/home/alicen/Projects/limiric/test_data/mouse_liver/velocyto/",
                         output_path = "/home/alicen/Projects/ReviewArticle"
)

mouse_lung <- preprocess(project_name = "mLung",
                        organism = "Mmus",
                        raw_path = "/home/alicen/Projects/limiric/test_data/mouse_lung/raw/",
                        filtered_path = "/home/alicen/Projects/limiric/test_data/mouse_lung/filtered/",
                        velocyto_path = "/home/alicen/Projects/limiric/test_data/mouse_lung/velocyto/",
                        output_path = "/home/alicen/Projects/ReviewArticle"
)

mouse_PBMC <- preprocess(project_name = "mPBMC",
                        organism = "Mmus",
                        raw_path = "/home/alicen/Projects/limiric/test_data/mouse_PBMC/raw/",
                        filtered_path = "/home/alicen/Projects/limiric/test_data/mouse_PBMC/filtered/",
                        velocyto_path = "/home/alicen/Projects/limiric/test_data/mouse_PBMC/velocyto/",
                        output_path = "/home/alicen/Projects/ReviewArticle"
)


# 3 Tumor samples 
tumor_ductal <- preprocess(project_name = "ductal",
                          organism = "Hsap",
                          raw_path = "/home/alicen/Projects/limiric/test_data/tumor_ductal/raw/",
                          filtered_path = "/home/alicen/Projects/limiric/test_data/tumor_ductal/filtered/",
                          velocyto_path = "/home/alicen/Projects/limiric/test_data/tumor_ductal/velocyto/",
                          output_path = "/home/alicen/Projects/ReviewArticle"
)

tumor_glio <- preprocess(project_name = "glio",
                        organism = "Hsap",
                        raw_path = "/home/alicen/Projects/limiric/test_data/tumor_glioblastoma/raw/",
                        filtered_path = "/home/alicen/Projects/limiric/test_data/tumor_glioblastoma/filtered/",
                        velocyto_path = "/home/alicen/Projects/limiric/test_data/tumor_glioblastoma/velocyto/",
                        output_path = "/home/alicen/Projects/ReviewArticle"
)

tumor_hodgkin <- preprocess(project_name = "hodgkin",
                           organism = "Hsap",
                           raw_path = "/home/alicen/Projects/limiric/test_data/tumor_hodgkin/raw/",
                           filtered_path = "/home/alicen/Projects/limiric/test_data/tumor_hodgkin/filtered/",
                           velocyto_path = "/home/alicen/Projects/limiric/test_data/tumor_hodgkin/velocyto/",
                           output_path = "/home/alicen/Projects/ReviewArticle"
)


# Ground truth cases
HEK293_control <- preprocess(project_name = "HEK293_control",
                             organism = "Hsap",
                             scale_nf_malat1 = FALSE,
                             raw_path = "/home/alicen/Projects/limiric/test_data/ground_truth/healthy/raw/",
                             filtered_path = "/home/alicen/Projects/limiric/test_data/ground_truth/healthy/filtered/",
                             velocyto_path = "/home/alicen/Projects/limiric/test_data/ground_truth/healthy/velocyto/",
                             output_path = "/home/alicen/Projects/ReviewArticle"
)
                             
HEK293_apoptotic <- preprocess(project_name = "HEK293_apoptotic",
                             organism = "Hsap",
                             scale_nf_malat1 = FALSE,
                             raw_path = "/home/alicen/Projects/limiric/test_data/ground_truth/apoptotic/raw/",
                             filtered_path = "/home/alicen/Projects/limiric/test_data/ground_truth/apoptotic/filtered/",
                             velocyto_path = "/home/alicen/Projects/limiric/test_data/ground_truth/apoptotic/velocyto/",
                             output_path = "/home/alicen/Projects/ReviewArticle"
)

HEK293_proapoptotic <- preprocess(project_name = "HEK293_proapoptotic",
                               organism = "Hsap",
                               SoupX = FALSE,
                               scale_nf_malat1 = FALSE,
                               raw_path = "/home/alicen/Projects/limiric/test_data/ground_truth/proapoptotic/raw/",
                               filtered_path = "/home/alicen/Projects/limiric/test_data/ground_truth/proapoptotic/filtered/",
                               velocyto_path = "/home/alicen/Projects/limiric/test_data/ground_truth/proapoptotic/velocyto/",
                               output_path = "/home/alicen/Projects/ReviewArticle"
)

GM18507_control_processed <- preprocess(project_name = "GM18507_control",
                              organism = "Hsap",
                              SoupX = FALSE,
                              scale_nf_malat1 = FALSE,
                              Seurat = GM18507_control,
                              output_path = "/home/alicen/Projects/ReviewArticle"
)

GM18507_dying_processed <- preprocess(project_name = "GM18507_dying",
                                      organism = "Hsap",
                                      SoupX = FALSE,
                                      scale_nf_malat1 = FALSE,
                                      Seurat = GM18507_dying,
                                      output_path = "/home/alicen/Projects/ReviewArticle"
)

GM18507_dead_processed <- preprocess(project_name = "GM18507_dead",
                                      organism = "Hsap",
                                      SoupX = FALSE,
                                     scale_nf_malat1 = FALSE,
                                      Seurat = GM18507_dead,
                                      output_path = "/home/alicen/Projects/ReviewArticle"
)

PDX_control_processed <- preprocess(project_name = "PDX_control",
                                    organism = "Hsap",
                                    SoupX = FALSE,
                                    scale_nf_malat1 = FALSE,
                                    Seurat = PDX_control,
                                    output_path = "/home/alicen/Projects/ReviewArticle"
)

PDX_dead_processed <- preprocess(project_name = "PDX_dead",
                                 organism = "Hsap",
                                 SoupX = FALSE,
                                 scale_nf_malat1 = FALSE,
                                 Seurat = PDX_dead,
                                 output_path = "/home/alicen/Projects/ReviewArticle"
)


### End
