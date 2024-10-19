# SCRIPT CONTEXT 
#
# To keep tool comparisons fair, the identical processed count matrix for each sample is 
# used as input into all tool testing. In R, the processed count matrices are stored
# in 'Seurat' objects. For tools not run in R, like ddqc and EnsembleKQC, count matrices 
# are converted to a csv. In all cases, matrices house identical counts. 
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
#
# Further groundtruth processing done: 
# 1. Converting Ensembl gene symbols (ENSG0000....) to standard HGNC symbols (CD79A) where needed
# 2. Merging and integrating damaged and control samples 
# 
# NOTE:
# (NB!) Ground truth cases are different to non ground truth cases. Each ground truth sample 
#       for testing damaged detection strategies includes two separately sequenced samples, 
#       one for treated & sorted dead cells, and the other for untreated, sorted live cells. 
#
#      For damaged cell detection, the count matrices of the damaged and control 
#      samples are merged and integrated to resemble 'one psuedosample'. This simply
#      concatenates the count matrices without losing their original labelss. The 
#      integration serves only for downstream clustering and visualisation, the 
#      original (processed) count matrices for each sample are used for tool testing. 
#
# Note: After running this script, the user will need to run the tool in python separately (ii - Run ddqc.ipynb)


#-------------------------------------------------------------------------------
# PREPARATIONS
#-------------------------------------------------------------------------------

# Load libraries -------

packages <- c("AnnotationHub", "biomaRt", "cowplot", "devtools", "dplyr", "DoubletFinder", "DoubletFinder", 
              "DropletQC", "ensembldb", "ggplot2", "glmGamPoi", "limiric", "Matrix", "miQC", "png", "presto", "robustbase", 
              "scuttle", "Seurat", "SingleCellExperiment", "SoupX", "tidyr", "valiDrops")

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
# PREPROCESSING FUNCTION DEFINED 
#-------------------------------------------------------------------------------


# Function for processing ------

preprocess <- function(
    project_name,        # string with dataset identifier 
    organism = "Hsap",   # Either Hsap or Mmus
    hemo_threshold = 50, # Adjustable for more stringent filtering of red blood cells 
    SoupX = TRUE,        # TRUE or FALSE 
    raw_path,            # path to .gz raw output of STARsolo
    filtered_path,       # path to .gz filtered output of STARsolo
    velocyto_path,       # path to .gz velocyto output of STARsolo
    output_path          # where all outputs are saved 
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
  else {
    
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
  
  
  # 3. Calculate nuclear fraction -----
  
  message("Calculating nuclear fraction scores...")
  
  if (!is.null(velocyto_path)){
    
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
  
  if (organism == "Hsap"){ seurat$nf_malat1 <- FetchData(seurat, vars = "MALAT1")}
  if (organism == "Mmus"){ seurat$nf_malat1 <- FetchData(seurat, vars = "Malat1")}
    
  # min-max normalization: scales values so the min value maps to 0 and max maps to 1 (make the expression values more like nf scores)
  seurat$nf_malat1 <- (seurat$nf_malat1 - min(seurat$nf_malat1)) / 
    (max(seurat$nf_malat1) - min(seurat$nf_malat1))
  
  
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
      filter(grepl("mt-", gene_name)) %>% 
      pull(gene_name)
    
    # Isolate ribosomal genes (RPS and RPL)
    rb_genes <- annotations %>%
      filter(grepl("^rsp|^rpl", gene_name)) %>%
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

  mt_genes_present <- intersect(mt_genes, rownames(seurat@assays$RNA))
  mt_gene_expression <- FetchData(seurat, vars = mt_genes_present)
  seurat$mt <- apply(mt_gene_expression, 1, mean)
  
  seurat$rb.percent <- PercentageFeatureSet(
    object   = seurat,
    features = intersect(rb_genes, rownames(seurat@assays$RNA)),
    assay    = "RNA"
  ) 
  
  rb_genes_present <- intersect(rb_genes, rownames(seurat@assays$RNA))
  rb_gene_expression <- FetchData(seurat, vars = rb_genes_present)
  seurat$rb <- apply(rb_gene_expression, 1, mean)
  
  
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
  write.csv(matrix, file = paste0(output_path, "/ddqc/", project_name, "_ddqc_matrix_data.csv"))
  message("\u2714  Input for ddqc prepared")
  
  # Save processed Seurat object 
  saveRDS(seurat, file = paste0(output_path, "/R_objects/", project_name, "_processed.rds"))

  return(seurat)
  
}

# Test 
test <- preprocess(project_name = "test", 
                   organism = "Hsap",
                   SoupX = FALSE,
                   raw_path = "/home/alicen/Projects/limiric/test_data/cellline_A549/raw/", 
                   filtered_path = "/home/alicen/Projects/limiric/test_data/cellline_A549/filtered/", 
                   velocyto_path = "/home/alicen/Projects/limiric/test_data/cellline_A549/velocyto/", 
                   output_path = "/home/alicen/Projects/ReviewArticle")

# Samples 
# A549 
cellline_A549 <- benchmark(project_name = "A549",
                           organism = "Hsap",
                           SoupX = TRUE,
                           raw_path = "/home/alicen/Projects/limiric/test_data/cellline_A549/raw/",
                           filtered_path = "/home/alicen/Projects/limiric/test_data/cellline_A549/filtered/",
                           velocyto_path = "/home/alicen/Projects/limiric/test_data/cellline_A549/velocyto/",
                           output_path = "/home/alicen/Projects/ReviewArticle"
)

# HCT-116 
cellline_HCT116 <- benchmark(project_name = "HCT116",
                             organism = "Hsap",
                             SoupX = FALSE,
                             raw_path = "/home/alicen/Projects/limiric/test_data/cellline_HCT-116/raw/",
                             filtered_path = "/home/alicen/Projects/limiric/test_data/cellline_HCT-116/filtered/",
                             velocyto_path = "/home/alicen/Projects/limiric/test_data/cellline_HCT-116/velocyto/",
                             output_path = "/home/alicen/Projects/ReviewArticle"
)


# Jurkat 
cellline_jurkat <- benchmark(project_name = "jurkat",
                             organism = "Hsap",
                             SoupX = FALSE,
                             raw_path = "/home/alicen/Projects/limiric/test_data/cellline_jurkat/raw/",
                             filtered_path = "/home/alicen/Projects/limiric/test_data/cellline_jurkat/filtered/",
                             velocyto_path = "/home/alicen/Projects/limiric/test_data/cellline_jurkat/velocyto/",
                             output_path = "/home/alicen/Projects/ReviewArticle"
)



# 3 Diseased tissue 

diseased_liver <- benchmark(project_name = "dLiver",
                            organism = "Hsap",
                            raw_path = "/home/alicen/Projects/limiric/test_data/diseased_liver/raw/",
                            filtered_path = "/home/alicen/Projects/limiric/test_data/diseased_liver/filtered/",
                            velocyto_path = "/home/alicen/Projects/limiric/test_data/diseased_liver/velocyto/",
                            output_path = "/home/alicen/Projects/ReviewArticle"
)


diseased_lung <- benchmark(project_name = "dLung",
                           organism = "Hsap",
                           raw_path = "/home/alicen/Projects/limiric/test_data/diseased_lung/raw/",
                           filtered_path = "/home/alicen/Projects/limiric/test_data/diseased_lung/filtered/",
                           velocyto_path = "/home/alicen/Projects/limiric/test_data/diseased_lung/velocyto/",
                           output_path = "/home/alicen/Projects/ReviewArticle"
)


diseased_PBMC <- benchmark(project_name = "dPBMC",
                           organism = "Hsap",
                           raw_path = "/home/alicen/Projects/limiric/test_data/diseased_PBMC/raw/",
                           filtered_path = "/home/alicen/Projects/limiric/test_data/diseased_PBMC/filtered/",
                           velocyto_path = "/home/alicen/Projects/limiric/test_data/diseased_PBMC/velocyto/",
                           output_path = "/home/alicen/Projects/ReviewArticle"
)


# 3 healthy tissues 
healthy_liver <- benchmark(project_name = "hLiver",
                           organism = "Hsap",
                           raw_path = "/home/alicen/Projects/limiric/test_data/healthy_liver/raw/",
                           filtered_path = "/home/alicen/Projects/limiric/test_data/healthy_liver/filtered/",
                           velocyto_path = "/home/alicen/Projects/limiric/test_data/healthy_liver/velocyto/",
                           output_path = "/home/alicen/Projects/ReviewArticle"
)

healthy_lung <- benchmark(project_name = "hLung",
                          organism = "Hsap",
                          SoupX = FALSE,=
                          raw_path = "/home/alicen/Projects/limiric/test_data/healthy_lung/raw/",
                          filtered_path = "/home/alicen/Projects/limiric/test_data/healthy_lung/filtered/",
                          velocyto_path = "/home/alicen/Projects/limiric/test_data/healthy_lung/velocyto/",
                          output_path = "/home/alicen/Projects/ReviewArticle"
)

healthy_PBMC <- benchmark(project_name = "hPBMC",
                          organism = "Hsap",
                          raw_path = "/home/alicen/Projects/limiric/test_data/healthy_PBMC/raw/",
                          filtered_path = "/home/alicen/Projects/limiric/test_data/healthy_PBMC/filtered/",
                          velocyto_path = "/home/alicen/Projects/limiric/test_data/healthy_PBMC/velocyto/",
                          output_path = "/home/alicen/Projects/ReviewArticle"
)


# 3 Mouse samples
mouse_liver <- benchmark(project_name = "mLiver",
                         organism = "Mmus",
                         raw_path = "/home/alicen/Projects/limiric/test_data/mouse_liver/raw/",
                         filtered_path = "/home/alicen/Projects/limiric/test_data/mouse_liver/filtered/",
                         velocyto_path = "/home/alicen/Projects/limiric/test_data/mouse_liver/velocyto/",
                         output_path = "/home/alicen/Projects/ReviewArticle"
)

mouse_lung <- benchmark(project_name = "mLung",
                        organism = "Mmus",
                        raw_path = "/home/alicen/Projects/limiric/test_data/mouse_lung/raw/",
                        filtered_path = "/home/alicen/Projects/limiric/test_data/mouse_lung/filtered/",
                        velocyto_path = "/home/alicen/Projects/limiric/test_data/mouse_lung/velocyto/",
                        output_path = "/home/alicen/Projects/ReviewArticle"
)


# PBMC
mouse_PBMC <- benchmark(project_name = "mPBMC",
                        organism = "Mmus",
                        raw_path = "/home/alicen/Projects/limiric/test_data/mouse_PBMC/raw/",
                        filtered_path = "/home/alicen/Projects/limiric/test_data/mouse_PBMC/filtered/",
                        velocyto_path = "/home/alicen/Projects/limiric/test_data/mouse_PBMC/velocyto/",
                        ddqc_path = "/home/alicen/Projects/limiric/benchmarking/ddqc_output/mouse_PBMC_ddqc_output.csv",
                        output_path = "/home/alicen/Projects/limiric/damage_left_behind_analysis"
)


# 3 Tumor samples -------
# Ductal 
tumor_ductal <- benchmark(project_name = "ductal",
                          organism = "Hsap",
                          cluster_ranks = 5,
                          raw_path = "/home/alicen/Projects/limiric/test_data/tumor_ductal/raw/",
                          filtered_path = "/home/alicen/Projects/limiric/test_data/tumor_ductal/filtered/",
                          velocyto_path = "/home/alicen/Projects/limiric/test_data/tumor_ductal/velocyto/",
                          ddqc_path = "/home/alicen/Projects/limiric/benchmarking/ddqc_output/tumor_ductal_ddqc_output.csv",
                          output_path = "/home/alicen/Projects/limiric/damage_left_behind_analysis")

# Glio
tumor_glio <- benchmark(project_name = "glio",
                        organism = "Hsap",
                        raw_path = "/home/alicen/Projects/limiric/test_data/tumor_glioblastoma/raw/",
                        filtered_path = "/home/alicen/Projects/limiric/test_data/tumor_glioblastoma/filtered/",
                        velocyto_path = "/home/alicen/Projects/limiric/test_data/tumor_glioblastoma/velocyto/",
                        ddqc_path = "/home/alicen/Projects/limiric/benchmarking/ddqc_output/tumor_glio_ddqc_output.csv",
                        output_path = "/home/alicen/Projects/limiric/damage_left_behind_analysis"
)

# Hodgkin
tumor_hodgkin <- benchmark(project_name = "hodgkin",
                           organism = "Hsap",
                           #  cluster_ranks = 2,
                           raw_path = "/home/alicen/Projects/limiric/test_data/tumor_hodgkin/raw/",
                           filtered_path = "/home/alicen/Projects/limiric/test_data/tumor_hodgkin/filtered/",
                           velocyto_path = "/home/alicen/Projects/limiric/test_data/tumor_hodgkin/velocyto/",
                           ddqc_path = "/home/alicen/Projects/limiric/benchmarking/ddqc_output/tumor_hodgkin_ddqc_output.csv",
                           output_path = "/home/alicen/Projects/limiric/damage_left_behind_analysis"
)






 ggplot(test@meta.data, aes_string(x = "nf", y = "nf_malat1")) +
    geom_point() +  # Plot points
    geom_smooth(method = "lm", color = "blue") +  # Add line of best fit
    theme_classic() 
  




#-------------------------------------------------------------------------------
# GENE SYMBOL FUNCTION DEFINED 
#-------------------------------------------------------------------------------

# Function for converting gene symbols ------


# Function to merge, integrate, and visualise groups 
mergeandintegrate <- function(input_list = SA928, 
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

