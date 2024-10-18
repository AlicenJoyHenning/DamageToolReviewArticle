# SCRIPT CONTEXT 
#
# Benchmarking damaged cell detection tools. First, input raw files are
# processed. Then, the ten damaged strategies are applied to each
# processed count matrix to obtain output damaged labels for each cell barcode.
# For all methods, this will either be 'damaged' or 'cell'.
#
# Notes: 
# 1. The function requires ddqc output from python to have been generated 
#     (i - Prepare ddqc input.R) and (ii - Run ddqc.ipynb)
#    
# 2. Ground truth cases are different to non ground truth cases. Each ground truth sample 
#    for testing damaged detection strategies includes two separately sequenced samples, 
#    one for treated & sorted dead cells and the other for untreated, sorted live cells. 
#
#    For damaged cell detection, the count matrices of the damaged and control 
#    samples are merged and integrated to resemble 'one psuedosample'. This simply
#    concatenates the count matrices without losing their original labels. The 
#    integration serves only for downstream clustering and visualisation, the 
#    original (processed) count matrices for each sample are used of tool testing. 
#
#    
# 3. Absent data from ground truth
#    The second ground truth dataset has only processed count matrices available,
#    so raw counts or velocyto counts are absent. This means SoupX ambient RNA
#    correction cannot be performed for the samples. It also means that there is 
#    no velocyto output for calculating nuclear fraction required by DropletQC. 
#    These limitations prevented the testing of DropletQC.


#-------------------------------------------------------------------------------
# PREPARATIONS
#-------------------------------------------------------------------------------

# Load libraries -------

packages <- c("cowplot", "devtools", "dplyr", "ggplot2", "glmGamPoi", "Matrix", "robustbase",
              "png", "Seurat", "SoupX", "DoubletFinder", "tidyr", "DoubletFinder", "biomaRt",
              "limiric", "miQC", "SingleCellExperiment",
              "scuttle", "presto", "valiDrops", "DropletQC")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    library(pkg)
  }
}


#-------------------------------------------------------------------------------
# FUNCTION DEFINED 
#-------------------------------------------------------------------------------

# Function for main benchmarking ------

benchmark <- function(
    project_name,        # string with dataset identifier 
    organism = "Hsap",   # Either Hsap or Mmus
    resolution = 0.5,    # Adjustable
    cluster_ranks = 1,   # Adjustable
    hemo_threshold = 50, # Adjustable
    SoupX = TRUE,        # TRUE or FALSE 
    model_method = NULL, # Adjustable, alternative to specify the model type "linear", "spline", "polynomial", or 
"one_dimensional"
    raw_path,            # path to .gz raw output of STARsolo
    filtered_path,       # path to .gz filtered output of STARsolo
    velocyto_path,       # path to .gz velocyto output of STARsolo
    ddqc_path,           # path to resulting ddqc
    output_path          # where all outputs are saved 
    
){
  
  
  message("\nBegin pre-processing for ", project_name, "...")
  
  
  # SoupX Correction -------
  
  if (SoupX)  {
    
    # Using Seurat read in the matrices from the STARsolo output for filtered (TOC) and raw (TOD) counts  (must be zipped input files)
    table_of_counts <- suppressWarnings(Read10X(filtered_path))
    
    # Only reading in if necessary  
    table_of_droplets <- suppressWarnings(Read10X(raw_path))
    
    # Create the soup channel (sc)
    sc <- SoupChannel(table_of_droplets, table_of_counts, calcSoupProfile = FALSE)
    
    # Estimate the contamination
    sc <- estimateSoup(sc)
    
    # Use Seurat to cluster the filtered matrix, although not essential it is recommended to get better 
estimations
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
    
    # With defined clusters, run the main SoupX function to calculate the contamination fraction rho where rho E 
(0, 1) and the closer to 1, the more contaminated
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
                                                  min.cells = 10,       
                                                  project = project_name))
    
    # Terminal output : 
    cell_number <- length(Cells(seurat))
    message("\u2714  Ambient correction complete for ", cell_number, " cells")
    message("\u2714  SoupX contamination estimate of ", rho_estimate)
    
  }
  
  else {
    
    # Using Seurat read in the matrices from the STARsolo output for filtered (TOC) and raw (TOD) counts  (must be 
zipped input files)
    table_of_counts <- suppressWarnings(Read10X(filtered_path))
    
    seurat <- CreateSeuratObject(counts = table_of_counts, 
                                 min.cells = 10,       
                                 project = project_name)
    
    cell_number <- length(Cells(seurat))
    
    message("\u2714 ", cell_number, " cells detected, ambient correction skipped")
    
  }
  
  # DoubletFinder filtering ------
  
  message("Doublet Finder running...")
  
  # Manual ground truth without raw data
  # seurat <- SA928_dead_live 
  # seurat <- SA928_dying_live
  # seurat <- SA604_dead_live 
  # seurat <- JoinLayers(seurat)
  
  # Prepare seurat object (required) for DropletQC
  seurat_DF <- NormalizeData(seurat, verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    FindNeighbors(dims = 1:10, verbose = FALSE) %>%
    FindClusters(verbose = FALSE) %>%
    RunUMAP(dims = 1:10, verbose = FALSE)
  
  # Open a connection to a temporary file for writing (output of DoubletFinder is RIDICULOUSLY verbose but lacks 
the ability to quieten it)
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
  nExp_poi <- round(0.075 * nrow(seurat_DF@meta.data))  # Assuming 7.5% doublet formation rate (can be changed for 
dataset specificity)
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
  
  message("Done")
  
  # Tool 1: DropletQC (needs original count matrix) ------
  
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
  
  # Now that nf has been added, remove doublets before running any of the detection tools 
  seurat <- subset(seurat, DF == "Singlet")
  
  seurat$DF <- NULL
  
  # For manual ground truth
  #seurat <- HEK293_apo_live
  #seurat <- HEK293_proapo_live
  #seurat <- JoinLayers(seurat)
  
  # Run DropletQC
  # Extract nf meta data & associated cell barcode from Seurat object
  edDf <- data.frame(nf = as.numeric(seurat$nf), umi = seurat$nCount_RNA)
  
  # Use droplet_qc function to identify empty droplets
  edresultsDf <- identify_empty_drops(edDf)
  
  # Identify damaged cells
  # Create vector of length same as cell number (runs default with no adjustment on groups of cells, need to fill 
this column requirement while getting results for that sample as a whole)
  n_cells <- length(Cells(seurat))
  cell_type <- rep(1, n_cells)
  
  dcDf <- data.frame(
    nf = as.numeric(seurat$nf),
    umi = seurat$nCount_RNA,
    cell_status = edresultsDf$cell_status,
    cell_type = cell_type
  )
  
  dcresultsDf <- identify_damaged_cells(nf_umi_ed_ct = dcDf)
  
  # Add to limiric_results 
  seurat$DropletQC <- dcresultsDf$df$cell_status[match(rownames(dcresultsDf$df), rownames(seurat@meta.data))]
  
  message("\u2714  DropletQC benchmark complete")
  
  
  # Tool 2: limiric -------
  
  # Run limiric
  limiric <- limiric(
    project_name  = project_name,
    organism      = organism,
    seurat_input  = seurat,
    filter_rbc    = FALSE,
    hemo_threshold = hemo_threshold,
    filter_output = FALSE,
    resolution = resolution,
    cluster_ranks = cluster_ranks,
    output_path  = output_path
  )
  
  # Subset seurat based. on those retained by limiric (non-rbc)
  non_rbc <- rownames(limiric@meta.data)
  seurat <- subset(seurat, cells = non_rbc)
  
  # Add output metadata to object 
  seurat$mito.ratio <- limiric$limiric.mi
  seurat$ribo.ratio <- limiric$limiric.ri
  seurat$limiric <- limiric$limiric # cell or damaged 
  

  # Tool 3: miQC -------
  
  # Create input
  counts <- seurat@assays$RNA
  sce <- SingleCellExperiment(list(counts = counts))
  rowData(sce)$featureType <- rownames(seurat@assays$RNA)
  colData(sce)$subsets_mito_percent <- seurat$mito.ratio
  mainExpName(sce) <- 'gene'
  colData(sce)$detected <- seurat$nFeature_RNA
  
  # Function to run miQC & generate stats (single cell experiment object)
  evaluate_model <- function(sce, model_type){
    
    # Run the miQC function 
    model <- mixtureModel(sce, model_type)
    
    # Calculate AIC and BIC for the fitted model
    aic_value <- AIC(model)
    bic_value <- BIC(model)
    
    return(list(model = model, 
                AIC = aic_value, 
                BIC = bic_value))
    
  }
  
  # Generate the different models 
  model_types <- c("linear", "spline", "polynomial", "one_dimensional")
  
  # Loop through each model 
  results <- list()
  for (model_type in model_types) {
    
    result <- evaluate_model(sce, model_type = model_type)
    
    results[[model_type]] <- result
    
  }
  
  # Extract AIC and BIC values for all models
  aic_values <- sapply(results, function(x) x$AIC)
  bic_values <- sapply(results, function(x) x$BIC)
  
  # Rank the models based on AIC and BIC (lowest is better)
  aic_ranks <- rank(aic_values)
  bic_ranks <- rank(bic_values)
  combined_ranks <- aic_ranks + bic_ranks
  
  # Select the best fit model (model with lowest combined rank)
  minimum_score <- min(combined_ranks) # For terminal output only 
  best_model_type <- names(which.min(combined_ranks))

  if (!is.null(model_method)) {best_model_type <- model_method}
  # best_model_type <- "linear"
  
  # Susbet based on model 
  sce_subset <- filterCells(sce, results[[best_model_type]]$model) 
  
  # Report results in Seurat object 
  miQC <- colnames(sce_subset) 
  seurat$miQC <- ifelse(rownames(seurat@meta.data) %in% miQC, "cell", "damaged")
  
  message("\n\u2714 miQC benchmark using ", best_model_type, " model")
  
  
  # Tool 4: valiDrops -------
  
  # Extract matrix from Seurat object
  expression_matrix <- GetAssayData(seurat, layer = "counts")
  
  if (organism == "Hsap") {species <- "human"}
  if (organism == "Mmus") {species <- "mouse"}
  
  # Run valiDrops internal function to get list of quality control metrics
  expression_metrics <- quality_metrics(counts = expression_matrix,
                                        species = species)
  
  
  # Run valiDrops damaged detection (label_dead valiDrops function edited to allow for this input, not previous 
ValiDrop input)
  valiDrops <- label_dead(counts = expression_matrix, 
                          metrics = expression_metrics$metrics)
  
  # Add to benchmark seurat object
  seurat$valiDrops <- valiDrops$label[match(rownames(seurat@meta.data), valiDrops$metrics$barcode)]
  
  message("\u2714  valiDrops damaged detection")
  
  
  # Tool 5: ddqc -------
  
  # Read in ddqc results 
  ddqc_results <- read.csv2(ddqc_path, sep = ',')
  ddqc_results <- ddqc_results[, c("barcodekey", "passed_qc")]
  
  # Add to Seurat 
  seurat$ddqc <- ddqc_results$passed_qc[match(ddqc_results$barcodekey, rownames(seurat@meta.data))]
  seurat$ddqc <- ifelse(is.na(seurat$ddqc), "Uncertain", seurat$ddqc)
  
  
  # Manual methods 
  
  # Make the labels comparable
  seurat$DropletQC <- ifelse(seurat$DropletQC == "damaged_cell", "damaged", seurat$DropletQC)
  seurat$valiDrops <- ifelse(seurat$valiDrops == "live", "cell", "damaged")
  seurat$ddqc <- ifelse(seurat$ddqc == "True", "cell", "damaged")
  
  # Calculate MALAT1 percentage expression 
  
  if (organism == "Hsap"){
  seurat$malat1.ratio <- PercentageFeatureSet(
    object   = seurat,
    features = c("MALAT1"),
    assay    = "RNA"
  )  
  }
  
  if (organism == "Mmus"){
    seurat$malat1.ratio <- PercentageFeatureSet(
      object   = seurat,
      features = c("Malat1"),
      assay    = "RNA"
    )  
  }
  
  
  # Extract df and add cell barcode column 
  seurat_df <- seurat@meta.data
  seurat_df$barcodes <- rownames(seurat_df)

  
  # Transform values for robust outlier detection 
  seurat_df$nCount_RNA <- log10(seurat_df$nCount_RNA)
  seurat_df$nFeature_RNA <- log10(seurat_df$nFeature_RNA)
  
  # Transform malat1 counts for cases where they are close to zero 
  if (project_name %in% c("hLiver")){
    
    seurat_df$malat1.ratio <- log10(seurat_df$malat1.ratio)
    
  }
  

  # Helper function for detecting multivariate outliers 
  identify_moutliers <- function(df){
    
    # Compute robust covariance matrix using the MCD (Minimum Covariance Determinant) method
    mcd_result <- covMcd(df)
    

    # Calculate robust Mahalanobis distances
    mahalanobis_distances <- mahalanobis(df, center = mcd_result$center, cov = mcd_result$cov)
    
    # Define a threshold for outliers (Chi-square distribution with degrees of freedom = # columns)
    threshold <- qchisq(0.975, df = ncol(df)) # 97.5% quantile of Chi-square distribution
    outliers <- mahalanobis_distances > (threshold + (threshold / 2) )
    
    # View the results
    df_results <- data.frame(barcodes = rownames(df), mahalanobis_distances, outliers)
    df_results$outliers <- ifelse(df_results$outliers == "TRUE", "damaged", "cell")
    
    return(df_results)
    
  }
    
  # Manual 1: All multivariate outliers (robustbase): UMI and feature counts & mitochondrial, ribosomal and MALAT1 
percentages -------
  
  # Define input df, run outlier detection, and add to Seurat meta data 
  manual1_df <- seurat_df[, c("nCount_RNA", "nFeature_RNA", "mito.ratio", "ribo.ratio", "malat1.ratio")]
  manual1_results <- identify_moutliers(manual1_df)
  seurat$manual_all <- manual1_results$outliers
  
  # Correct for NA if needed
  seurat$manual_all <- ifelse(is.na(seurat$manual_all), "cell", seurat$manual_all)
  
  # Manual 2: Mito/Ribo multivariate outliers (robustbase): UMI and feature counts & mitochondrial and ribosomal 
percentages -------
  
  manual2_df <- seurat_df[, c("nCount_RNA", "nFeature_RNA", "mito.ratio", "ribo.ratio")]
  manual2_results <- identify_moutliers(manual2_df)
  seurat$manual_mito_ribo <- manual2_results$outliers
 
  
  # Manual 3: Mito multivariate outliers (robustbase): UMI and feature counts & mitochondrial percentage -------
  
  manual3_df <- seurat_df[, c("nCount_RNA", "nFeature_RNA", "mito.ratio")]
  manual3_results <- identify_moutliers(manual3_df)
  seurat$manual_mito <- manual3_results$outliers
  
  
  # Manual 4: MALAT1 multivariate outliers (robustbase): UMI and feature counts & MALAT1 percentage -------
  
  manual4_df <- seurat_df[, c("nCount_RNA", "nFeature_RNA", "malat1.ratio")]
  manual4_results <- identify_moutliers(manual4_df)
  seurat$manual_malat1 <- manual4_results$outliers
  
  # Correct for NA if needed
  seurat$manual_malat1 <- ifelse(is.na(seurat$manual_malat1), "cell", seurat$manual_malat1)
  
  
  # Manual 5: Mito univariate outliers (MAD) ------
  
  manual5_df  <- seurat_df[, c("mito.ratio")]
  
  # Calculate the median and MAD 
  median_mito <- median(manual5_df, na.rm = TRUE)
  mad_mito <- mad(manual5_df, constant = 1, na.rm = TRUE)  # Use constant = 1 to follow the typical MAD definition
  
  # Define the threshold for outliers as 3 times the MAD 
  threshold_upper <- median_mito + 3 * mad_mito
  outliers_mito <- manual5_df > threshold_upper 
  

  manual5_df  <- as.data.frame(seurat_df[, c("mito.ratio")])
  rownames(manual5_df) <- rownames(seurat_df)
  manual5_results <- data.frame(barcodes = rownames(manual5_df), outliers = outliers_mito)
  manual5_results$outliers <- ifelse(manual5_results$outliers == "TRUE", "damaged", "cell")

  # Add to seurat object 
  seurat$manual_mito_isolated <- manual5_results$outliers
  
  
  # Return output -------
  
  # For manual methods where DropletQC couldn't be performed
  # seurat$DropletQC <- "cell"
  
  final_df <- summarise_results(seurat)

  write.csv(final_df, 
            paste0(output_path, "/", project_name, "_results.csv"), 
            quote = FALSE, 
            row.names = FALSE) 
  
  saveRDS(seurat, 
          paste0(output_path, "/", project_name, ".rds"), )
  
  
  # Plot the output -------
  
  # Call helper function
  plot <- BenchPlot(seurat, 
                    output_file = paste0(output_path, "/", project_name, ".png"), 
                    organism)
  
  
  # Global environment output 
  
  return(list(
    seurat_obj = seurat,
    output_df  = final_df
  ))
  
}

#-------------------------------------------------------------------------------
# FUNCTION RUN 
#-------------------------------------------------------------------------------

# Run the function on 15 non-groundtruth datasets 
# Note differences in the ground truth cases! 

# 3 Cell line -------

# A549 
cellline_A549 <- benchmark(project_name = "A549",
                           organism = "Hsap",
                           cluster_ranks = 1,
                           model_method = "linear",
                           raw_path = "/home/alicen/Projects/limiric/test_data/cellline_A549/raw/",
                           filtered_path = "/home/alicen/Projects/limiric/test_data/cellline_A549/filtered/",
                           velocyto_path = "/home/alicen/Projects/limiric/test_data/cellline_A549/velocyto/",
                           ddqc_path = 
"/home/alicen/Projects/limiric/benchmarking/ddqc_output/cellline_A549_ddqc_output.csv",
                           output_path = "/home/alicen/Projects/limiric/damage_left_behind_analysis"
)

# HCT-116 
cellline_HCT116 <- benchmark(project_name = "HCT116",
                             organism = "Hsap",
                             SoupX = FALSE,
                             raw_path = "/home/alicen/Projects/limiric/test_data/cellline_HCT-116/raw/",
                             filtered_path = "/home/alicen/Projects/limiric/test_data/cellline_HCT-116/filtered/",
                             velocyto_path = "/home/alicen/Projects/limiric/test_data/cellline_HCT-116/velocyto/",
                             ddqc_path = 
"/home/alicen/Projects/limiric/benchmarking/ddqc_output/cellline_HCT-116_ddqc_output.csv",
                             output_path = "/home/alicen/Projects/limiric/damage_left_behind_analysis"
)


# Jurkat 
cellline_jurkat <- benchmark(project_name = "jurkat",
                             organism = "Hsap",
                             SoupX = FALSE,
                             resolution = 1,
                             cluster_ranks = 2,
                             raw_path = "/home/alicen/Projects/limiric/test_data/cellline_jurkat/raw/",
                             filtered_path = "/home/alicen/Projects/limiric/test_data/cellline_jurkat/filtered/",
                             velocyto_path = "/home/alicen/Projects/limiric/test_data/cellline_jurkat/velocyto/",
                             ddqc_path = 
"/home/alicen/Projects/limiric/benchmarking/ddqc_output/cellline_jurkat_ddqc_output.csv",
                             output_path = "/home/alicen/Projects/limiric/damage_left_behind_analysis"
)



# 3 Diseased tissue -------

diseased_liver <- benchmark(project_name = "dLiver",
                            organism = "Hsap",
                            resolution = 0.8,
                            cluster_ranks = 3,
                            raw_path = "/home/alicen/Projects/limiric/test_data/diseased_liver/raw/",
                            filtered_path = "/home/alicen/Projects/limiric/test_data/diseased_liver/filtered/",
                            velocyto_path = "/home/alicen/Projects/limiric/test_data/diseased_liver/velocyto/",
                            ddqc_path = 
"/home/alicen/Projects/limiric/benchmarking/ddqc_output/diseased_liver_ddqc_output.csv",
                            output_path = "/home/alicen/Projects/limiric/damage_left_behind_analysis"
)


# Diseased lung
diseased_lung <- benchmark(project_name = "dLung",
                           organism = "Hsap",
                           raw_path = "/home/alicen/Projects/limiric/test_data/diseased_lung/raw/",
                           filtered_path = "/home/alicen/Projects/limiric/test_data/diseased_lung/filtered/",
                           velocyto_path = "/home/alicen/Projects/limiric/test_data/diseased_lung/velocyto/",
                           ddqc_path = 
"/home/alicen/Projects/limiric/benchmarking/ddqc_output/diseased_lung_ddqc_output.csv",
                           output_path = "/home/alicen/Projects/limiric/damage_left_behind_analysis"
)


# Diseased PBMC
diseased_PBMC <- benchmark(project_name = "dPBMC",
                           organism = "Hsap",
                           resolution = 1.2, 
                           cluster_ranks = 1,
                           raw_path = "/home/alicen/Projects/limiric/test_data/diseased_PBMC/raw/",
                           filtered_path = "/home/alicen/Projects/limiric/test_data/diseased_PBMC/filtered/",
                           velocyto_path = "/home/alicen/Projects/limiric/test_data/diseased_PBMC/velocyto/",
                           ddqc_path = 
"/home/alicen/Projects/limiric/benchmarking/ddqc_output/diseased_PBMC_ddqc_output.csv",
                           output_path = "/home/alicen/Projects/limiric/damage_left_behind_analysis"
)


# 3 healthy tissues -------

# Problem here : 

healthy_liver <- benchmark(project_name = "hLiver",
                           organism = "Hsap",
                           resolution = 1,
                           cluster_ranks = 7,
                           model_method = "polynomial",
                           raw_path = "/home/alicen/Projects/limiric/test_data/healthy_liver/raw/",
                           filtered_path = "/home/alicen/Projects/limiric/test_data/healthy_liver/filtered/",
                           velocyto_path = "/home/alicen/Projects/limiric/test_data/healthy_liver/velocyto/",
                           ddqc_path = 
"/home/alicen/Projects/limiric/benchmarking/ddqc_output/healthy_liver_ddqc_output.csv",
                           output_path = "/home/alicen/Projects/limiric/damage_left_behind_analysis"
)

# healthy lung
healthy_lung <- benchmark(project_name = "hLung",
                          organism = "Hsap",
                          SoupX = FALSE,
                          hemo_threshold = 1,
                          resolution = 1, 
                          cluster_ranks = 3,
                          raw_path = "/home/alicen/Projects/limiric/test_data/healthy_lung/raw/",
                          filtered_path = "/home/alicen/Projects/limiric/test_data/healthy_lung/filtered/",
                          velocyto_path = "/home/alicen/Projects/limiric/test_data/healthy_lung/velocyto/",
                          ddqc_path = 
"/home/alicen/Projects/limiric/benchmarking/ddqc_output/healthy_lung_ddqc_output.csv",
                          output_path = "/home/alicen/Projects/limiric/damage_left_behind_analysis")



# Healthy PBMC
healthy_PBMC <- benchmark(project_name = "hPBMC",
                          organism = "Hsap",
                          raw_path = "/home/alicen/Projects/limiric/test_data/healthy_PBMC/raw/",
                          filtered_path = "/home/alicen/Projects/limiric/test_data/healthy_PBMC/filtered/",
                          velocyto_path = "/home/alicen/Projects/limiric/test_data/healthy_PBMC/velocyto/",
                          ddqc_path = 
"/home/alicen/Projects/limiric/benchmarking/ddqc_output/healthy_PBMC_ddqc_output.csv",
                          output_path = "/home/alicen/Projects/limiric/damage_left_behind_analysis"
)


# 3 Mouse samples -------

# Mouse liver 
mouse_liver <- benchmark(project_name = "mLiver",
                         organism = "Mmus",
                         resolution = 1.8, 
                         cluster_ranks = 1,
                         raw_path = "/home/alicen/Projects/limiric/test_data/mouse_liver/raw/",
                         filtered_path = "/home/alicen/Projects/limiric/test_data/mouse_liver/filtered/",
                         velocyto_path = "/home/alicen/Projects/limiric/test_data/mouse_liver/velocyto/",
                         ddqc_path = 
"/home/alicen/Projects/limiric/benchmarking/ddqc_output/mouse_liver_ddqc_output.csv",
                         output_path = "/home/alicen/Projects/limiric/damage_left_behind_analysis"
)


# Mouse lung 
mouse_lung <- benchmark(project_name = "mLung",
                        organism = "Mmus",
                        raw_path = "/home/alicen/Projects/limiric/test_data/mouse_lung/raw/",
                        filtered_path = "/home/alicen/Projects/limiric/test_data/mouse_lung/filtered/",
                        velocyto_path = "/home/alicen/Projects/limiric/test_data/mouse_lung/velocyto/",
                        ddqc_path = 
"/home/alicen/Projects/limiric/benchmarking/ddqc_output/mouse_lung_ddqc_output.csv",
                        output_path = "/home/alicen/Projects/limiric/damage_left_behind_analysis"
)


# PBMC
mouse_PBMC <- benchmark(project_name = "mPBMC",
                        organism = "Mmus",
                        raw_path = "/home/alicen/Projects/limiric/test_data/mouse_PBMC/raw/",
                        filtered_path = "/home/alicen/Projects/limiric/test_data/mouse_PBMC/filtered/",
                        velocyto_path = "/home/alicen/Projects/limiric/test_data/mouse_PBMC/velocyto/",
                        ddqc_path = 
"/home/alicen/Projects/limiric/benchmarking/ddqc_output/mouse_PBMC_ddqc_output.csv",
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
                          ddqc_path = 
"/home/alicen/Projects/limiric/benchmarking/ddqc_output/tumor_ductal_ddqc_output.csv",
                          output_path = "/home/alicen/Projects/limiric/damage_left_behind_analysis")

# Glio
tumor_glio <- benchmark(project_name = "glio",
                        organism = "Hsap",
                        raw_path = "/home/alicen/Projects/limiric/test_data/tumor_glioblastoma/raw/",
                        filtered_path = "/home/alicen/Projects/limiric/test_data/tumor_glioblastoma/filtered/",
                        velocyto_path = "/home/alicen/Projects/limiric/test_data/tumor_glioblastoma/velocyto/",
                        ddqc_path = 
"/home/alicen/Projects/limiric/benchmarking/ddqc_output/tumor_glio_ddqc_output.csv",
                        output_path = "/home/alicen/Projects/limiric/damage_left_behind_analysis"
)

# Hodgkin
tumor_hodgkin <- benchmark(project_name = "hodgkin",
                           organism = "Hsap",
                          #  cluster_ranks = 2,
                           raw_path = "/home/alicen/Projects/limiric/test_data/tumor_hodgkin/raw/",
                           filtered_path = "/home/alicen/Projects/limiric/test_data/tumor_hodgkin/filtered/",
                           velocyto_path = "/home/alicen/Projects/limiric/test_data/tumor_hodgkin/velocyto/",
                           ddqc_path = 
"/home/alicen/Projects/limiric/benchmarking/ddqc_output/tumor_hodgkin_ddqc_output.csv",
                           output_path = "/home/alicen/Projects/limiric/damage_left_behind_analysis"
)



# Ground truth datasets -----

# HEK293 : apoptotic & healthy 
# Preprocess (SoupX, DoubletFinder, nf scores) then merge, then create ddqc matrix, then run benchmark with seurat 
as input
project_name = "apoptotic"
organism = "Hsap"
raw_path = "/home/alicen/Projects/limiric/test_data/ground_truth/apoptotic/raw/"
filtered_path = "/home/alicen/Projects/limiric/test_data/ground_truth/apoptotic/filtered/"
velocyto_path = "/home/alicen/Projects/limiric/test_data/ground_truth/apoptotic/velocyto/"
apoptotic <- seurat

project_name = "HEK293_control"
organism = "Hsap"
raw_path = "/home/alicen/Projects/limiric/test_data/ground_truth/healthy/raw/"
filtered_path = "/home/alicen/Projects/limiric/test_data/ground_truth/healthy/filtered/"
velocyto_path = "/home/alicen/Projects/limiric/test_data/ground_truth/healthy/velocyto/"
HEK293_control <- seurat

project_name = "pro_apoptotic"
organism = "Hsap"
raw_path = "/home/alicen/Projects/limiric/test_data/ground_truth/proapoptotic/raw/"
filtered_path = "/home/alicen/Projects/limiric/test_data/ground_truth/proapoptotic/filtered/"
velocyto_path = "/home/alicen/Projects/limiric/test_data/ground_truth/proapoptotic/velocyto/"
pro_apoptotic <- seurat

# Getting damaged cells in isolation first 

# Preprocess (SoupX, DoubletFinder, nf scores) 
project_name <- "apoptotic"
organism <- "Hsap"
raw_path <- "/home/alicen/Projects/limiric/test_data/ground_truth/apoptotic/raw/"
apoptotic <- seurat
saveRDS(apoptotic, 
        "/home/alicen/Projects/limiric/damage_left_behind_analysis/groundtruth/R_objects/apoptotic_isolated.rds")


project_name <- "pro_apoptotic"
organism <- "Hsap"
filtered_path <- "/home/alicen/Projects/limiric/test_data/ground_truth/proapoptotic/filtered/"
pro_apoptotic <- seurat
saveRDS(pro_apoptotic, 
"/home/alicen/Projects/limiric/damage_left_behind_analysis/groundtruth/R_objects/pro_apoptotic_isolated.rds")



# Run the Seurat merge function for apoptotic and live control
HEK293_apo_live <- merge(
  x = apoptotic,
  y = HEK293_control,
  add.cell.ids = c("apoptotic", "control"),  
  project = "HEK293"
) # Run through benchmark with 3 cluster ranks for limiric, 1D miQC
ddqc_path = "/home/alicen/Projects/limiric/damage_left_behind_analysis/ddqc/output/apoptotic_ddqc_output.csv"

# Same for pro-apoptotic and control sample
HEK293_proapo_live <- merge(
  x = pro_apoptotic,
  y = HEK293_control,
  add.cell.ids = c("pro_apoptotic", "control"),  
  project = "HEK293"
) # Linear miQC
ddqc_path = "/home/alicen/Projects/limiric/damage_left_behind_analysis/ddqc/output/pro_apoptotic_ddqc_output.csv"


# Other ground truth dataset -----

# SA928 ground truth dataset
dead_SA928 <- 
readRDS("/home/alicen/Projects/limiric/test_data/tumour_groundtruth/batchfx-for-zenodo/TENX049_SA928_003_sceset_v3_raw.rds")
dead_SA928 <- as.Seurat(dead_SA928) 
dying_SA928 <- 
readRDS("/home/alicen/Projects/limiric/test_data/tumour_groundtruth/batchfx-for-zenodo/TENX049_SA928_002_sceset_v3_raw.rds")
dying_SA928 <- as.Seurat(dying_SA928)
live_SA928 <- 
readRDS("/home/alicen/Projects/limiric/test_data/tumour_groundtruth/batchfx-for-zenodo/TENX049_SA928_001_sceset_v3_raw.rds")
live_SA928 <- as.Seurat(live_SA928)


# SA604 ground truth dataset
live_SA604_rep2 <- 
readRDS("/home/alicen/Projects/limiric/test_data/tumour_groundtruth/batchfx-for-zenodo/TENX019_SA604X7XB02089_002_sceset_v3_raw.rds")
live_SA604_rep2 <- as.Seurat(live_SA604_rep2)
dead_SA604 <- 
readRDS("/home/alicen/Projects/limiric/test_data/tumour_groundtruth/batchfx-for-zenodo/TENX019_SA604X7XB02089_004_sceset_v3_raw.rds")
dead_SA604 <- as.Seurat(dead_SA604)


# Above processed files (Zenodo) have non-conventional gene names that need to be changed
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
live_SA604_rep2 <- edit_metadata(live_SA604_rep2)
dead_SA604 <- edit_metadata(dead_SA604)


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
SA928_livendead <- list(live_SA928, dead_SA928)
SA928_livendying <- list(live_SA928, dying_SA928)
SA604 <- list(live_SA604_rep2, dead_SA604)

# Run
SA928_livendead <- merge_and_see(SA928_livendead, mtrb_reduce = FALSE, c("live", "dead"), "SA928_livendead")
SA928_livendying <- merge_and_see(SA928_livendying, mtrb_reduce = FALSE, c("live", "dying"), "SA928_livendying")
SA604 <- merge_and_see(SA604, mtrb_reduce = FALSE, c("live", "dead"), "SA604_livendead")

# check & save
SA928_livendead$plot
DimPlot(SA928_livendead$seurat)
saveRDS(SA928_livendead$seurat, 
"/home/alicen/Projects/limiric/damage_left_behind_analysis/groundtruth/R_objects/SA928_dead_live.rds")

SA928_livendying$plot
saveRDS(SA928_livendying$seurat, 
"/home/alicen/Projects/limiric/damage_left_behind_analysis/groundtruth/R_objects/SA928_dying_live.rds")


SA604$plot
saveRDS(SA604$seurat, "/home/alicen/Projects/limiric/damage_left_behind_analysis/groundtruth/R_objects/SA604.rds")




