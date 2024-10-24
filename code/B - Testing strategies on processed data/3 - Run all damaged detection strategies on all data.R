# SCRIPT CONTEXT 
#
# Benchmarking damaged cell detection tools. First, input raw files are processed. 
# Then, the ten damaged strategies are applied to each processed count matrix
# to obtain output damaged labels for each cell barcode.
# For all methods, this will either be 'damaged' or 'cell'.
#
# Tools tested : 
# 1. ddqc (loading results)
# 2. ensembleKQC (loading results)
# 3. DropletQC : package functions 
# 4. miQC : select best model
# 5. valiDrops : modified function 
# + scater PCA 
#
# NB: 
# - Output from ddqc and ensembleKQC must be done
#     (2.1 - Run ddqc.ipynb) and (2.2 - Run ensembleKQC.sh)
# - Helper functions must be in global environment 
#     - 3.1 Helper plotting function 
#     - 3.2 Helper summarising function 
#     - 3.3 Helper valiDrops function


#-------------------------------------------------------------------------------
# PREPARATIONS
#-------------------------------------------------------------------------------

# Load libraries -------

packages <- c("cowplot", "devtools", "dplyr", "ggplot2", "glmGamPoi", "Matrix", "robustbase",
              "png", "Seurat", "tidyr", 
              "limiric", "miQC", "SingleCellExperiment",
              "scuttle", "presto", "valiDrops", "DropletQC")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    library(pkg)
  }
}

# Load processed datasets ----

# Cell lines 
A549 <- readRDS("/home/alicen/Projects/ReviewArticle/R_objects/preprocess_output/A549_processed.rds")
HCT116 <- readRDS("/home/alicen/Projects/ReviewArticle/R_objects/preprocess_output/HCT116_processed.rds")
Jurkat <- readRDS("/home/alicen/Projects/ReviewArticle/R_objects/preprocess_output/jurkat_processed.rds")

# Healthy, human tissue extracts
hLiver <- readRDS("/home/alicen/Projects/ReviewArticle/R_objects/preprocess_output/hLiver_processed.rds")
hLung <- readRDS("/home/alicen/Projects/ReviewArticle/R_objects/preprocess_output/hLung_processed.rds")
hPBMC <- readRDS("/home/alicen/Projects/ReviewArticle/R_objects/preprocess_output/hPBMC_processed.rds")

# Diseased, human tissue extracts
dLiver <- readRDS("/home/alicen/Projects/ReviewArticle/R_objects/preprocess_output/dLiver_processed.rds")
dLung <- readRDS("/home/alicen/Projects/ReviewArticle/R_objects/preprocess_output/dLung_processed.rds")
dPBMC <- readRDS("/home/alicen/Projects/ReviewArticle/R_objects/preprocess_output/dPBMC_processed.rds")

# Mouse tissue extracts 
mLiver <- readRDS("/home/alicen/Projects/ReviewArticle/R_objects/preprocess_output/mLiver_processed.rds")
mLung <- readRDS("/home/alicen/Projects/ReviewArticle/R_objects/preprocess_output/mLung_processed.rds")
mPBMC <- readRDS("/home/alicen/Projects/ReviewArticle/R_objects/preprocess_output/mPBMC_processed.rds")

# Tumour isolates 
ductal <- readRDS("/home/alicen/Projects/ReviewArticle/R_objects/preprocess_output/ductal_processed.rds")
glio <- readRDS("/home/alicen/Projects/ReviewArticle/R_objects/preprocess_output/glio_processed.rds")
hodgkin <- readRDS("/home/alicen/Projects/ReviewArticle/R_objects/preprocess_output/hodgkin_processed.rds")

# Ground truth 
GM18507_dead <- readRDS("/home/alicen/Projects/ReviewArticle/R_objects/preprocess_groundtruth/GM18507_dead_control.rds")
GM18507_dying <- readRDS("/home/alicen/Projects/ReviewArticle/R_objects/preprocess_groundtruth/GM18507_dying_control.rds")
HEK293_apo <- readRDS("/home/alicen/Projects/ReviewArticle/R_objects/preprocess_groundtruth/HEK293_apo_control.rds")
HEK293_pro <- readRDS("/home/alicen/Projects/ReviewArticle/R_objects/preprocess_groundtruth/HEK293_pro_control.rds")
PDX_dead <- readRDS("/home/alicen/Projects/ReviewArticle/R_objects/preprocess_groundtruth/PDX_dead_control.rds")




#-------------------------------------------------------------------------------
# FUNCTION DEFINED 
#-------------------------------------------------------------------------------


# Function for main benchmarking ------

benchmark <- function(
    seurat, 
    project_name,        # string with dataset identifier 
    organism = "Hsap",   # Either Hsap or Mmus
    model_method = NULL, # Adjustable, alternative to specify the model type "linear", "spline", "polynomial", or "one_dimensional"
    ddqc_path,           # path to output ddqc labels 
    ensembleKQC_path,    # path to output ensembleKQC labels 
    output_path          # where all outputs are saved 
    
){
  
  message("\nBegin testing ", project_name, "...")
  
  # Tool 1: ddqc -------
  
  message("Tool 1: ddqc ...")
  
  # Read in ddqc results 
  ddqc_results <- read.csv2(ddqc_path, sep = ',')
  ddqc_results <- ddqc_results[, c("barcodekey", "passed_qc")]
  
  # Add to Seurat 
  seurat$ddqc <- ddqc_results$passed_qc[match(ddqc_results$barcodekey, rownames(seurat@meta.data))]
  seurat$ddqc <- ifelse(seurat$ddqc == "False", "damaged", "cell")
  seurat$ddqc <- ifelse(is.na(seurat$ddqc), "cell", seurat$ddqc)
  
  
  # Tool 2: ensembleKQC -------
  
  message("Tool 2: EnsembleKQC ...")
  
  # Read in output and add to Seurat meta data 
  ensembleKQC_results <- read.csv2(ensembleKQC_path, sep = ',', col.names = "EnsembleKQC")
  seurat$ensembleKQC <- ensembleKQC_results$EnsembleKQC[match(rownames(ensembleKQC_results), rownames(seurat@meta.data))]
  seurat$ensembleKQC <- ifelse(is.na(seurat$ensembleKQC), "cell", seurat$ensembleKQC)
  
  # Tool 3: DropletQC ------
  
  message("Tool 3: DropletQC ...")
  
  if ("nf" %in% colnames(seurat@meta.data)){
    
    # Extract nf meta data & associated cell barcode from Seurat object
    edDf <- data.frame(nf = as.numeric(seurat$nf), umi = seurat$nCount_RNA)
    nf_col <- as.numeric(seurat$nf)
    
  } else {
    
    # Extract nf alternative nf_malat1 if nf absent & associated cell barcode from Seurat object
    edDf <- data.frame(nf = as.numeric(seurat$nf_malat1), umi = seurat$nCount_RNA)
    nf_col <- as.numeric(seurat$nf_malat1)
    
  }
  
  # Use droplet_qc function to identify empty droplets
  edresultsDf <- identify_empty_drops(edDf)
  
  # Identify damaged cells

  # Generate rough cluster labels 
  temp <- NormalizeData(seurat, verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    FindNeighbors(dims = 1:10, verbose = FALSE) %>%
    FindClusters(verbose = FALSE, resolution = 0.1)
  
  # Extract cluster assignments 
  cell_type <- temp$seurat_clusters 
  
  # Create input data frame 
  dcDf <- data.frame(
    nf = nf_col,
    umi = seurat$nCount_RNA,
    cell_status = edresultsDf$cell_status,
    cell_type = cell_type
  )
  
  # Run damaged detection & ignore empty droplets 
  dcresultsDf <- identify_damaged_cells(nf_umi_ed_ct = dcDf)
  seurat$DropletQC <- dcresultsDf$df$cell_status[match(rownames(dcresultsDf$df), rownames(seurat@meta.data))]
  seurat$DropletQC <- ifelse(seurat$DropletQC == "damaged_cell", "damaged", "cell")
  
  
  # Tool 4: miQC -------
  
  # Create input
  counts <- seurat@assays$RNA
  sce <- SingleCellExperiment(list(counts = counts))
  rowData(sce)$featureType <- rownames(seurat@assays$RNA)
  colData(sce)$subsets_mito_percent <- seurat$mt.percent
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
  
  
  message("\nTool 4: miQC ...    ", best_model_type) 
  
  
  # Tool 5: valiDrops -------
  
  message("Tool 5: valiDrops ...")
  
  # Extract matrix from Seurat object
  expression_matrix <- GetAssayData(seurat, layer = "counts")
  
  if (organism == "Hsap") {species <- "human"}
  if (organism == "Mmus") {species <- "mouse"}
  
  # Run valiDrops internal function to get list of quality control metrics
  expression_metrics <- quality_metrics(counts = expression_matrix,
                                        species = species)
  
  
  # Run valiDrops damaged detection (label_dead valiDrops function edited to allow for this input, not previous ValiDrop input)
  valiDrops <- label_dead(counts = expression_matrix, 
                          metrics = expression_metrics$metrics)
  
  # Add to benchmark seurat object
  seurat$valiDrops <- valiDrops$label[match(rownames(seurat@meta.data), valiDrops$metrics$barcode)]
  seurat$valiDrops <- ifelse(seurat$valiDrops == "dead", "damaged", "cell")
  
  
  # Manual methods 
  
  # Extract df and add cell barcode column 
  seurat_df <- seurat@meta.data
  seurat_df$barcodes <- rownames(seurat_df)

  # Transform values for robust outlier detection 
  seurat_df$nCount_RNA <- log10(seurat_df$nCount_RNA)
  seurat_df$nFeature_RNA <- log10(seurat_df$nFeature_RNA)
  
  # Transform malat1 counts for cases where they are close to zero 
  if (project_name %in% c("hLiver")){
    
    seurat_df$malat1.percent <- log10(seurat_df$malat1.percent)
    
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
    
  # Manual 1: All multivariate outliers (robustbase): UMI and feature counts & mitochondrial, ribosomal and MALAT1 percentages -------
  
  # Define input df, run outlier detection, and add to Seurat meta data 
  manual1_df <- seurat_df[, c("nCount_RNA", "nFeature_RNA", "mt.percent", "rb.percent", "malat1.percent")]
  manual1_results <- identify_moutliers(manual1_df)
  seurat$manual_all <- manual1_results$outliers
  
  # Correct for NA if needed
  seurat$manual_all <- ifelse(is.na(seurat$manual_all), "cell", seurat$manual_all)
  
  # Manual 2: Mito/Ribo multivariate outliers (robustbase): UMI and feature counts & mitochondrial and ribosomal percentages -------
  
  manual2_df <- seurat_df[, c("nCount_RNA", "nFeature_RNA", "mt.percent", "rb.percent")]
  manual2_results <- identify_moutliers(manual2_df)
  seurat$manual_mito_ribo <- manual2_results$outliers
 
  
  # Manual 3: Mito multivariate outliers (robustbase): UMI and feature counts & mitochondrial percentage -------
  
  manual3_df <- seurat_df[, c("nCount_RNA", "nFeature_RNA", "mt.percent")]
  manual3_results <- identify_moutliers(manual3_df)
  seurat$manual_mito <- manual3_results$outliers
  
  # Manual 4: MALAT1 multivariate outliers (robustbase): UMI and feature counts & MALAT1 percentage -------
  
  manual4_df <- seurat_df[, c("nCount_RNA", "nFeature_RNA", "malat1.percent")]
  manual4_results <- identify_moutliers(manual4_df)
  seurat$manual_malat1 <- manual4_results$outliers
  
  # Correct for NA if needed
  seurat$manual_malat1 <- ifelse(is.na(seurat$manual_malat1), "cell", seurat$manual_malat1)
  
  
  # Manual 5: Mito univariate outliers (MAD) ------
  
  manual5_df  <- seurat_df[, c("mt.percent")]
  
  # Calculate the median and MAD 
  median_mito <- median(manual5_df, na.rm = TRUE)
  mad_mito <- mad(manual5_df, constant = 1, na.rm = TRUE)  # Use constant = 1 to follow the typical MAD definition
  
  # Define the threshold for outliers as 3 times the MAD 
  threshold_upper <- median_mito + 3 * mad_mito
  outliers_mito <- manual5_df > threshold_upper 
  
  manual5_df  <- as.data.frame(seurat_df[, c("mt.percent")])
  rownames(manual5_df) <- rownames(seurat_df)
  manual5_results <- data.frame(barcodes = rownames(manual5_df), outliers = outliers_mito)
  manual5_results$outliers <- ifelse(manual5_results$outliers == "TRUE", "damaged", "cell")

  # Add to seurat object 
  seurat$manual_mito_isolated <- manual5_results$outliers
  
  
  # Return output -------
  
  # Call helper function 
  
  final_df <- summarise_results(seurat = seurat)

  write.csv(final_df, 
            paste0(output_path, "/", project_name, "_summary.csv"), 
            quote = FALSE, 
            row.names = FALSE) 
  
  write.csv(seurat@meta.data, 
            paste0(output_path, "/", project_name, ".csv"), 
            quote = FALSE, 
            row.names = TRUE) 
  
  
  # Plot the output -------
  
  # Call helper function
  plot <- BenchPlot(seurat, 
                    output_file = paste0(output_path, "/", project_name, ".png"), 
                    organism)
  
  
  # Global environment output 
  
  return(seurat)
  
}

#-------------------------------------------------------------------------------
# FUNCTION RUN 
#-------------------------------------------------------------------------------

# Cell lines 
A549 <- benchmark(seurat = A549,
                  project_name = "A549",
                  model_method = "linear",
                  ddqc_path = "/home/alicen/Projects/ReviewArticle/python/ddqc_output/A549.csv",
                  ensembleKQC_path = "/home/alicen/Projects/ReviewArticle/python/EnsembleKQC_output/A549_input.csv",
                  output_path = "/home/alicen/Projects/ReviewArticle/benchmark_results")

HCT116 <- benchmark(seurat = HCT116,
                  project_name = "HCT116",
                  ddqc_path = "/home/alicen/Projects/ReviewArticle/python/ddqc_output/HCT116.csv",
                  ensembleKQC_path = "/home/alicen/Projects/ReviewArticle/python/EnsembleKQC_output/HCT116_input.csv",
                  output_path = "/home/alicen/Projects/ReviewArticle/benchmark_results")

Jurkat <- benchmark(seurat = Jurkat,
                  project_name = "Jurkat",
                  # model_method = "linear",
                  ddqc_path = "/home/alicen/Projects/ReviewArticle/python/ddqc_output/jurkat.csv",
                  ensembleKQC_path = "/home/alicen/Projects/ReviewArticle/python/EnsembleKQC_output/jurkat_input.csv",
                  output_path = "/home/alicen/Projects/ReviewArticle/benchmark_results")


# Diseased 
dLiver <- benchmark(seurat = dLiver,
                    project_name = "dLiver",
                    # model_method = "linear",
                    ddqc_path = "/home/alicen/Projects/ReviewArticle/python/ddqc_output/dLiver.csv",
                    ensembleKQC_path = "/home/alicen/Projects/ReviewArticle/python/EnsembleKQC_output/dLiver_input.csv",
                    output_path = "/home/alicen/Projects/ReviewArticle/benchmark_results")


dLung <- benchmark(seurat = dLung,
                    project_name = "dLung",
                    # model_method = "linear",
                    ddqc_path = "/home/alicen/Projects/ReviewArticle/python/ddqc_output/dLung.csv",
                    ensembleKQC_path = "/home/alicen/Projects/ReviewArticle/python/EnsembleKQC_output/dLung_input.csv",
                    output_path = "/home/alicen/Projects/ReviewArticle/benchmark_results")

dPBMC <- benchmark(seurat = dPBMC,
                    project_name = "dPBMC",
                    # model_method = "linear",
                    ddqc_path = "/home/alicen/Projects/ReviewArticle/python/ddqc_output/dPBMC.csv",
                    ensembleKQC_path = "/home/alicen/Projects/ReviewArticle/python/EnsembleKQC_output/dPBMC_input.csv",
                    output_path = "/home/alicen/Projects/ReviewArticle/benchmark_results")

# Healthy 
hLiver <- benchmark(seurat = hLiver,
                    project_name = "hLiver",
                    model_method = "one-dimensional",
                    ddqc_path = "/home/alicen/Projects/ReviewArticle/python/ddqc_output/hLiver.csv",
                    ensembleKQC_path = "/home/alicen/Projects/ReviewArticle/python/EnsembleKQC_output/hLiver_input.csv",
                    output_path = "/home/alicen/Projects/ReviewArticle/benchmark_results")


hLung <- benchmark(seurat = hLung,
                   project_name = "hLung",
                   # model_method = "linear",
                   ddqc_path = "/home/alicen/Projects/ReviewArticle/python/ddqc_output/hLung.csv",
                   ensembleKQC_path = "/home/alicen/Projects/ReviewArticle/python/EnsembleKQC_output/hLung_input.csv",
                   output_path = "/home/alicen/Projects/ReviewArticle/benchmark_results")

hPBMC <- benchmark(seurat = hPBMC,
                   project_name = "hPBMC",
                   # model_method = "linear",
                   ddqc_path = "/home/alicen/Projects/ReviewArticle/python/ddqc_output/hPBMC.csv",
                   ensembleKQC_path = "/home/alicen/Projects/ReviewArticle/python/EnsembleKQC_output/hPBMC_input.csv",
                   output_path = "/home/alicen/Projects/ReviewArticle/benchmark_results")


# Mouse 
mLiver <- benchmark(seurat = mLiver,
                    project_name = "mLiver",
                    # model_method = "linear",
                    organism = "Mmus",
                    ddqc_path = "/home/alicen/Projects/ReviewArticle/python/ddqc_output/mLiver.csv",
                    ensembleKQC_path = "/home/alicen/Projects/ReviewArticle/python/EnsembleKQC_output/mLiver.csv",
                    output_path = "/home/alicen/Projects/ReviewArticle/benchmark_results")


mLung <- benchmark(seurat = mLung,
                   project_name = "mLung",
                   # model_method = "linear",
                   organism = "Mmus",
                   ddqc_path = "/home/alicen/Projects/ReviewArticle/python/ddqc_output/mLung.csv",
                   ensembleKQC_path = "/home/alicen/Projects/ReviewArticle/python/EnsembleKQC_output/mLung.csv",
                   output_path = "/home/alicen/Projects/ReviewArticle/benchmark_results")

mPBMC <- benchmark(seurat = mPBMC,
                   project_name = " PBMC",
                   # model_method = "linear",
                   organism = "Mmus",
                   ddqc_path = "/home/alicen/Projects/ReviewArticle/python/ddqc_output/mPBMC.csv",
                   ensembleKQC_path = "/home/alicen/Projects/ReviewArticle/python/EnsembleKQC_output/mPBMC.csv",
                   output_path = "/home/alicen/Projects/ReviewArticle/benchmark_results")


# Tumour isolates 
ductal <- benchmark(seurat = ductal,
                   project_name = "ductal",
                   model_method = "linear",
                   ddqc_path = "/home/alicen/Projects/ReviewArticle/python/ddqc_output/ductal.csv",
                   ensembleKQC_path = "/home/alicen/Projects/ReviewArticle/python/EnsembleKQC_output/ductal_input.csv",
                   output_path = "/home/alicen/Projects/ReviewArticle/benchmark_results")

glio <- benchmark(seurat = glio,
                    project_name = "glio",
                    # model_method = "linear",
                    ddqc_path = "/home/alicen/Projects/ReviewArticle/python/ddqc_output/glio.csv",
                    ensembleKQC_path = "/home/alicen/Projects/ReviewArticle/python/EnsembleKQC_output/glio_input.csv",
                    output_path = "/home/alicen/Projects/ReviewArticle/benchmark_results")

hodgkin <- benchmark(seurat = hodgkin,
                    project_name = "hodgkin",
                    model_method = "linear",
                    ddqc_path = "/home/alicen/Projects/ReviewArticle/python/ddqc_output/hodgkin.csv",
                    ensembleKQC_path = "/home/alicen/Projects/ReviewArticle/python/EnsembleKQC_output/hodgkin_input.csv",
                    output_path = "/home/alicen/Projects/ReviewArticle/benchmark_results")


# Ground truth 
GM18507_dead <- benchmark(seurat = GM18507_dead, 
                          project_name = "GM18507_dead",
                          ddqc_path = "/home/alicen/Projects/ReviewArticle/python/ddqc_output/GM18507_dead.csv", 
                          ensembleKQC_path = "/home/alicen/Projects/ReviewArticle/python/EnsembleKQC_output/GM18507dead.csv",
                          output_path = "/home/alicen/Projects/ReviewArticle/benchmark_results")
  
GM18507_dying <- benchmark(seurat = GM18507_dying, 
                          project_name = "GM18507_dying",
                          ddqc_path = "/home/alicen/Projects/ReviewArticle/python/ddqc_output/GM18507_dying.csv", 
                          ensembleKQC_path = "/home/alicen/Projects/ReviewArticle/python/EnsembleKQC_output/GM18507dying.csv",
                          output_path = "/home/alicen/Projects/ReviewArticle/benchmark_results")


HEK293_apo <-  benchmark(seurat = HEK293_apo , 
                         project_name = "HEK293_apo",
                         model_method = "polynomial",
                         ddqc_path = "/home/alicen/Projects/ReviewArticle/python/ddqc_output/HEK293_apo.csv", 
                         ensembleKQC_path = "/home/alicen/Projects/ReviewArticle/python/EnsembleKQC_output/HEK293apo.csv",
                         output_path = "/home/alicen/Projects/ReviewArticle/benchmark_results")

  
HEK293_pro <-  benchmark(seurat = HEK293_pro, 
                         project_name = "HEK293_pro",
                         model_method = "one-dimensional",
                         ddqc_path = "/home/alicen/Projects/ReviewArticle/python/ddqc_output/HEK293_pro.csv", 
                         ensembleKQC_path = "/home/alicen/Projects/ReviewArticle/python/EnsembleKQC_output/HEK293pro.csv",
                         output_path = "/home/alicen/Projects/ReviewArticle/benchmark_results")


PDX_dead <- benchmark(seurat = PDX_dead, 
                      project_name = "PDX_dead",
                      model_method = "linear",
                      ddqc_path = "/home/alicen/Projects/ReviewArticle/python/ddqc_output/PDX_dead.csv", 
                      ensembleKQC_path = "/home/alicen/Projects/ReviewArticle/python/EnsembleKQC_output/PDX.csv",
                      output_path = "/home/alicen/Projects/ReviewArticle/benchmark_results")


### End 