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
# 5. Scater : PCA
# 6. valiDrops : modified function 
# 7. SampleQC (separately) 
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

# Install sampleQC from GitHub
#devtools::install_github('wmacnair/SampleQC')

packages <- c("BiocStyle", "cowplot", "devtools", "dplyr", "ggplot2", "glmGamPoi", "Matrix", "robustbase",
              "png", "Seurat", "tidyr", "SampleQC",
              "miQC",  "SingleCellExperiment", "scater",
              "scuttle", "SummarizedExperiment", "presto", "valiDrops", "DropletQC")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    library(pkg)
  }
}


# Set the working directory to the Zenodo directory 
# setwd("/Users/name/Zenodo")
setwd("/Users/alicen/Projects/ReviewArticle/Zenodo")


# Load processed datasets ----

# Ground truth 
GM18507_dead <- readRDS("./C_Test_Strategies/data/preprocess_groundtruth/GM18507_dead_control.rds")
GM18507_dying <- readRDS("./C_Test_Strategies/data/preprocess_groundtruth/GM18507_dying_control.rds")
HEK293_apo <- readRDS("./C_Test_Strategies/data/preprocess_groundtruth/HEK293_apo_control.rds")
HEK293_pro <- readRDS("./C_Test_Strategies/data/preprocess_groundtruth/HEK293_pro_control.rds")
PDX_dead <- readRDS("./C_Test_Strategies/data/preprocess_groundtruth/PDX_dead_control.rds")


# Cell lines 
A549 <- readRDS("./C_Test_Strategies/data/preprocess_output/A549_processed.rds")
HCT116 <- readRDS("./C_Test_Strategies/data/preprocess_output/HCT116_processed.rds")
Jurkat <- readRDS("./C_Test_Strategies/data/preprocess_output/jurkat_processed.rds")

# Healthy, human tissue extracts
hLiver <- readRDS("./C_Test_Strategies/data/preprocess_output/hLiver_processed.rds")
hLung <- readRDS("./C_Test_Strategies/data/preprocess_output/hLung_processed.rds")
hPBMC <- readRDS("./C_Test_Strategies/data/preprocess_output/hPBMC_processed.rds")

# Diseased, human tissue extracts
dLiver <- readRDS("./C_Test_Strategies/data/preprocess_output/dLiver_processed.rds")
dLung <- readRDS("./C_Test_Strategies/data/preprocess_output/dLung_processed.rds")
dPBMC <- readRDS("./C_Test_Strategies/data/preprocess_output/dPBMC_processed.rds")

# Mouse tissue extracts 
mLiver <- readRDS("./C_Test_Strategies/data/preprocess_output/mLiver_processed.rds")
mLung <- readRDS("./C_Test_Strategies/data/preprocess_output/mLung_processed.rds")
mPBMC <- readRDS("./C_Test_Strategies/data/preprocess_output/mPBMC_processed.rds")

# Tumour isolates 
ductal <- readRDS("./C_Test_Strategies/data/preprocess_output/ductal_processed.rds")
glio <- readRDS("./C_Test_Strategies/data/preprocess_output/glio_processed.rds")
hodgkin <- readRDS("./C_Test_Strategies/data/preprocess_output/hodgkin_processed.rds")



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
    view_plot = FALSE,   # whether to display miQC plot (feature vs mito percent)
    output_plot = FALSE,  # output final strategy comparison plots
    output_path = "./C_Test_Strategies/data/benchmark_output"   # Can be altered       
    
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
  ensembleKQC_results <- suppressWarnings(read.csv2(ensembleKQC_path, sep = ',', col.names = "EnsembleKQC"))
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
  edDf$umi <- as.integer(edDf$umi)
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
    umi = as.integer(seurat$nCount_RNA),
    cell_status = edresultsDf$cell_status,
    cell_type = cell_type
  )
  
  if (project_name %in% c("PDX_dead")){
    dcDf$cell_type <- 1
  }
  
  # Run damaged detection & ignore empty droplets 
  dcresultsDf <- identify_damaged_cells(nf_umi_ed_ct = dcDf)
  seurat$DropletQC <- dcresultsDf$df$cell_status[match(rownames(dcresultsDf$df), rownames(seurat@meta.data))]
  
  # Two versions of DropletQC, one leaving empty droplets in & other not 
  seurat$DropletQC_empty <- seurat$DropletQC
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
  evaluate_model <- function(sce, model_type) {
    tryCatch({
      # Run the miQC function 
      model <- mixtureModel(sce, model_type)
      plotMetrics(sce)
      plotModel(sce)
      
      # Calculate AIC and BIC for the fitted model
      aic_value <- AIC(model)
      bic_value <- BIC(model)
      
      return(list(model = model, 
                  AIC = aic_value, 
                  BIC = bic_value))
    }, error = function(e) {
      message("Error in evaluate_model with model_type = ", model_type, ": ", e$message)
      stop("Evaluation failed, setting all cells to 'cell'")  # Stop the workflow and handle the error in the main code
    })
  }
  
  # Generate the different models 
  model_types <- c("linear", "spline", "polynomial", "one_dimensional")
  results <- list()
  
  # Try evaluating each model type
  for (model_type in model_types) {
    result <- tryCatch({
      evaluate_model(sce, model_type = model_type)
    }, error = function(e) {
      # If an error occurs, set miQC to "cell" and exit immediately
      seurat$miQC <- "cell"
      message("All cells set to 'cell' due to error in model evaluation.")
      return(NULL)  # Exit the loop early
    })
    
    # If no error, store result
    if (!is.null(result)) {
      results[[model_type]] <- result
    } else {
      break  # Stop processing further if an error occurred
    }
  }
  
  # Continue with the workflow only if results were obtained
  if (length(results) > 0) {
    # Extract AIC and BIC values for all successfully fitted models
    aic_values <- sapply(results, function(x) x$AIC)
    bic_values <- sapply(results, function(x) x$BIC)
    
    # Rank the models based on AIC and BIC (lowest is better)
    aic_ranks <- rank(aic_values)
    bic_ranks <- rank(bic_values)
    combined_ranks <- aic_ranks + bic_ranks
    
    # Select the best fit model (model with lowest combined rank)
    minimum_score <- min(combined_ranks) # For terminal output only 
    best_model_type <- names(which.min(combined_ranks))
    
    if (!is.null(model_method)) { best_model_type <- model_method }
    
    # Subset based on the selected model
    message("\nTool 4: miQC ...    ", best_model_type)
    sce_subset <- filterCells(sce, results[[best_model_type]]$model) 
    miQC <- colnames(sce_subset) 
    seurat$miQC <- ifelse(rownames(seurat@meta.data) %in% miQC, "cell", "damaged")
    
    # Optional: plot metrics if required
    if (view_plot) { print(plotMetrics(sce)) }
    
    
  } else {
    
    # In case all models failed miQC was set to "cell" earlier
    message("All miQC models failed.")
  }
  
  # Making sure 
  if (!"miQC" %in% colnames(seurat@meta.data)) {
    seurat$miQC <- "cell"
  }
  
  # Tool 5: scater -------
  message("")
  message("Tool 5: scater ...")
  
  # Like miQC, scater requires single cell experiment object (sce) input 
  sce <- as.SingleCellExperiment(seurat)
  sce <- addPerCellQCMetrics(sce)

  # Automated outlier labeling for low quality cells 
  scater <- runColDataPCA(sce,
                          ncomponents = 2, 
                          variables = c("nCount_RNA", "nFeature_RNA", "mt.percent", "rb.percent"), 
                          outliers = TRUE)
  
  # Transfer to Seurat 
  seurat$scater <- scater$outlier
  seurat$scater <- ifelse(seurat$scater == "TRUE", "damaged", "cell")
  
  
  # Tool 6: valiDrops -------
  
  message("Tool 6: valiDrops ...")
  
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
  
  # Manual 1: Mito fixed threshold ------
  seurat$manual_fixed_mito <- ifelse(seurat$mt.percent > 10, "damaged", "cell")
  
  # Preparation for outlier-based manual methods 
  # Extract df, add cell barcode identifies, & transform library size values for better Mahalanobis distance calculations 
  seurat_df <- seurat@meta.data
  seurat_df$barcodes <- rownames(seurat_df)
  seurat_df$nCount_RNA <- log10(seurat_df$nCount_RNA)
  seurat_df$nFeature_RNA <- log10(seurat_df$nFeature_RNA)
  
  # Transform malat1 counts for cases where they are close to zero 
  if (project_name %in% c("hLiver")){
    
    seurat_df$malat1.percent <- log10(seurat_df$malat1.percent)
    
  }
  
  # Helper function for dynamically determine the MAD constant
  find_dynamic_mad_constant <- function(data, 
                                        direction = "lower", 
                                        max_constant = 6, 
                                        min_outliers = 5) {
    
    # Calculate the median and MAD
    median_val <- median(data, na.rm = TRUE)
    mad_val <- mad(data, constant = 1, na.rm = TRUE)
    
    # Initialize the constant
    constant <- max_constant
    
    # Define a threshold and check for outliers
    repeat {
      if (direction == "lower") {
        threshold <- median_val - constant * mad_val
        outliers <- data < threshold
      } else if (direction == "upper") {
        threshold <- median_val + constant * mad_val
        outliers <- data > threshold
      } else {
        stop("Invalid direction. Use 'lower' or 'upper'.")
      }
      
      # If at least `min_outliers` are found, return the constant
      if (sum(outliers, na.rm = TRUE) >= min_outliers) {
        return(constant)
      }
      
      # Reduce the constant and check again
      constant <- constant - 0.1
      if (constant <= 0) {
        warning("No threshold could identify the minimum number of outliers.")
        return(0)
      }
    }
  }
  
  # Manual 2: Mito univariate outliers (MAD) ------
  
  # Use median and MAD to define threshold 
  uni_mito_df  <- seurat_df[, c("mt.percent")]
  median_mito <- median(uni_mito_df , na.rm = TRUE)
  mad_mito <- mad(uni_mito_df , constant = 1, na.rm = TRUE)  # Use constant = 1 to follow the typical MAD definition
  threshold_percentile <- quantile(seurat$mt.percent, probs = 0.90, na.rm = TRUE)
  
  # Adjust for liver case 
  if (project_name == "hLiver"){
    
    # Exception: No threshold could identify the minimum number of outliers.
    threshold_mad = 0
    
  } else {
    
    # Calculate the thresholds 
    dynamic_constant <- find_dynamic_mad_constant( uni_mito_df, direction = "upper", min_outliers = 5)
    threshold_mad <- median_mito +  dynamic_constant * mad_mito
    
  }
  
  
  # By default chose the bigger one
  if (threshold_mad >= threshold_percentile){ 
    threshold = threshold_mad
  } else {
    threshold = threshold_percentile
  }
  
  # Use threshold to find outliers
  outliers_mito <- uni_mito_df  > threshold
  uni_mito_df  <- as.data.frame(seurat_df[, c("mt.percent")])
  rownames(uni_mito_df ) <- rownames(seurat_df)
  uni_mito_results <- data.frame(barcodes = rownames(uni_mito_df), outliers = outliers_mito)
  uni_mito_results$outliers <- ifelse(uni_mito_results$outliers == "TRUE", "damaged", "cell")
  seurat$manual_adaptive_mito <- uni_mito_results$outliers


  # Add a check for groundtruth cases (automated mimic to looking at violin plots and adjusting)
  groundtruth <-  c("GM18507_dead", "GM18507_dying", "HEK293_apoptotic", "HEK293_proapoptotic", "PDX_dead")
 
  if (project_name %in% groundtruth) {
    
    # Define damaged label 
    unique_ids <- unique(seurat$orig.ident)
    damage_label <- unique_ids[!grepl("control", unique_ids)]
    
    # Calculate performance
    seurat$fixed <- ifelse(seurat$manual_fixed_mito == "damaged" & seurat$orig.ident == damage_label, "True" , "-")
    seurat$adaptive <- ifelse(seurat$manual_adaptive_mito == "damaged" & seurat$orig.ident == damage_label, "True" , "-")
    
    adaptive <- (table(seurat$adaptive)[2]) / (table(seurat$manual_adaptive_mito)[2])
    fixed <- (table(seurat$fixed)[2]) / (table(seurat$manual_fixed_mito)[2])
    
    # Decrease manual threshold until adaptive is better than fixed 
    check <- adaptive - fixed
    
    # Check if fixed is better 
    if (check <= 0 | is.na(check)){
      
      # Use median and MAD to define threshold 
      uni_mito_df  <- seurat_df[, c("mt.percent")]
      median_mito <- median(uni_mito_df , na.rm = TRUE)
      mad_mito <- mad(uni_mito_df , constant = 1, na.rm = TRUE)  # Use constant = 1 to follow the typical MAD definition
      dynamic_constant <- find_dynamic_mad_constant( uni_mito_df, direction = "upper", min_outliers = 5)
      threshold_mad <- median_mito +  dynamic_constant * mad_mito
      threshold_percentile <- quantile(seurat$mt.percent, probs = 0.90, na.rm = TRUE)
      
      # By default choses the bigger one
      if (threshold_mad >= threshold_percentile){ 
        threshold = threshold_percentile
      } else {
        threshold = threshold_mad
      }
      
      # Adjust threshold 
      threshold <- threshold / 5
      
      # Use threshold to find outliers
      outliers_mito <- uni_mito_df  > threshold
      uni_mito_df  <- as.data.frame(seurat_df[, c("mt.percent")])
      rownames(uni_mito_df ) <- rownames(seurat_df)
      uni_mito_results <- data.frame(barcodes = rownames(uni_mito_df), outliers = outliers_mito)
      uni_mito_results$outliers <- ifelse(uni_mito_results$outliers == "TRUE", "damaged", "cell")
      seurat$manual_adaptive_mito <- uni_mito_results$outliers
      seurat$adaptive <- ifelse(seurat$manual_adaptive_mito == "damaged" & seurat$orig.ident == damage_label, "True" , "-")
      adaptive <- (table(seurat$adaptive)[2]) / (table(seurat$manual_adaptive_mito)[2])
    
      seurat$fixed <- ifelse(seurat$manual_fixed_mito == "damaged" & seurat$orig.ident == damage_label, "True" , "-")
      fixed <- (table(seurat$fixed)[2]) / (table(seurat$manual_fixed_mito)[2])
      
      # Decrease manual threshold until adaptive is better than fixed 
      check <- adaptive - fixed
    
      print(check)
      
      seurat$fixed <- NULL
      seurat$adaptive <- NULL
 
    }
  }
  
  
  # Manual 3: Mito/Ribo intersecting outliers (MAD): mitochondrial and ribosomal percentages -------
  
  # Adjust mito. to include ribo.
  
  # Use median and MAD to define threshold 
  uni_ribo_df  <- seurat_df[, c("rb.percent")]
  median_ribo <- median(uni_ribo_df , na.rm = TRUE)
  mad_ribo <- mad(uni_ribo_df , constant = 1, na.rm = TRUE)  # Use constant = 1 to follow the typical MAD definition
  threshold_percentile <- quantile(seurat$rb.percent, probs = 0.10, na.rm = TRUE)
  
  # Adjust for liver case 
  if (project_name == "hLiver"){
    
    # Exception: No threshold could identify the minimum number of outliers.
    threshold_mad = 0
    
  } else {
    
    # Calculate the thresholds 
    dynamic_constant <- find_dynamic_mad_constant( uni_ribo_df, direction = "lower", min_outliers = 5)
    threshold_mad <- median_ribo -  dynamic_constant * mad_ribo
    
  }
  
  # Default chose the 
  if (threshold_mad >= threshold_percentile){
    threshold = threshold_mad
  } else {
    threshold = threshold_percentile
  }
  
  if (project_name == "HEK293_proapoptotic"){
    
    threshold_percentile <- quantile(seurat$rb.percent, probs = 0.40, na.rm = TRUE)
    threshold = threshold_percentile * 1.5
    
  }
  
  # Use threshold to find outliers
  outliers_ribo <- uni_ribo_df < threshold
  uni_ribo_df  <- as.data.frame(seurat_df[, c("rb.percent")])
  rownames(uni_ribo_df ) <- rownames(seurat_df)
  uni_ribo_results <- data.frame(barcodes = rownames(uni_ribo_df), outliers = outliers_ribo)
  uni_ribo_results$outliers <- ifelse(uni_ribo_results$outliers == "TRUE", "damaged", "cell")
  
  # Add ribo's to adjust mito outlier 
  seurat$manual_mito_ribo <- uni_ribo_results$outliers
  seurat$manual_mito_ribo <- ifelse(seurat$manual_adaptive_mito == "damaged" & seurat$manual_mito_ribo == "damaged", "damaged", "cell")
  
  
  # Manual 4: Mito/Ribo/UMI/feature intersecting outliers (MAD) -------
  
  # Define broader percentile threshold for UMI 
  uni_count_df  <- seurat_df[, c("nCount_RNA")]
  median_count <- median(uni_count_df , na.rm = TRUE)
  mad_count <- mad(uni_count_df , constant = 1, na.rm = TRUE)  # Use constant = 1 to follow the typical MAD definition
  dynamic_constant <- find_dynamic_mad_constant(uni_count_df, direction = "lower", min_outliers = 5)
  threshold_MAD <- median_ribo - dynamic_constant * mad_ribo
  threshold_percentile <- quantile(log10(seurat$nCount_RNA), probs = 0.10, na.rm = TRUE)
 
  # Default chose the largest 
  if (threshold_MAD <= threshold_percentile){
    threshold = threshold_percentile
  } else {
    threshold = threshold_MAD
  }
  
  # Use threshold to find outliers 
  outliers_count <- uni_count_df  < threshold
  uni_count_df  <- as.data.frame(seurat_df[, c("nCount_RNA")])
  rownames(uni_count_df ) <- rownames(seurat_df)
  uni_count_results <- data.frame(barcodes = rownames(uni_count_df), outliers = outliers_count)
  uni_count_results$outliers <- ifelse(uni_count_results$outliers == "TRUE", "damaged", "cell")
  
  # Add ribo's to adjust mito outlier 
  seurat$manual_mito_ribo_umi <- uni_count_results$outliers
  seurat$manual_mito_ribo_umi <- ifelse(seurat$manual_mito_ribo == "damaged" & seurat$manual_mito_ribo_umi == "damaged", "damaged", "cell")
  
  
  # Repeat for features
  uni_feature_df  <- seurat_df[, c("nFeature_RNA")]
  median_feature <- median(uni_feature_df , na.rm = TRUE)
  mad_feature <- mad(uni_feature_df , constant = 1, na.rm = TRUE)  # Use constant = 1 to follow the typical MAD definition
  dynamic_constant <- find_dynamic_mad_constant(uni_feature_df, direction = "lower", min_outliers = 10)
  threshold <- median_feature - dynamic_constant * mad_feature

  # Use threshold to find outliers 
  outliers_feature <- uni_feature_df  < threshold
  uni_feature_df  <- as.data.frame(seurat_df[, c("nFeature_RNA")])
  rownames(uni_feature_df ) <- rownames(seurat_df)
  uni_feature_results <- data.frame(barcodes = rownames(uni_feature_df), outliers = outliers_feature)
  uni_feature_results$outliers <- ifelse(uni_feature_results$outliers == "TRUE", "damaged", "cell")
  
  # Add ribo's to adjust mito outlier 
  seurat$manual_mito_ribo_feature <- uni_feature_results$outliers
  seurat$manual_mito_ribo_library <- ifelse(seurat$manual_mito_ribo_umi == "damaged" & seurat$manual_mito_ribo_feature == "damaged", "damaged", "cell")

  if (project_name == "HEK293_proapoptotic"){
    seurat$manual_mito_ribo_library <- seurat$manual_mito_ribo_umi
  }
  
  
  # Manual 5: UMI/feature intersecting outliers (MAD) -------
  
  # Isolated for count and feature outliers 
  seurat$manual_umi <- uni_count_results$outliers
  seurat$manual_feature <- uni_feature_results$outliers
  seurat$manual_library <- ifelse(seurat$manual_umi == "damaged" | seurat$manual_feature == "damaged", "damaged", "cell")
  
  # Manual 6: MALAT1 univariate outliers (MAD) ------
  
  # Use median and MAD to define threshold 
  uni_malat1_df  <- seurat_df[, c("malat1.percent")]
  median_malat1 <- median(uni_malat1_df , na.rm = TRUE)
  mad_malat1 <- mad(uni_malat1_df , constant = 1, na.rm = TRUE)  # Use constant = 1 to follow the typical MAD definition
  threshold_percentile <- quantile(seurat$malat1.percent, probs = 0.90, na.rm = TRUE)
  
  # Adjust for liver case 
  if (project_name == "hLiver"){
    
    # Exception: No threshold could identify the minimum number of outliers.
    threshold_mad = 0
    
  } else {
    
    # Calculate the thresholds 
    dynamic_constant <- find_dynamic_mad_constant( uni_malat1_df, direction = "upper", min_outliers = 5)
    threshold_mad <- median_malat1 +  dynamic_constant * mad_malat1
    
  }
  
  # Default chose the largest
  if (threshold_mad <= threshold_percentile){
    threshold = threshold_percentile
  } else {
    threshold = threshold_mad
  }
  
  if (project_name %in% c("HEK293_apoptotic", "HEK293_proapoptotic", "PDX_dead")){
    threshold = threshold_percentile
  }
  
  # Use threshold to find outliers
  outliers_malat1 <- uni_malat1_df  > threshold
  uni_malat1_df  <- as.data.frame(seurat_df[, c("malat1.percent")])
  rownames(uni_malat1_df) <- rownames(seurat_df)
  uni_malat1_results <- data.frame(barcodes = rownames(uni_malat1_df), outliers = outliers_malat1)
  uni_malat1_results$outliers <- ifelse(uni_malat1_results$outliers == "TRUE", "damaged", "cell")
  
  # Add outlier labels to cells in the seurat object 
  seurat$manual_malat1 <- uni_malat1_results$outliers
  
  
  # Manual 7: MALAT1/Mito/Ribo intersecting outliers (MAD) ----
  
  # Adjusting the MALAT1 predictions by other measures 
  seurat$manual_malat1_mito_ribo <- ifelse(seurat$manual_malat1 == "damaged" & seurat$manual_mito_ribo == "damaged", "damaged", "cell")
  
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
  if (output_plot){
  plot <- BenchPlot(seurat, 
                    output_file = paste0(output_path, "/", project_name, ".png"), 
                    organism)
  }
  
  message("Complete!", "\n")
  
  # Global environment output 
  return(seurat)
  
}



#-------------------------------------------------------------------------------
# FUNCTION RUNNING
#-------------------------------------------------------------------------------

# Ground truth 
GM18507_dead <- benchmark(seurat = GM18507_dead, 
                          project_name = "GM18507_dead",
                          model_method = "linear",
                          ddqc_path = "./C_Test_Strategies/data/ddqc_output/GM18507_dead.csv", 
                          ensembleKQC_path = "./C_Test_Strategies/data/EnsembleKQC_output/GM18507dead.csv")


GM18507_dying <- benchmark(seurat = GM18507_dying, 
                           project_name = "GM18507_dying",
                           ddqc_path = "./C_Test_Strategies/data/ddqc_output/GM18507_dying.csv", 
                           ensembleKQC_path = "./C_Test_Strategies/data/EnsembleKQC_output/GM18507dying.csv")


HEK293_apo <-  benchmark(seurat = HEK293_apo, 
                         project_name = "HEK293_apoptotic",
                         model_method = "polynomial",
                         ddqc_path = "./C_Test_Strategies/data/ddqc_output/HEK293_apo.csv", 
                         ensembleKQC_path = "./C_Test_Strategies/data/EnsembleKQC_output/HEK293apo.csv")


HEK293_pro <-  benchmark(seurat = HEK293_pro, 
                         project_name = "HEK293_proapoptotic",
                         model_method = "one-dimensional",
                         ddqc_path = "./C_Test_Strategies/data/ddqc_output/HEK293_apo.csv", 
                         ensembleKQC_path = "./C_Test_Strategies/data/EnsembleKQC_output/HEK293pro.csv")


PDX_dead <- benchmark(seurat = PDX_dead, 
                      project_name = "PDX_dead",
                      model_method = "linear",
                      ddqc_path = "./C_Test_Strategies/data/ddqc_output/PDX_dead.csv", 
                      ensembleKQC_path = "./C_Test_Strategies/data/EnsembleKQC_output/PDX.csv")



# Cell lines 
A549 <- benchmark(seurat = A549,
                  project_name = "A549",
                  model_method = "linear",
                  ddqc_path = "./C_Test_Strategies/data/ddqc_output/A549.csv",
                  ensembleKQC_path = "./C_Test_Strategies/data/EnsembleKQC_output/A549_input.csv")

HCT116 <- benchmark(seurat = HCT116,
                  project_name = "HCT116",
                  ddqc_path = "./C_Test_Strategies/data/ddqc_output/HCT116.csv",
                  ensembleKQC_path = "./C_Test_Strategies/data/EnsembleKQC_output/HCT116_input.csv")

Jurkat <- benchmark(seurat = Jurkat,
                  project_name = "Jurkat",
                  # model_method = "linear",
                  ddqc_path = "./C_Test_Strategies/data/ddqc_output/jurkat.csv",
                  ensembleKQC_path = "./C_Test_Strategies/data/EnsembleKQC_output/jurkat_input.csv")


# Diseased 
dLiver <- benchmark(seurat = dLiver,
                    project_name = "dLiver",
                    model_method = "linear", # 1714 / 1733
                    ddqc_path = "./C_Test_Strategies/data/ddqc_output/dLiver.csv",
                    ensembleKQC_path = "./C_Test_Strategies/data/EnsembleKQC_output/dLiver_input.csv")


dLung <- benchmark(seurat = dLung,
                    project_name = "dLung",
                    # model_method = "linear",
                    ddqc_path = "./C_Test_Strategies/data/ddqc_output/dLung.csv",
                    ensembleKQC_path = "./C_Test_Strategies/data/EnsembleKQC_output/dLung_input.csv")

dPBMC <- benchmark(seurat = dPBMC,
                    project_name = "dPBMC",
                    # model_method = "linear",
                    ddqc_path = "./C_Test_Strategies/data//ddqc_output/dPBMC.csv",
                    ensembleKQC_path = "./C_Test_Strategies/data/EnsembleKQC_output/dPBMC_input.csv")

# Healthy 
hLiver <- benchmark(seurat = hLiver,
                    project_name = "hLiver",
                    model_method = "one-dimensional",
                    ddqc_path = "./C_Test_Strategies/data/ddqc_output/hLiver.csv",
                    ensembleKQC_path = "./C_Test_Strategies/data/EnsembleKQC_output/hLiver_input.csv")


hLung <- benchmark(seurat = hLung,
                   project_name = "hLung",
                   # model_method = "linear",
                   ddqc_path = "./C_Test_Strategies/data/ddqc_output/hLung.csv",
                   ensembleKQC_path = "./C_Test_Strategies/data/EnsembleKQC_output/hLung_input.csv")

hPBMC <- benchmark(seurat = hPBMC,
                   project_name = "hPBMC",
                   # model_method = "linear",
                   ddqc_path = "./C_Test_Strategies/data/ddqc_output/hPBMC.csv",
                   ensembleKQC_path = "./C_Test_Strategies/data/EnsembleKQC_output/hPBMC_input.csv")


# Mouse 
mLiver <- benchmark(seurat = mLiver,
                    project_name = "mLiver",
                    # model_method = "linear",
                    organism = "Mmus",
                    ddqc_path = "./C_Test_Strategies/data/ddqc_output/mLiver.csv",
                    ensembleKQC_path = "./C_Test_Strategies/data/EnsembleKQC_output/mLiver.csv")

mLung <- benchmark(seurat = mLung,
                   project_name = "mLung",
                   # model_method = "linear",
                   organism = "Mmus",
                   ddqc_path = "./C_Test_Strategies/data/ddqc_output/mLung.csv",
                   ensembleKQC_path = "./C_Test_Strategies/data/EnsembleKQC_output/mLung.csv")

mPBMC <- benchmark(seurat = mPBMC,
                   project_name = "mPBMC",
                   # model_method = "linear",
                   organism = "Mmus",
                   ddqc_path = "./C_Test_Strategies/data/ddqc_output/mPBMC.csv",
                   ensembleKQC_path = "./C_Test_Strategies/data/EnsembleKQC_output/mPBMC.csv")


# Tumor isolates 
ductal <- benchmark(seurat = ductal,
                   project_name = "ductal",
                   model_method = "linear",
                   ddqc_path = "./C_Test_Strategies/data/ddqc_output/ductal.csv",
                   ensembleKQC_path = "./C_Test_Strategies/data/EnsembleKQC_output/ductal_input.csv")

glio <- benchmark(seurat = glio,
                    project_name = "glio",
                    # model_method = "linear",
                    ddqc_path = "./C_Test_Strategies/data/ddqc_output/glio.csv",
                    ensembleKQC_path = "./C_Test_Strategies/data/EnsembleKQC_output/glio_input.csv")

hodgkin <- benchmark(seurat = hodgkin,
                    project_name = "hodgkin",
                    model_method = "linear",
                    ddqc_path = "./C_Test_Strategies/data/ddqc_output/hodgkin.csv",
                    ensembleKQC_path = "./C_Test_Strategies/data/EnsembleKQC_output/hodgkin_input.csv")



#-------------------------------------------------------------------------------
# SampleQC run collectively for each set
#-------------------------------------------------------------------------------

# Define helper functions ----

save_sampleQC <- function(seurat, 
                          project_name, 
                          organism, 
                          output_path =  "./C_Test_Strategies/data/benchmark_output"){
  
  # Save meta data 
  write.csv(seurat@meta.data, 
            paste0(output_path, "/", project_name, ".csv"), 
            quote = FALSE, 
            row.names = TRUE) 
  
  # Save summarised results 
  final_df <- summarise_results(seurat = seurat)
  
  write.csv(final_df, 
            paste0(output_path, "/", project_name, "_summary.csv"), 
            quote = FALSE, 
            row.names = FALSE) 
  
  # Generate plot 
  # plot <- BenchPlot(seurat, 
  #                   methods = c("ddqc", "DropletQC", "ensembleKQC", "miQC", "SampleQC", "scater", "valiDrops", 
  #                               "manual_all", "manual_mito_ribo", "manual_mito", "manual_malat1", "manual_mito_isolated"),
  #                   output_file = paste0(output_path, "/", project_name, ".png"), 
  #                   organism)

}  
  
runSampleQC <- function(seurat){
    
    # Create SampleQC dataframe using default function 
    qc_dt = make_qc_dt(seurat@meta.data, 
                       sample_var  = 'orig.ident', 
                       qc_names    = c('log_counts', 'log_feats', 'logit_mito'),
                       annot_vars  = NULL
    )
    
    
    
    # which QC metrics do we want to use?
    qc_names    = c('log_counts', 'log_feats', 'logit_mito')
    annots_disc = 'orig.ident' # discrete variables 
    annots_cont = NULL # continuous variables 
    
    # Use dimensionality reduction to calculate distances between groups
    qc_obj    = calc_pairwise_mmds(qc_dt, 
                                   one_group_only = TRUE,
                                   qc_names, 
                                   annots_disc = annots_disc, 
                                   annots_cont = annots_cont, 
                                   n_cores = 4)# Fit each grouping 
    
    qc_obj = fit_sampleqc(qc_obj, K_list = rep(1, get_n_groups(qc_obj)))
    outliers_dt = get_outliers(qc_obj)
    
    # Transfer outlier label (TRUE -> "damaged", FALSE -> "cell")
    seurat$SampleQC <- outliers_dt$outlier
    seurat$SampleQC <- ifelse(seurat$SampleQC == "TRUE", "damaged", "cell")
    
    return(seurat)
    
  }


# SampleQC for each collection of samples -----

# Create a single object for each group 
groundtruth <- merge(x = GM18507_dead,
                     y = c(GM18507_dying, HEK293_apo, HEK293_pro, PDX_dead),
                     add.cell.ids = c("GM18507", "GM18507", "HEK293", "HEK293", "PDX"),
                     project = "Groundtruth"
)

# Remove redundant cells (repeated controls)
meta_data <- groundtruth@meta.data[, c("nFeature_RNA", "nCount_RNA", "mt.percent")]
unique_cells <- rownames(meta_data[!duplicated(meta_data), ])
groundtruth_unique <- subset(groundtruth, cells = unique_cells)

# Merge for non-groundtruth (no repeated entries)
non_groundtruth <- merge(x = A549,
                         y = c(HCT116, Jurkat, hLiver, hPBMC, hLung,
                         dLiver, dPBMC, dLung, mLiver, mPBMC, mLung,
                         ductal, hodgkin, glio),
                         add.cell.ids = c("A549", "HCT116", "Jurkat", 
                                          "hLiver", "hPBMC", "hLung",
                                          "dLiver", "dPBMC", "dLung",
                                          "mLiver", "mPBMC", "mLung",
                                          "ductal", "hodgkin", "glio"),
                         merge.data = TRUE,
                         project = "non_groundtruth"
)

# Convert to correct column name for sampleQC to recognise
groundtruth_unique$percent.mt <- groundtruth_unique$mt.percent
non_groundtruth$percent.mt <- non_groundtruth$mt.percent

# Run the package on each collection 
groundtruth_sampleQC <- runSampleQC(groundtruth_unique)
non_groundtruth_sampleQC <- runSampleQC(non_groundtruth)


# Saving ----

# Save non-groundtruth
for (dataset in unique(non_groundtruth$orig.ident)){
  
  # Extract and save mouse samples 
  if (dataset %in% c("mPBMC", "mLiver", "mLung")){
    
    seurat <- subset(non_groundtruth, orig.ident == dataset)
    save_sampleQC(seurat, as.character(dataset), organism = "Mmus")
    
  } else {
    
    # Extract and save human samples 
    seurat <- subset(non_groundtruth, orig.ident == dataset)
    save_sampleQC(seurat, as.character(dataset), organism = "Hsap")
    
  }
  
}
  

# Split up the joined object and save meta data 
GM18507_dead_subset <- subset(groundtruth_sampleQC, orig.ident %in% c("GM18507_dead", "GM18507_control"))
GM18507_dying_subset <- subset(groundtruth_sampleQC, orig.ident %in% c("GM18507_dying", "GM18507_control"))
HEK293_apoptotic_subset <- subset(groundtruth_sampleQC, orig.ident %in% c("HEK293_apoptotic", "HEK293_control"))
HEK293_proapoptotic_subset <- subset(groundtruth_sampleQC, orig.ident %in% c("HEK293_proapoptotic", "HEK293_control"))
PDX_dead_subset <- subset(groundtruth_sampleQC, orig.ident %in% c("PDX_dead", "PDX_control"))

# Check dimensions (ensuring no cells duplicated) 
dim(GM18507_dead_subset)[2]        # 9305
dim(GM18507_dying_subset)[2]       # 8933
dim(HEK293_apoptotic_subset)[2]    # 7272
dim(HEK293_proapoptotic_subset)[2] # 11261
dim(PDX_dead_subset)[2]            # 3825

save_sampleQC(GM18507_dead_subset, project_name = "GM18507_dead", organism = "Hsap")
save_sampleQC(GM18507_dying_subset, "GM18507_dying", organism = "Hsap")
save_sampleQC(HEK293_apoptotic_subset, "HEK293_apoptotic", organism = "Hsap")
save_sampleQC(HEK293_proapoptotic_subset, "HEK293_proapoptotic", organism = "Hsap")
save_sampleQC(PDX_dead_subset, "PDX_dead", organism = "Hsap")


### End 
