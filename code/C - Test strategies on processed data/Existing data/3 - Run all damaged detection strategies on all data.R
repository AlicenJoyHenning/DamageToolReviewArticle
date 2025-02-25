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
              "png", "Seurat", "tidyr", "SampleQC", "ddqcR",
              "miQC",  "SingleCellExperiment", "scater",
              "scuttle", "SummarizedExperiment", "presto", "valiDrops", "DropletQC")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    library(pkg)
  }
}


# Load processed datasets ----
# Read in all .rds files & save as name in list 
directory <- "/Users/alicen/Projects/Damage_analsyis/data"

# Extract filenames without the .rds extension to use as list names
file_paths <- list.files(path = directory, pattern = "\\.rds$", full.names = TRUE)
file_names <- gsub(pattern = "\\.rds$", replacement = "", x = basename(file_paths))

# Read each .rds file and store in a named list
datasets <- setNames(lapply(file_paths, readRDS), file_names)


#-------------------------------------------------------------------------------
# FUNCTION DEFINED 
#-------------------------------------------------------------------------------


# Function for main benchmarking ------
# Testing
input_data <- datasets$CC_R1_DL10_CT2_CN2000
project_name <- "CC_R1_DL10_CT2_CN2000"

input_data <- datasets$A549_processed
project_name <- "CC_R1_DL10_CT2_CN2000"

benchmark <- function(
    input_data, 
    project_name,        # string with dataset identifier 
    organism = "Hsap",   # Either Hsap or Mmus
    model_method = NULL, # Adjustable, alternative to specify the model type "linear", "spline", "polynomial", or "one_dimensional"
    output_path = "/Users/alicen/Projects/Damage_analsyis/test_strategies/strategy_output"   # Can be altered      
){
  
  message("\nBegin testing ", project_name, "...")
  
  # Tool 1: DamageDetective -------
  message("Tool 1: DamageDetective ...")
  DamageDetective_start <- Sys.time()
  
  limiric_output <- limiric(
    project_name = project_name,
    filter_rbc   = FALSE,
    seurat_input = input_data,
    filter_output = FALSE,
    output_path  = tempdir()
  )
  
  # Add to Seurat 
  input_data$DamageDetective <- limiric_output$limiric
  
  DamageDetective_end <- Sys.time()
  
  # Record time 
  DamageDetective_duration <- as.numeric(difftime(DamageDetective_end, DamageDetective_start, units = "secs"))
  
  # Tool 2: ddqc -------
  message("Tool 2: ddqc ...")
  ddqc_start <- Sys.time()
  
  # Requires seurat object with named mito and ribosomal percentage columns
  input_data$percent.mt <- input_data$mt.percent
  input_data$percent.rb <- input_data$rb.percent
  ddqc_results <- ddqc.metrics(input_data) 

  # Add to Seurat 
  input_data$ddqc <- ddqc_results$passed.qc[match(rownames(input_data@meta.data), rownames(ddqc_results))]
  input_data$ddqc <- ifelse(input_data$ddqc == "FALSE", "damaged", "cell")

  ddqc_end <- Sys.time()
  
  # Record time 
  ddqc_duration <- as.numeric(difftime(ddqc_end, ddqc_start, units = "secs"))
  
  # Tool 3: DropletQC ------
  message("Tool 3: DropletQC ...")
  
  DropletQC_start_1 <- Sys.time()

  if ("nf" %in% colnames(input_data@meta.data)){
    
    # Extract nf meta data & associated cell barcode from Seurat object
    edDf <- data.frame(nf = as.numeric(input_data$nf), umi = input_data$nCount_RNA)
    nf_col <- as.numeric(input_data$nf)
    
  } else {
    
    # Extract nf alternative nf_malat1 if nf absent & associated cell barcode from Seurat object
    edDf <- data.frame(nf = as.numeric(input_data$nf_malat1), umi = input_data$nCount_RNA)
    nf_col <- as.numeric(input_data$nf_malat1)
    
  }
  
  # Use droplet_qc function to identify empty droplets
  edDf$umi <- as.integer(edDf$umi)
  edresultsDf <- identify_empty_drops(edDf)
  DropletQC_stop_1 <- Sys.time()
  

  # Generate rough cluster labels 
  temp <- NormalizeData(input_data, verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    FindNeighbors(dims = 1:10, verbose = FALSE) %>%
    FindClusters(verbose = FALSE, resolution = 0.1) %>% 
    RunUMAP(dims = 1:10)
  
  # If clusters are too small, mark to combine with others
  celltype_table <- table(temp$seurat_clusters)
  combine <- as.integer(names(celltype_table)[celltype_table < 200])
  remaining_clusters <- setdiff(as.integer(names(celltype_table)), combine)
  combine_with <- remaining_clusters[which.min(celltype_table[as.character(remaining_clusters)])]
  temp$seurat_clusters <- ifelse(temp$seurat_clusters %in% combine, combine_with, temp$seurat_clusters)

  # Extract cluster assignments 
  DropletQC_start_2 <- Sys.time()
  cell_type <- temp$seurat_clusters 
  
  # Create input data frame 
  dcDf <- data.frame(
    nf = nf_col,
    umi = as.integer(input_data$nCount_RNA),
    cell_status = edresultsDf$cell_status,
    cell_type = cell_type
  )
  
  # Run damaged detection 
  dcresultsDf <- identify_damaged_cells(nf_umi_ed_ct = dcDf)
  input_data$DropletQC <- dcresultsDf$df$cell_status[match(rownames(dcresultsDf$df), rownames(input_data@meta.data))]
  
  # Two versions of DropletQC, one leaving empty droplets in & other not 
  input_data$DropletQC_empty <- input_data$DropletQC
  input_data$DropletQC <- ifelse(input_data$DropletQC == "damaged_cell", "damaged", "cell")
  DropletQC_stop_2 <- Sys.time()
  
  # Calculate timing
  DropletQC_duration_1 <- as.numeric(difftime(DropletQC_stop_1, DropletQC_start_1, units = "secs"))
  DropletQC_duration_2 <- as.numeric(difftime(DropletQC_stop_2, DropletQC_start_2, units = "secs"))
  DropletQC_duration <- DropletQC_duration_1 + DropletQC_duration_2
  
  # Tool 4: miQC -------
  message("Tool 4: miQC ...")
  # Start timing the miQC process
  miQC_start <- Sys.time()
  
  # Extract counts from the input data
  counts <- input_data@assays$RNA
  
  # Create a SingleCellExperiment object
  sce <- SingleCellExperiment(list(counts = counts))
  
  # Assign feature types and mitochondrial percentages
  rowData(sce)$featureType <- rownames(input_data@assays$RNA)
  colData(sce)$subsets_mito_percent <- input_data$mt.percent
  mainExpName(sce) <- 'gene'
  colData(sce)$detected <- input_data$nFeature_RNA
  
  # Function to run miQC and count the number of cells marked as 'damaged'
  evaluate_model <- function(sce, model_type) {
    tryCatch({
      # Fit the mixture model
      model <- mixtureModel(sce, model_type)
      
      # Filter cells based on the model
      sce_filtered <- filterCells(sce, model)
      
      # Determine the number of cells marked as 'damaged'
      damaged_cells <- setdiff(colnames(sce), colnames(sce_filtered))
      num_damaged <- length(damaged_cells)
      
      return(list(model = model, num_damaged = num_damaged))
    }, error = function(e) {
      message("Error in evaluate_model with model_type = ", model_type, ": ", e$message)
      stop("Evaluation failed, setting all cells to 'cell'")
    })
  }
  
  # Define the model types to evaluate
  model_types <- c("linear", "spline", "polynomial", "one_dimensional")
  
  # Evaluate each model type
  results <- list()
  for (model_type in model_types) {
    result <- tryCatch({
      evaluate_model(sce, model_type = model_type)
    }, error = function(e) {
      # If an error occurs, set miQC to "cell" and exit immediately
      input_data$miQC <- "cell"
      message("All cells set to 'cell' due to error in model evaluation.")
      return(NULL)
    })
    
    # If no error, store the result
    if (!is.null(result)) {
      results[[model_type]] <- result
    } else {
      break  # Stop processing further if an error occurred
    }
  }
  
  # Proceed only if results were obtained
  if (length(results) > 0) {
    # Extract the number of damaged cells for each model
    num_damaged_cells <- sapply(results, function(x) x$num_damaged)
    
    # Select the model that removes the fewest cells
    best_model_type <- names(which.min(num_damaged_cells))
    
    # If a specific model method is provided, override the selection
    if (!is.null(model_method)) {
      best_model_type <- model_method
    }
    
    # Apply the best model to filter cells
    message("\nTool 4: miQC ...    ", best_model_type)
    best_model <- results[[best_model_type]]$model
    sce_subset <- filterCells(sce, best_model)
    miQC_cells <- colnames(sce_subset)
    
    # Update the input data with miQC results
    input_data$miQC <- ifelse(rownames(input_data@meta.data) %in% miQC_cells, "cell", "damaged")
    
    # Optional: plot metrics if required
    if (view_plot) {
      print(plotModel(sce))
    }
    
  } else {
    # In case all models failed, set miQC to "cell"
    message("All miQC models failed.")
    input_data$miQC <- "cell"
  }
  
  
  
  miQC_end <- Sys.time()
  miQC_duration <- as.numeric(difftime(miQC_end, miQC_start, units = "secs"))
  
  
  # Tool 5 & 6: scater -------
  message("Tool 5: scater ...")

  # Like miQC, scater requires single cell experiment object (sce) input 
  scater_PCA_start <- Sys.time()
  sce <- as.SingleCellExperiment(input_data)
  sce <- addPerCellQCMetrics(sce)
  
  # Automated outlier labeling for low quality cells 
  scater <- runColDataPCA(sce,
                          ncomponents = 2, 
                          variables = c("nCount_RNA", "nFeature_RNA", "mt.percent", "rb.percent"), 
                          outliers = TRUE)
  
  input_data$scater_PCA <- scater$outlier
  input_data$scater_PCA <- ifelse(input_data$scater_PCA == "TRUE", "damaged", "cell")
  
  # Record time
  scater_PCA_end <- Sys.time()
  scater_PCA_duration <- as.numeric(difftime(scater_PCA_end, scater_PCA_start, units = "secs"))
    
  
  # Alternative scater function
  scater_Filter_start <- Sys.time()
  sce <- as.SingleCellExperiment(input_data)
  sce <- addPerCellQCMetrics(sce)
  perCellQCFilters <- perCellQCFilters(sce, sub.fields = c("nCount_RNA", "nFeature_RNA", "mt.percent", "rb.percent"))
  input_data$scater_Filter <- perCellQCFilters$discard
  input_data$scater_Filter <- ifelse(input_data$scater_Filter == "TRUE", "damaged", "cell")
  
  
  # Record time
  scater_Filter_end <- Sys.time()
  scater_Filter_duration <- as.numeric(difftime(scater_Filter_end, scater_Filter_start, units = "secs"))
  
  
  
  # Tool 7: valiDrops -------
  message("Tool 6: valiDrops ...")
  
  # Extract matrix from Seurat object
  valiDrops_start <- Sys.time()
  expression_matrix <- GetAssayData(input_data, layer = "counts")
  
  if (organism == "Hsap") {species <- "human"}
  if (organism == "Mmus") {species <- "mouse"}
  
  # Run valiDrops internal function to get list of quality control metrics
  expression_metrics <- quality_metrics(counts = expression_matrix,
                                        species = species)
  
  
  # Run valiDrops damaged detection (label_dead valiDrops function edited to allow for this input, not previous ValiDrop input)
  valiDrops <- label_dead(counts = expression_matrix, 
                          metrics = expression_metrics$metrics)
  
  # Add to benchmark seurat object
  input_data$valiDrops <- valiDrops$label[match(rownames(input_data@meta.data), valiDrops$metrics$barcode)]
  input_data$valiDrops <- ifelse(input_data$valiDrops == "dead", "damaged", "cell")
  
  # Record time
  valiDrops_end <- Sys.time()
  valiDrops_duration <- as.numeric(difftime(valiDrops_end, valiDrops_start, units = "secs"))
  
  # Manual methods ----
  
  # Manual 1: Mito fixed threshold ------
  input_data$manual_fixed_mito <- ifelse(input_data$mt.percent > 10, "damaged", "cell")
  
  # Preparation for outlier-based manual methods 
  # Extract df, add cell barcode identifies, & transform library size values for better Mahalanobis distance calculations 
  input_data_df <- input_data@meta.data
  input_data_df$barcodes <- rownames(input_data_df)
  input_data_df$nCount_RNA <- log10(input_data_df$nCount_RNA)
  input_data_df$nFeature_RNA <- log10(input_data_df$nFeature_RNA)
  
  # Transform malat1 counts for cases where they are close to zero 
  if (project_name %in% c("hLiver")){
    
    input_data_df$malat1.percent <- log10(input_data_df$malat1.percent)
    
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
  uni_mito_df  <- input_data_df[, c("mt.percent")]
  median_mito <- median(uni_mito_df , na.rm = TRUE)
  mad_mito <- mad(uni_mito_df , constant = 1, na.rm = TRUE)  # Use constant = 1 to follow the typical MAD definition
  # threshold_percentile <- quantile(input_data$mt.percent, probs = 0.90, na.rm = TRUE)
  
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
  # if (threshold_mad >= threshold_percentile){ 
  #   threshold = threshold_mad
  # } else {
  #   threshold = threshold_percentile
  # }
  
  # Use threshold to find outliers
  outliers_mito <- uni_mito_df  > threshold_mad
  uni_mito_df  <- as.data.frame(input_data_df[, c("mt.percent")])
  rownames(uni_mito_df ) <- rownames(input_data_df)
  uni_mito_results <- data.frame(barcodes = rownames(uni_mito_df), outliers = outliers_mito)
  uni_mito_results$outliers <- ifelse(uni_mito_results$outliers == "TRUE", "damaged", "cell")
  input_data$manual_adaptive_mito <- uni_mito_results$outliers
  
  # Manual 3: Mito/Ribo intersecting outliers (MAD): mitochondrial and ribosomal percentages -------
  
  # Adjust mito. to include ribo.
  
  # Use median and MAD to define threshold 
  uni_ribo_df  <- input_data_df[, c("rb.percent")]
  median_ribo <- median(uni_ribo_df , na.rm = TRUE)
  mad_ribo <- mad(uni_ribo_df , constant = 1, na.rm = TRUE)  # Use constant = 1 to follow the typical MAD definition
  # threshold_percentile <- quantile(input_data$rb.percent, probs = 0.10, na.rm = TRUE)
  
  # Adjust for liver case 
  if (project_name == "hLiver"){
    
    # Exception: No threshold could identify the minimum number of outliers.
    threshold_mad = 0
    
  } else {
    
    # Calculate the thresholds 
    dynamic_constant <- find_dynamic_mad_constant( uni_ribo_df, direction = "lower", min_outliers = 5)
    threshold_mad <- median_ribo -  dynamic_constant * mad_ribo
    
  }
  
  # Use threshold to find outliers
  outliers_ribo <- uni_ribo_df < threshold_mad
  uni_ribo_df  <- as.data.frame(input_data_df[, c("rb.percent")])
  rownames(uni_ribo_df ) <- rownames(input_data_df)
  uni_ribo_results <- data.frame(barcodes = rownames(uni_ribo_df), outliers = outliers_ribo)
  uni_ribo_results$outliers <- ifelse(uni_ribo_results$outliers == "TRUE", "damaged", "cell")
  
  # Add ribo's to adjust mito outlier 
  input_data$manual_mito_ribo <- uni_ribo_results$outliers
  input_data$manual_mito_ribo <- ifelse(input_data$manual_adaptive_mito == "damaged" & input_data$manual_mito_ribo == "damaged", "damaged", "cell")
  
  
  # Manual 4: Mito/Ribo/UMI/feature intersecting outliers (MAD) -------
  
  # Define broader percentile threshold for UMI 
  uni_count_df  <- input_data_df[, c("nCount_RNA")]
  median_count <- median(uni_count_df , na.rm = TRUE)
  mad_count <- mad(uni_count_df , constant = 1, na.rm = TRUE)  # Use constant = 1 to follow the typical MAD definition
  dynamic_constant <- find_dynamic_mad_constant(uni_count_df, direction = "lower", min_outliers = 5)
  threshold_MAD <- median_ribo - dynamic_constant * mad_ribo
  
  
  # Use threshold to find outliers 
  outliers_count <- uni_count_df  < threshold_MAD
  uni_count_df  <- as.data.frame(input_data_df[, c("nCount_RNA")])
  rownames(uni_count_df ) <- rownames(input_data_df)
  uni_count_results <- data.frame(barcodes = rownames(uni_count_df), outliers = outliers_count)
  uni_count_results$outliers <- ifelse(uni_count_results$outliers == "TRUE", "damaged", "cell")
  
  # Add ribo's to adjust mito outlier 
  input_data$manual_mito_ribo_umi <- uni_count_results$outliers
  input_data$manual_mito_ribo_umi <- ifelse(input_data$manual_mito_ribo == "damaged" & input_data$manual_mito_ribo_umi == "damaged", "damaged", "cell")
  
  
  # Repeat for features
  uni_feature_df  <- input_data_df[, c("nFeature_RNA")]
  median_feature <- median(uni_feature_df , na.rm = TRUE)
  mad_feature <- mad(uni_feature_df , constant = 1, na.rm = TRUE)  # Use constant = 1 to follow the typical MAD definition
  dynamic_constant <- find_dynamic_mad_constant(uni_feature_df, direction = "lower", min_outliers = 10)
  threshold <- median_feature - dynamic_constant * mad_feature
  
  # Use threshold to find outliers 
  outliers_feature <- uni_feature_df  < threshold
  uni_feature_df  <- as.data.frame(input_data_df[, c("nFeature_RNA")])
  rownames(uni_feature_df ) <- rownames(input_data_df)
  uni_feature_results <- data.frame(barcodes = rownames(uni_feature_df), outliers = outliers_feature)
  uni_feature_results$outliers <- ifelse(uni_feature_results$outliers == "TRUE", "damaged", "cell")
  
  # Add ribo's to adjust mito outlier 
  input_data$manual_mito_ribo_feature <- uni_feature_results$outliers
  input_data$manual_mito_ribo_library <- ifelse(input_data$manual_mito_ribo_umi == "damaged" & input_data$manual_mito_ribo_feature == "damaged", "damaged", "cell")
  
  
  # Manual 5: UMI/feature intersecting outliers (MAD) -------
  
  # Isolated for count and feature outliers 
  input_data$manual_umi <- uni_count_results$outliers
  input_data$manual_feature <- uni_feature_results$outliers
  input_data$manual_library <- ifelse(input_data$manual_umi == "damaged" | input_data$manual_feature == "damaged", "damaged", "cell")
  
  # Manual 6: MALAT1 univariate outliers (MAD) ------
  
  # Use median and MAD to define threshold 
  uni_malat1_df  <- input_data_df[, c("malat1")]
  median_malat1 <- median(uni_malat1_df , na.rm = TRUE)
  mad_malat1 <- mad(uni_malat1_df , constant = 1, na.rm = TRUE)  # Use constant = 1 to follow the typical MAD definition
  
  # Adjust for liver case 
  if (project_name == "hLiver"){
    
    # Exception: No threshold could identify the minimum number of outliers.
    threshold_mad = 0
    
  } else {
    
    # Calculate the thresholds 
    dynamic_constant <- find_dynamic_mad_constant( uni_malat1_df, direction = "upper", min_outliers = 5)
    threshold_mad <- median_malat1 +  dynamic_constant * mad_malat1
    
  }
  
  # Use threshold to find outliers
  outliers_malat1 <- uni_malat1_df  > threshold_mad
  uni_malat1_df  <- as.data.frame(input_data_df[, c("malat1")])
  rownames(uni_malat1_df) <- rownames(input_data_df)
  uni_malat1_results <- data.frame(barcodes = rownames(uni_malat1_df), outliers = outliers_malat1)
  uni_malat1_results$outliers <- ifelse(uni_malat1_results$outliers == "TRUE", "damaged", "cell")
  
  # Add outlier labels to cells in the input_data object 
  input_data$manual_malat1 <- uni_malat1_results$outliers
  
  # Manual 7: MALAT1/Mito/Ribo intersecting outliers (MAD) ----
  
  # Adjusting the MALAT1 predictions by other measures 
  input_data$manual_malat1_mito <- ifelse(input_data$manual_malat1 == "damaged" & input_data$manual_adaptive_mito == "damaged", "damaged", "cell")
  
  # Return output -------
  
  # Call helper function 
  final_df <- summarise_results(seurat = input_data)
  
  write.csv(final_df, 
            paste0(output_path, "/", project_name, "_summary.csv"), 
            quote = FALSE, 
            row.names = FALSE) 
  
  write.csv(input_data@meta.data, 
            paste0(output_path, "/", project_name, ".csv"), 
            quote = FALSE, 
            row.names = TRUE) 
  
  
  message("Complete!", "\n")
  
  # Global environment output 
  return(input_data)
  
}


#-------------------------------------------------------------------------------
# FUNCTION RUNNING
#-------------------------------------------------------------------------------

# Define mouse samples
mouse_samples <- c("mLiver_processed", "mLung_processed", "mPBMC_processed")

# Apply the benchmark function to each dataset in the results list
benchmark_results <- lapply(names(datasets), function(sample_name) {
  
  # Retrieve the dataset
  input_data <- datasets[[sample_name]]
  
  # Determine the organism based on the sample name
  organism <- if (sample_name %in% mouse_samples) "Mmus" else "Hsap"
  
  # Run the benchmark function
  benchmark(
    input_data = input_data,
    project_name = sample_name,
    organism = organism
  )
  
})



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
