# SCRIPT CONTEXT 
#
# Running damaged cell detection on simulated datasets. This process is identical
# to that in part B of the analysis, the only difference being the simulated nature 
# of the input data. 
#
#
# Tools tested : 
# 1. ddqc (loading results)
# 2. ensembleKQC (loading results)
# 3. DropletQC : package functions 
# 4. miQC : select best model
# 5. scater PCA 
# 6. valiDrops : modified function 
# 7. SampleQC (run collectively on all)
#
#
# NB: 
# - Output from ddqc and ensembleKQC must be done
#     - 2.1 - Run ddqc.ipynb
#     - 2.2 - Run ensembleKQC.sh
# - All previous functions must be in global environment from previous scripts, they are not defined here.
#     - 3. Run all damaged detection strategies on all data
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

# Read in 30 simulated datasets using loop ----

# Parent directories (save output in list)
parent_directory <- "/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/R_tool_input/"
conditions <- c("control_sim", "stimulated_sim")
percentages <- c("2.5", "5", "10", "15", "20")
reps <- 1:3

# Store the Seurat objects in a list 
data_list <- list()

# Loop through conditions, percentages, and replicates
for (condition in conditions) {
  for (percentage in percentages) {
    for (rep in reps) {
      # Construct the file name
      file_name <- paste0(condition, "_", rep, "_", percentage, ".rds")
      file_path <- file.path(parent_directory, file_name)
      
      # Check if the file exists before reading
      if (file.exists(file_path)) {
        
        # Read the seurat object 
        data <- readRDS(file_path)
        
        # Store the data in the list with a meaningful name
        data_list[[paste0(condition, "_", rep, "_", percentage)]] <- data
      }
    }
  }
}


#-------------------------------------------------------------------------------
# RUN FUNCTION
#-------------------------------------------------------------------------------

# Tool testing on each dataset ----

# Parent directories 
ddqc_parent_path <- "/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/python_tool_output/ddqc_output/"
ensemble_parent_path <- "/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/python_tool_output/EnsembleKQC_simulated_output/"
output_parent_path <- "/home/alicen/Projects/ReviewArticle/damage_perturbation/benchmark_results/"

# Define sample types and damage percentages
sample_types <- c("control_sim", "stimulated_sim")
replicates <- 1:3
damage_percentages <- c(2.5, 5, 10, 15, 20)

# Initialize an empty list to store the results
simulated <- list()

# Loop over each sample type, replicate, and damage percentage
for (sample_type in sample_types) {
  for (replicate in replicates) {
    for (damage in damage_percentages) {
      # Construct the project name
      project_name <- paste0(sample_type, "_", replicate, "_", damage)
      
      # Construct the file paths
      ddqc_path <- paste0(ddqc_parent_path, project_name, ".csv")
      ensembleKQC_path <- paste0(ensemble_parent_path, project_name, "_input.csv")
      
      # Call the benchmark function
      result <- benchmark(
        seurat = data_list[[project_name]], 
        project_name = project_name,
        ddqc_path = ddqc_path, 
        ensembleKQC_path = ensembleKQC_path,
        output_path = output_parent_path
      )
      
      # Save the result to the list
      simulated[[project_name]] <- result
    }
  }
}



# SampleQC for each collection of samples -----

# Create a single object for all samples 
simulated_merged <- merge(x = simulated$control_sim_1_2.5,
                          y = c(simulated$control_sim_1_5, simulated$control_sim_1_10, simulated$control_sim_1_15, simulated$control_sim_1_20,
                                simulated$control_sim_2_2.5, simulated$control_sim_2_5, simulated$control_sim_2_10, simulated$control_sim_2_15, simulated$control_sim_2_20,
                                simulated$control_sim_3_2.5, simulated$control_sim_3_5, simulated$control_sim_3_10, simulated$control_sim_3_15, simulated$control_sim_3_20,
                                simulated$stimulated_sim_1_2.5, simulated$stimulated_sim_1_5, simulated$stimulated_sim_1_10, simulated$stimulated_sim_1_15, simulated$stimulated_sim_1_20, 
                                simulated$stimulated_sim_2_2.5, simulated$stimulated_sim_2_5, simulated$stimulated_sim_2_10, simulated$stimulated_sim_2_15, simulated$stimulated_sim_2_20, 
                                simulated$stimulated_sim_3_2.5, simulated$stimulated_sim_3_5, simulated$stimulated_sim_3_10, simulated$stimulated_sim_3_15, simulated$stimulated_sim_3_20),
                          add.cell.ids = names(simulated),
                     project = "Simulated"
)

# Convert to correct column name for sampleQC to recognise
simulated_merged$percent.mt <- simulated_merged$mt.percent


# Create SampleQC dataframe using default function 
qc_dt = make_qc_dt(simulated_merged@meta.data, 
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
simulated_merged$SampleQC <- outliers_dt$outlier
simulated_merged$SampleQC <- ifelse(simulated_merged$SampleQC == "TRUE", "damaged", "cell")


# Saving ----

for (dataset in names(simulated_merged)){
    
    # Extract and save human samples 
    seurat <- subset(non_groundtruth, orig.ident == dataset)
    save_sampleQC(seurat, as.character(dataset), organism = "Hsap")

}


### End 