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
# (Isolated to allow for adjustments if necessary)

# Parent directories 
ddqc_path <- "/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/python_tool_output/ddqc_output/"
ensemble_path <- "/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/python_tool_output/EnsembleKQC_simulated_output/"
output_path <- "/home/alicen/Projects/ReviewArticle/damage_perturbation/benchmark_results/"


# 2.5 % damaged 
benchmark(seurat = data_list$control_sim_1_2.5, 
          project_name = "control_sim_1_2.5",
          ddqc_path = paste0(ddqc_path, "control_sim_1_2.5.csv"), 
          ensembleKQC_path = paste0(ensemble_path, "control_sim_1_2.5_input.csv"),
          output_path = output_path)

benchmark(seurat = data_list$control_sim_2_2.5, 
          project_name = "control_sim_2_2.5",
          ddqc_path = paste0(ddqc_path, "control_sim_2_2.5.csv"), 
          ensembleKQC_path = paste0(ensemble_path, "control_sim_2_2.5_input.csv"),
          output_path = output_path)

benchmark(seurat = data_list$control_sim_3_2.5, 
          project_name = "control_sim_3_2.5",
          ddqc_path = paste0(ddqc_path, "control_sim_3_2.5.csv"), 
          ensembleKQC_path = paste0(ensemble_path, "control_sim_3_2.5_input.csv"),
          output_path = output_path)

benchmark(seurat = data_list$stimulated_sim_1_2.5, 
          project_name = "stimulated_sim_1_2.5",
          ddqc_path = paste0(ddqc_path, "stimulated_sim_1_2.5.csv"), 
          ensembleKQC_path = paste0(ensemble_path, "stimulated_sim_1_2.5_input.csv"),
          output_path = output_path)

benchmark(seurat = data_list$stimulated_sim_2_2.5, 
          project_name = "stimulated_sim_2_2.5",
          ddqc_path = paste0(ddqc_path, "stimulated_sim_2_2.5.csv"), 
          ensembleKQC_path = paste0(ensemble_path, "stimulated_sim_2_2.5_input.csv"),
          output_path = output_path)

benchmark(seurat = data_list$stimulated_sim_3_2.5, 
          project_name = "stimulated_sim_3_2.5",
          ddqc_path = paste0(ddqc_path, "stimulated_sim_3_2.5.csv"), 
          ensembleKQC_path = paste0(ensemble_path, "stimulated_sim_3_2.5_input.csv"),
          output_path = output_path)



# 5 % damaged 
benchmark(seurat = data_list$control_sim_1_5, 
          project_name = "control_sim_1_5",
          ddqc_path = paste0(ddqc_path, "control_sim_1_5.csv"), 
          ensembleKQC_path = paste0(ensemble_path, "control_sim_1_5_input.csv"),
          output_path = output_path)

benchmark(seurat = data_list$control_sim_2_5, 
          project_name = "control_sim_2_5",
          ddqc_path = paste0(ddqc_path, "control_sim_2_5.csv"), 
          ensembleKQC_path = paste0(ensemble_path, "control_sim_2_5_input.csv"),
          output_path = output_path)

benchmark(seurat = data_list$control_sim_3_5, 
          project_name = "control_sim_3_5",
          ddqc_path = paste0(ddqc_path, "control_sim_3_5.csv"), 
          ensembleKQC_path = paste0(ensemble_path, "control_sim_3_5_input.csv"),
          output_path = output_path)

benchmark(seurat = data_list$stimulated_sim_1_5, 
          project_name = "stimulated_sim_1_5",
          ddqc_path = paste0(ddqc_path, "stimulated_sim_1_5.csv"), 
          ensembleKQC_path = paste0(ensemble_path, "stimulated_sim_1_5_input.csv"),
          output_path = output_path)

benchmark(seurat = data_list$stimulated_sim_2_5, 
          project_name = "stimulated_sim_2_5",
          ddqc_path = paste0(ddqc_path, "stimulated_sim_2_5.csv"), 
          ensembleKQC_path = paste0(ensemble_path, "stimulated_sim_2_5_input.csv"),
          output_path = output_path)

benchmark(seurat = data_list$stimulated_sim_3_5, 
          project_name = "stimulated_sim_3_5",
          ddqc_path = paste0(ddqc_path, "stimulated_sim_3_5.csv"), 
          ensembleKQC_path = paste0(ensemble_path, "stimulated_sim_3_5_input.csv"),
          output_path = output_path)


# 10 % damaged 
benchmark(seurat = data_list$control_sim_1_10, 
          project_name = "control_sim_1_10",
          ddqc_path = paste0(ddqc_path, "control_sim_1_10.csv"), 
          ensembleKQC_path = paste0(ensemble_path, "control_sim_1_10_input.csv"),
          output_path = output_path)

benchmark(seurat = data_list$control_sim_2_10, 
          project_name = "control_sim_2_10",
          ddqc_path = paste0(ddqc_path, "control_sim_2_10.csv"), 
          ensembleKQC_path = paste0(ensemble_path, "control_sim_2_10_input.csv"),
          output_path = output_path)

benchmark(seurat = data_list$control_sim_3_10, 
          project_name = "control_sim_3_10",
          ddqc_path = paste0(ddqc_path, "control_sim_3_10.csv"), 
          ensembleKQC_path = paste0(ensemble_path, "control_sim_3_10_input.csv"),
          output_path = output_path)

benchmark(seurat = data_list$stimulated_sim_1_10, 
          project_name = "stimulated_sim_1_10",
          ddqc_path = paste0(ddqc_path, "stimulated_sim_1_10.csv"), 
          ensembleKQC_path = paste0(ensemble_path, "stimulated_sim_1_10_input.csv"),
          output_path = output_path)

benchmark(seurat = data_list$stimulated_sim_2_10, 
          project_name = "stimulated_sim_2_10",
          ddqc_path = paste0(ddqc_path, "stimulated_sim_2_10.csv"), 
          ensembleKQC_path = paste0(ensemble_path, "stimulated_sim_2_10_input.csv"),
          output_path = output_path)

benchmark(seurat = data_list$stimulated_sim_3_10, 
          project_name = "stimulated_sim_3_10",
          ddqc_path = paste0(ddqc_path, "stimulated_sim_3_10.csv"), 
          ensembleKQC_path = paste0(ensemble_path, "stimulated_sim_3_10_input.csv"),
          output_path = output_path)

# 15 % damaged 
benchmark(seurat = data_list$control_sim_1_15, 
          project_name = "control_sim_1_15",
          ddqc_path = paste0(ddqc_path, "control_sim_1_15.csv"), 
          ensembleKQC_path = paste0(ensemble_path, "control_sim_1_15_input.csv"),
          output_path = output_path)

benchmark(seurat = data_list$control_sim_2_15, 
          project_name = "control_sim_2_15",
          ddqc_path = paste0(ddqc_path, "control_sim_2_15.csv"), 
          ensembleKQC_path = paste0(ensemble_path, "control_sim_2_15_input.csv"),
          output_path = output_path)

benchmark(seurat = data_list$control_sim_3_15, 
          project_name = "control_sim_3_15",
          ddqc_path = paste0(ddqc_path, "control_sim_3_15.csv"), 
          ensembleKQC_path = paste0(ensemble_path, "control_sim_3_15_input.csv"),
          output_path = output_path)

benchmark(seurat = data_list$stimulated_sim_1_15, 
          project_name = "stimulated_sim_1_15",
          ddqc_path = paste0(ddqc_path, "stimulated_sim_1_15.csv"), 
          ensembleKQC_path = paste0(ensemble_path, "stimulated_sim_1_15_input.csv"),
          output_path = output_path)

benchmark(seurat = data_list$stimulated_sim_2_15, 
          project_name = "stimulated_sim_2_15",
          ddqc_path = paste0(ddqc_path, "stimulated_sim_2_15.csv"), 
          ensembleKQC_path = paste0(ensemble_path, "stimulated_sim_2_15_input.csv"),
          output_path = output_path)

benchmark(seurat = data_list$stimulated_sim_3_15, 
          project_name = "stimulated_sim_3_15",
          ddqc_path = paste0(ddqc_path, "stimulated_sim_3_15.csv"), 
          ensembleKQC_path = paste0(ensemble_path, "stimulated_sim_3_15_input.csv"),
          output_path = output_path)


# 20 % damaged 
benchmark(seurat = data_list$control_sim_1_20, 
          project_name = "control_sim_1_20",
          ddqc_path = paste0(ddqc_path, "control_sim_1_20.csv"), 
          ensembleKQC_path = paste0(ensemble_path, "control_sim_1_20_input.csv"),
          output_path = output_path)

benchmark(seurat = data_list$control_sim_2_20, 
          project_name = "control_sim_2_20",
          ddqc_path = paste0(ddqc_path, "control_sim_2_20.csv"), 
          ensembleKQC_path = paste0(ensemble_path, "control_sim_2_20_input.csv"),
          output_path = output_path)

benchmark(seurat = data_list$control_sim_3_20, 
          project_name = "control_sim_3_20",
          ddqc_path = paste0(ddqc_path, "control_sim_3_20.csv"), 
          ensembleKQC_path = paste0(ensemble_path, "control_sim_3_20_input.csv"),
          output_path = output_path)

benchmark(seurat = data_list$stimulated_sim_1_20, 
          project_name = "stimulated_sim_1_20",
          ddqc_path = paste0(ddqc_path, "stimulated_sim_1_20.csv"), 
          ensembleKQC_path = paste0(ensemble_path, "stimulated_sim_1_20_input.csv"),
          output_path = output_path)

benchmark(seurat = data_list$stimulated_sim_2_20, 
          project_name = "stimulated_sim_2_20",
          ddqc_path = paste0(ddqc_path, "stimulated_sim_2_20.csv"), 
          ensembleKQC_path = paste0(ensemble_path, "stimulated_sim_2_20_input.csv"),
          output_path = output_path)

benchmark(seurat = data_list$stimulated_sim_3_20, 
          project_name = "stimulated_sim_3_20",
          ddqc_path = paste0(ddqc_path, "stimulated_sim_3_20.csv"), 
          ensembleKQC_path = paste0(ensemble_path, "stimulated_sim_3_20_input.csv"),
          output_path = output_path)


### End 