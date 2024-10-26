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

# Load processed datasets ----

# 2.5 % damaged 
control_sim_1_2.5 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/R_tool_input/control_sim_1_2.5.rds")
control_sim_2_2.5 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/R_tool_input/control_sim_2_2.5.rds")
control_sim_2_2.5 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/R_tool_input/control_sim_2_2.5.rds")
stimulated_sim_1_2.5 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/R_tool_input/stimulated_sim_1_2.5.rds")
stimulated_sim_2_2.5 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/R_tool_input/stimulated_sim_2_2.5.rds")
stimulated_sim_3_2.5 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/R_tool_input/stimulated_sim_3_2.5.rds")

# 5 % damaged 
control_sim_1_5 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/R_tool_input/control_sim_1_5.rds")
control_sim_2_5 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/R_tool_input/control_sim_2_5.rds")
control_sim_2_5 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/R_tool_input/control_sim_2_5.rds")
stimulated_sim_1_5 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/R_tool_input/stimulated_sim_1_5.rds")
stimulated_sim_2_5 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/R_tool_input/stimulated_sim_2_5.rds")
stimulated_sim_3_5 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/R_tool_input/stimulated_sim_3_5.rds")

# 10 % damaged 
control_sim_1_10 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/R_tool_input/control_sim_1_10.rds")
control_sim_2_10 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/R_tool_input/control_sim_2_10.rds")
control_sim_2_10 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/R_tool_input/control_sim_2_10.rds")
stimulated_sim_1_10 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/R_tool_input/stimulated_sim_1_10.rds")
stimulated_sim_2_10 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/R_tool_input/stimulated_sim_2_10.rds")
stimulated_sim_3_10 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/R_tool_input/stimulated_sim_3_10.rds")

# 15 % damaged 
control_sim_1_15 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/R_tool_input/control_sim_1_15.rds")
control_sim_2_15 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/R_tool_input/control_sim_2_15.rds")
control_sim_2_15 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/R_tool_input/control_sim_2_15.rds")
stimulated_sim_1_15 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/R_tool_input/stimulated_sim_1_15.rds")
stimulated_sim_2_15 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/R_tool_input/stimulated_sim_2_15.rds")
stimulated_sim_3_15 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/R_tool_input/stimulated_sim_3_15.rds")

# 20 % damaged 
control_sim_1_20 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/R_tool_input/control_sim_1_20.rds")
control_sim_2_20 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/R_tool_input/control_sim_2_20.rds")
control_sim_2_20 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/R_tool_input/control_sim_2_20.rds")
stimulated_sim_1_20 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/R_tool_input/stimulated_sim_1_20.rds")
stimulated_sim_2_20 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/R_tool_input/stimulated_sim_2_20.rds")
stimulated_sim_3_20 <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2/R_tool_input/stimulated_sim_3_20.rds")



#-------------------------------------------------------------------------------
# TEST THE TOOLS 
#-------------------------------------------------------------------------------


# Use predefined functions to run the testing ----





