# SCRIPT CONTEXT 
#
# Benchmarking damaged cell detection strategies on ground truth dataset. This script begins with the strategies
# already run on the ground truth cases. The labelled output objects (3 - Run remaining detection strategies.R)
# are loaded at the start of the script. From here, the script uses the labels stored in the objects to 
# calculate confusion metrics (TP, TN, FP, FN) for each strategy, followed by sensitivity and F1 scores.
# The script then plots the sensitivity and F1 scores for each strategy
# 
# Still do : 
# - automated rather than manual confusion metric caculations 
# - automated precision, recall, PR-AUC calculations 
# - finalised plots (new aesthetics, labels, etc)


#-------------------------------------------------------------------------------
# PREPARATIONS
#-------------------------------------------------------------------------------

# Load necessary libraries
library(ggplot2)
library(dplyr)


# Load datasets (5 cases) -----

# HEK293 data
apoptotic <- readRDS("/home/alicen/Projects/limiric/damage_left_behind_analysis/groundtruth/apoptotic.rds")
pro_apoptotic <- readRDS("/home/alicen/Projects/limiric/damage_left_behind_analysis/groundtruth/pro_apoptotic.rds")

# SA928 
SA928_dead <- readRDS("/home/alicen/Projects/limiric/damage_left_behind_analysis/groundtruth/SA928_dead_live.rds")
SA928_dying <- readRDS("/home/alicen/Projects/limiric/damage_left_behind_analysis/groundtruth/SA928_dying_live.rds")

# SA604 
SA604_dead <- readRDS("/home/alicen/Projects/limiric/damage_left_behind_analysis/groundtruth/SA604_dead_live.rds")


#-------------------------------------------------------------------------------
# PERFORMANCE METRICS
#-------------------------------------------------------------------------------

# Confusion metrics  ----


# Adding measures of TP, TN, FP, FN
seurat <- SA604_dead
seurat$DropletQC <- ifelse(seurat$DropletQC == "empty_droplet", "cell", seurat$DropletQC)

# for SA928 
seurat$orig.ident <- rownames(seurat@meta.data)
seurat$orig.ident <- sub("__.*", "", seurat$orig.ident)


# String inputs 
TP <- c("dead")
TN <- c("live")


# List of methods
method <- c("manual_malat1")
# "ddqc", "DropletQC", "limiric", "miQC", "valiDrops", 
# "manual_all", "manual_mito_ribo", "manual_mito", "manual_malat1", "manual_mito_isolated")


seurat$method_outcome <- "na"
seurat$method_outcome <- ifelse(seurat$orig.ident == TP & seurat[[method]] == "damaged", "true_positive", seurat$method_outcome)
seurat$method_outcome <- ifelse(seurat$orig.ident == TP & seurat[[method]] == "cell", "false_negative", seurat$method_outcome)
seurat$method_outcome <- ifelse(seurat$orig.ident == TN & seurat[[method]] == "cell", "true_negative", seurat$method_outcome)
seurat$method_outcome <- ifelse(seurat$orig.ident == TN & seurat[[method]] == "damaged", "false_positive", seurat$method_outcome)

# Summarise for statistics 
outcome <- seurat@meta.data %>% 
  group_by(method_outcome) %>%
  summarise(Count = n())

outcome 

# Plot sensitivity & F1 --------


# Create the data frame
data <- data.frame(
  Dataset = c("HEK293_apoptotic", "HEK293_apoptotic", "HEK293_apoptotic", "HEK293_apoptotic", "HEK293_apoptotic", 
              "HEK293_apoptotic", "HEK293_apoptotic", "HEK293_apoptotic", "HEK293_apoptotic", "HEK293_apoptotic", 
              "HEK293_pro_apoptotic", "HEK293_pro_apoptotic", "HEK293_pro_apoptotic", "HEK293_pro_apoptotic", 
              "HEK293_pro_apoptotic", "HEK293_pro_apoptotic", "HEK293_pro_apoptotic", "HEK293_pro_apoptotic", 
              "HEK293_pro_apoptotic", "HEK293_pro_apoptotic", "SA928_dead", "SA928_dead", "SA928_dead", "SA928_dead", 
              "SA928_dead", "SA928_dead", "SA928_dead", "SA928_dead", "SA928_dead", "SA928_dead", "SA928_dying", 
              "SA928_dying", "SA928_dying", "SA928_dying", "SA928_dying", "SA928_dying", "SA928_dying", "SA928_dying", 
              "SA928_dying", "SA928_dying", "SA604_dead", "SA604_dead", "SA604_dead", "SA604_dead", "SA604_dead", 
              "SA604_dead", "SA604_dead", "SA604_dead", "SA604_dead", "SA604_dead"),
  Tool = c("ddqc", "DropletQC", "limirc", "miQC", "valiDrops", "all", "mito_ribo", "mito", 
           "mito_isolated", "malat1", "ddqc", "DropletQC", "limirc", "miQC", "valiDrops", "all", 
           "mito_ribo", "mito", "mito_isolated", "malat1", "ddqc", "DropletQC", "limirc", 
           "miQC", "valiDrops", "all", "mito_ribo", "mito", "mito_isolated", "malat1", 
           "ddqc", "DropletQC", "limirc", "miQC", "valiDrops", "all", "mito_ribo", "mito", 
           "mito_isolated", "malat1", "ddqc", "DropletQC", "limirc", "miQC", "valiDrops", "all", 
           "mito_ribo", "mito", "mito_isolated", "malat1"),
  Sensitivity = c(0.067956583, 0, 0.50023596, 0.5261916, 0.084945729, 0.555450684, 0.55261916, 0.489853705, 0.299669655, 
                  0.290231241, 0.084509353, 0.404988513, 0.24762061, 0.64292747, 0.001805054, 0.109123728, 0.270922219, 
                  0.29881851, 0.018214637, 0.171972432, 0.070691025, 0, 0.525814138, 0.673550437, 0.218427323, 0.5, 
                  0.889594917, 0.883240667, 0.659253376, 0.853852264, 0.049219688, 0, 0.204081633, 0.241296519, 
                  0.022809124, 0.356542617, 0.337334934, 0.326530612, 0.235294118, 0.2899106, 0.283333333, 0, 0.55, 
                  0.561111111, 0.561111111, 0.786111111, 0.775, 0.730555556, 0.575, 0.694444444),
  F1 = c(0.101802757, 0, 0.492679526, 0.375484088, 0.155709343, 0.422848931, 0.422896352, 0.38430211, 0.317897372, 
         0.302583026, 0.147670251, 0.457884972, 0.314801293, 0.625478927, 0.003503743, 0.172794595, 0.34203439, 
         0.3679903, 0.033474065, 0.252165544, 0.088778055, 0, 0.474551971, 0.427850656, 0.357607282, 0.374313758, 
         0.457236171, 0.472688629, 0.438573316, 0.466478629, 0.061148397, 0, 0.189732143, 0.137624101, 0.042316258, 
         0.151918159, 0.149547632, 0.155074116, 0.137737175, 0.134399053, 0.131274131, 0, 0.709677419, 0.58045977, 
         0.690598291, 0.33993994, 0.366863905, 0.352310784, 0.419452888, 0.441306267)
)

# Create box plots
sensitivity_plot <- ggplot(data, aes(x = Tool, y = Sensitivity)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, aes(color = Dataset)) +
  labs(title = "Sensitivity Scores", y = "Sensitivity") + 
  theme_classic() +  
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  )


f1_plot <- ggplot(data, aes(x = Tool, y = F1)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, aes(color = Dataset)) +
  labs(title = "F1 Scores", y = "F1 Score")  +
  theme_classic() +  
  theme(
    #legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  )

f1_plot + sensitivity_plot

# bar plots for median 

# Calculate median scores for each tool
median_scores <- data %>%
  group_by(Tool) %>%
  summarise(Median_Sensitivity = median(Sensitivity), Median_F1 = median(F1))

# Create bar plots for median scores
sensitivity_bar_plot <- ggplot(median_scores, aes(x = Tool, y = Median_Sensitivity)) +
  geom_bar(stat = "identity") +
  labs(title = "Median Sensitivity Scores", y = "Median Sensitivity") + 
  theme_classic() +  
  theme(
    #legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    axis.text = element_text(size = 16),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  )

f1_bar_plot <- ggplot(median_scores, aes(x = Tool, y = Median_F1)) +
  geom_bar(stat = "identity") +
  labs(title = "Median F1 Scores", y = "Median F1 Score")+ 
  theme_classic() +  
  theme(
    #legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    axis.text = element_text(size = 16),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  )

# Print the plots
print(sensitivity_bar_plot)
print(f1_bar_plot)


# False negative rate 

# Create the data frame
data <- data.frame(
  Dataset = c("HEK293_apoptotic", "HEK293_apoptotic", "HEK293_apoptotic", "HEK293_apoptotic", "HEK293_apoptotic", 
              "HEK293_apoptotic", "HEK293_apoptotic", "HEK293_apoptotic", "HEK293_apoptotic", "HEK293_apoptotic", 
              "HEK293_pro_apoptotic", "HEK293_pro_apoptotic", "HEK293_pro_apoptotic", "HEK293_pro_apoptotic", 
              "HEK293_pro_apoptotic", "HEK293_pro_apoptotic", "HEK293_pro_apoptotic", "HEK293_pro_apoptotic", 
              "HEK293_pro_apoptotic", "HEK293_pro_apoptotic", "SA928_dead", "SA928_dead", "SA928_dead", "SA928_dead", 
              "SA928_dead", "SA928_dead", "SA928_dead", "SA928_dead", "SA928_dead", "SA928_dead", "SA928_dying", 
              "SA928_dying", "SA928_dying", "SA928_dying", "SA928_dying", "SA928_dying", "SA928_dying", "SA928_dying", 
              "SA928_dying", "SA928_dying", "SA604_dead", "SA604_dead", "SA604_dead", "SA604_dead", "SA604_dead", 
              "SA604_dead", "SA604_dead", "SA604_dead", "SA604_dead", "SA604_dead"),
  Tool = c("ddqc", "DropletQC", "limirc", "miQC", "valiDrops", "all", "mito_ribo", "mito", 
           "mito_isolated", "malat1", "ddqc", "DropletQC", "limirc", "miQC", "valiDrops", "all", 
           "mito_ribo", "mito", "mito_isolated", "malat1", "ddqc", "DropletQC", "limirc", 
           "miQC", "valiDrops", "all", "mito_ribo", "mito", "mito_isolated", "malat1", 
           "ddqc", "DropletQC", "limirc", "miQC", "valiDrops", "all", "mito_ribo", "mito", 
           "mito_isolated", "malat1", "ddqc", "DropletQC", "limirc", "miQC", "valiDrops", "all", 
           "mito_ribo", "mito", "mito_isolated", "malat1"),
  FNR = c(0.932043417, 1, 0.49976404, 0.4738084, 0.915054271, 0.444549316, 0.44738084, 0.510146295, 0.700330345, 
          0.709768759, 0.915490647, 0.595011487, 0.75237939, 0.35707253, 0.998194946, 0.890876272, 0.729077781, 
          0.70118149, 0.981785363, 0.828027568, 0.929308975, 1, 0.474185862, 0.326449563, 0.781572677, 0.5, 
          0.110405083, 0.116759333, 0.340746624, 0.146147736, 0.950780312, 1, 0.795918367, 0.758703481, 0.977190876, 
          0.643457383, 0.662665066, 0.673469388, 0.764705882, 0.7100894, 0.716666667, 1, 0.45, 0.438888889, 
          0.438888889, 0.213888889, 0.225, 0.269444444, 0.425, 0.305555556)
)


# Calculate median FNR scores for each tool
median_scores <- data %>%
  group_by(Tool) %>%
  summarise(Median_FNR = median(FNR))

# Create bar plot for median FNR scores
fnr_bar_plot <- ggplot(median_scores, aes(x = Tool, y = Median_FNR)) +
  geom_bar(stat = "identity") +
  labs(title = "Median FNR", y = "Median FNR")+ 
  theme_classic() +  
  theme(
    #legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    axis.text = element_text(size = 16),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  )

fnr_bar_plot
