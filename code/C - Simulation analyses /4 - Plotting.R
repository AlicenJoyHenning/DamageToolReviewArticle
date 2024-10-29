# SCRIPT CONTEXT 
# 
# Scripts to generate plots from the simulated data.
#
# 1. Baseline HVG analysis 
# Control (undamaged) vs test (damaged added but unfiltered) cases: Overlap analysis to find mean intersection and isolated cases from 
# which the GO categories related to the isolated cases can be found.
#
# 2. HVG similarity for damaged cell strategies 
# Bar plots showing the similarity (jaccard index) of the HVGs 
# of the filtered cases compared to those of the control (undamaged) cases. 
#  
# 3. Baseline DEG correctness  
# Confusion matrices for control (undamaged) vs test (damaged added but unfiltered)cases 
#
# 4. DEG correctness for damaged cell strategies 
# Bar plots showing the correctness (F1 score) of celltype-specific DEGs 
# of the filtered cases compared to those of the control (undamaged) cases


#-------------------------------------------------------------------------------
# PREPARATIONS
#-------------------------------------------------------------------------------

# Load libraries -------

packages <- c("biomaRt", "clusterProfiler", "dplyr", "ggplot2", "glmGamPoi", "Matrix", "robustbase",
              "png", "Seurat", "tidyr", "SingleCellExperiment", "scuttle", "presto")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    library(pkg)
  }
}


#-------------------------------------------------------------------------------
# BASELINE HVG
#-------------------------------------------------------------------------------

# 1. Baseline HVG ----

HVGs <- readRDS("/home/alicen/Projects/ReviewArticle/damage_perturbation/analysis_results/HVGs_unique.rds")
HVG_overlap <- read.csv("/home/alicen/Projects/ReviewArticle/damage_perturbation/analysis_results/parent_v_unfiltered.csv")

# Plot correlation between damage percentage and jaccard scores ----
correlation <- cor(HVG_overlap$damage, HVG_overlap$jaccard_index)

# Create the plot
correlation_plot <- ggplot(HVG_overlap, aes(x = damage, y = jaccard_index)) +
  geom_point() +  
  geom_smooth(method = "lm", se = FALSE, color = "#99CDF6") +  
  labs(title = "",
       x = "Damage added (%)",
       y = "Jaccard Index") +
  annotate("text", x = 10, y = 0.95, 
           label = paste("R =", round(correlation, 2)), 
           size = 5, color = "black") +  
  theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill = NA))

ggsave("/home/alicen/Projects/ReviewArticle/damage_perturbation/analysis_results/correlation.png",
       plot = correlation_plot,
       width = 3, height = 3, units = "in")


# Venn ----
# Find median values for each extent of damage
hvg_parent_v_unfiltered$prefix <- sub("^(.*?_.*?_.*?)_.*$", "\\1", hvg_parent_v_unfiltered$case)

hvg_median_values <- hvg_parent_v_unfiltered %>%
  group_by(damage) %>%
  summarize(
    median_intersection = median(intersection, na.rm = TRUE),
    median_control_unique = median(control_unique, na.rm = TRUE),
    median_damage_unique = median(damage_unique, na.rm = TRUE),
    median_jaccard_index = median(jaccard_index, na.rm = TRUE)
  )

hvg_median_values$damage <- as.numeric(hvg_median_values$damage)
write.csv(hvg_median_values, 
          "/home/alicen/Projects/ReviewArticle/damage_perturbation/analysis_results/median_parent_v_unfiltered.csv",
          quote = FALSE, row.names = FALSE)

# Enriched GO terms in HVG sets ----

# Get GO terms for overlapping genes
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# damage
gene_lists <- list(
  HVGs$control_sim_3_20$damage_HVGs, 
  HVGs$control_sim_2_20$damage_HVGs, 
  HVGs$control_sim_1_20$damage_HVGs, 
  HVGs$stimulated_sim_3_20$damage_HVGs, 
  HVGs$stimulated_sim_2_20$damage_HVGs, 
  HVGs$stimulated_sim_1_20$damage_HVGs)

gene_lists <- list(HVGs$control_sim_3_20$control_HVGs, 
  HVGs$control_sim_2_20$control_HVGs, 
  HVGs$control_sim_1_20$control_HVGs, 
  HVGs$stimulated_sim_3_20$control_HVGs, 
  HVGs$stimulated_sim_2_20$control_HVGs, 
  HVGs$stimulated_sim_1_20$control_HVGs
)

# Combine gene lists into one vector and count occurrences
combined_genes <- unlist(gene_lists)
gene_counts <- table(combined_genes)

# Filter for genes that appear in at least three lists
filtered_genes <- names(gene_counts[gene_counts >= 4])
length(filtered_genes)
# Convert gene symbols to Ensembl IDs
gene_conversion <- getBM(attributes = c("external_gene_name", "ensembl_gene_id", "gene_biotype"),
                         filters = "external_gene_name", 
                         values = filtered_genes, 
                         mart = mart)

# Isolate protein-coding genes
protein_coding_genes <- gene_conversion %>%
  filter(gene_biotype == "protein_coding") %>%
  pull(ensembl_gene_id) %>%
  unique()

# Retrieve GO terms and descriptions for the filtered protein-coding gene set
go_terms <- getBM(attributes = c("ensembl_gene_id", "go_id", "name_1006"),  
                  filters = "ensembl_gene_id", 
                  values = protein_coding_genes, 
                  mart = mart)

# Count occurrences of each GO term
go_counts <- go_terms %>%
  group_by(go_id, name_1006) %>%  # Group by GO ID and description
  summarise(count = n(), .groups = 'drop') %>%
  arrange(desc(count)) %>%
  filter(go_id != "", name_1006 != "")

# Get the top 5 GO terms
top_go_terms <- head(go_counts, 10)

# Create bar plot with descriptions
damage_specific <- ggplot(top_go_terms, aes(x = reorder(name_1006, count), y = count)) +  # Use GO term names for x-axis
  geom_bar(stat = "identity", fill = "#D0E7FB") +
  labs(title = "",
       x = "GO Term",
       y = "Count") +
  coord_flip() +  # Flip coordinates for better readability
  theme_classic() + 
  theme(panel.border = element_rect(fill = NA, colour = "white"))

control_specific <- ggplot(top_go_terms, aes(x = reorder(name_1006, count), y = count)) +  # Use GO term names for x-axis
  geom_bar(stat = "identity", fill = "lightgrey") +
  labs(title = "",
       x = "GO Term",
       y = "Count") +
  coord_flip() +  # Flip coordinates for better readability
  theme_classic() + 
  theme(panel.border = element_rect(fill = NA, colour = "white"))


