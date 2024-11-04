# SCRIPT CONTEXT 
# 
# Scripts to generate plots from the simulated data.
#
# 1. HVG similarity for damaged cell strategies 
# Bar plots showing the similarity (jaccard index) of the HVGs 
# of the filtered cases compared to those of the control (undamaged) cases. 
#  
# 2. DEG correctness for damaged cell strategies 
# Bar plots showing the correctness (F1 score) of celltype-specific DEGs 
# of the filtered cases compared to those of the control (undamaged) cases
#
# 3. Baseline HVG analysis 
# Control (undamaged) vs test (damaged added but unfiltered) cases: Overlap analysis to find mean intersection and isolated cases from 
# which the GO categories related to the isolated cases can be found.


#-------------------------------------------------------------------------------
# PREPARATIONS
#-------------------------------------------------------------------------------

packages <- c("biomaRt", "clusterProfiler", "dplyr", "ggplot2", "glmGamPoi", "Matrix", "robustbase",
              "png", "Seurat", "tidyr", "SingleCellExperiment", "scuttle", "presto")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    library(pkg)
  }
}


#-------------------------------------------------------------------------------
# HVG for filtered cases 
#-------------------------------------------------------------------------------

# 1. HVG similarity -----

HVG_jaccard <- read.csv("/home/alicen/Projects/ReviewArticle/damage_perturbation/analysis_results/hvg_similarity.csv")
View(HVG_jaccard)

# Rank the tools by increasing similarity ----

# Identify columns that start with 'median_' and do not end in 'proportion'
median_columns <- grep("^median_\\d+(\\.\\d+)?$", colnames(HVG_jaccard), value = TRUE)

# Rank each method within each median column and store results in a new data frame
HVG_jaccard_ranked <- HVG_jaccard %>%
  select(method, all_of(median_columns)) %>%  # Keep only method and median columns
  mutate(across(all_of(median_columns), ~ rank(-., ties.method = "min"), .names = "rank_{.col}")) %>%
  select(method, starts_with("rank_"))

# Rename columns to match the specified output format
HVG_jaccard_ranked <- HVG_jaccard_ranked %>%
  rename_with(~ gsub("rank_median_", "median_", .x), starts_with("rank_"))

# Add final ranking by summing across datasets (lowest rank is the best performing)
HVG_jaccard_ranked$Final_rank <- HVG_jaccard_ranked$median_2.5 + HVG_jaccard_ranked$median_5 + HVG_jaccard_ranked$median_10 + HVG_jaccard_ranked$median_15 + HVG_jaccard_ranked$median_20

View(HVG_jaccard_ranked)

# Save and view final results
write.csv(HVG_jaccard_ranked, 
          file = "/home/alicen/Projects/ReviewArticle/damage_perturbation/analysis_results/HVG_similarity_ranked.csv",
          quote = FALSE, 
          row.names = FALSE
)


# Plot the similarity bars ----
# Plot the similarity bars with confidence intervals ----
# Define the colors
strategy_colours <- c(
  "ddqc" = "#CE9DBA",
  "DropletQC" = "#A799C9",
  "ensembleKQC" = "#E7E4F6", 
  "miQC" = "#808C98", 
  "SampleQC" = "#CED5DB",
  "scater" = "#88A1BD", 
  "valiDrops" = "#D3E2F6",
  "manual_all" = "#4F618F",
  "manual_mito_isolated" = "#A6BEAE", 
  "manual_mito" = "#DCECE2", 
  "manual_mito_ribo" = "#9DBD78",
  "manual_malat1" = "#D7E7BE"
)

# Define plot theme with error bars
plot_similarity <- function(damage, data, strategy_colours) {
  
  # Extract DropletQC value for the current damage percentage
  dropletqc_value <- data %>%
    filter(method == "DropletQC" & damage_percent == damage) %>%
    pull(median_score)
  
  # Create the plot with error bars
  p <- ggplot(data %>% filter(damage_percent == damage), 
              aes(y = method, x = median_score, fill = method)) +
    geom_bar(stat = "identity") +
    geom_errorbarh(aes(xmin = lowest_score, xmax = highest_score), height = 0.3, color = "gray40") +
    geom_vline(xintercept = dropletqc_value, color = "black", linetype = "dashed", linewidth = 0.5) +
    scale_fill_manual(values = strategy_colours) +
    coord_cartesian(xlim = c(0.65, 1)) +
    labs(title = "",
         x = "Jaccard Index",
         y = "") +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 16, vjust = -0.5),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_text(size = 16, vjust = -1.5),
          panel.border = element_rect(fill = NA, color = "black"))
  
  return(p)
}

# Define the order of strategies
strategy_order <- c(
  "ddqc", "DropletQC", "ensembleKQC", "miQC", "SampleQC", "scater", "valiDrops",
  "manual_all", "manual_mito_isolated", "manual_mito", "manual_mito_ribo", "manual_malat1"
)

# Reshape and order for plotting, including lowest and highest columns for error bars
HVG_jaccard <- HVG_jaccard %>%
  pivot_longer(cols = starts_with("median_"), names_to = "damage_percent", values_to = "median_score") %>%
  mutate(
    damage_percent = as.numeric(sub("median_", "", damage_percent)),
    method = factor(method, levels = strategy_order)  # Ensure the method column has the specified order
  ) %>%
  # Add lowest and highest scores for error bars
  mutate(
    lowest_score = case_when(damage_percent == 2.5 ~ lowest_2.5,
                             damage_percent == 5 ~ lowest_5,
                             damage_percent == 10 ~ lowest_10,
                             damage_percent == 15 ~ lowest_15,
                             damage_percent == 20 ~ lowest_20),
    highest_score = case_when(damage_percent == 2.5 ~ highest_2.5,
                              damage_percent == 5 ~ highest_5,
                              damage_percent == 10 ~ highest_10,
                              damage_percent == 15 ~ highest_15,
                              damage_percent == 20 ~ highest_20)
  )

# Create the plots using the function
plots <- lapply(unique(HVG_jaccard$damage_percent), function(damage) {
  plot_similarity(damage, HVG_jaccard, strategy_colours)
})

# Name the plots
names(plots) <- paste0("damage_", unique(HVG_jaccard$damage_percent))

# Compile and save plots
HVG_similarity_plot <- plots$damage_2.5 | plots$damage_5 | plots$damage_10 | plots$damage_15 | plots$damage_20

ggsave(filename = "/home/alicen/Projects/ReviewArticle/damage_perturbation/analysis_results/HVG_similarity_plot.png",
       plot = HVG_similarity_plot,
       width = 17, height = 6.4, units = "in")



#-------------------------------------------------------------------------------
# DEG correctness of filtered cases
#-------------------------------------------------------------------------------

# 2. DEG correctness -----
DEG_F1 <- read.csv("/home/alicen/Projects/ReviewArticle/damage_perturbation/analysis_results/deg_correctness.csv")
View(DEG_F1)

# Define the order of strategies
strategy_order <- c(
  "ddqc", "DropletQC", "ensembleKQC", "miQC", "SampleQC", "scater", "valiDrops",
  "manual_all", "manual_mito_isolated", "manual_mito", "manual_mito_ribo", "manual_malat1"
)

# Reshape and order for plotting, including minimum and maximum columns for error bars
DEG_F1 <- DEG_F1 %>%
  pivot_longer(cols = starts_with("median_F1_"), names_to = "damage_percent", values_to = "median_F1") %>%
  mutate(
    damage_percent = as.numeric(sub("median_F1_", "", damage_percent)),  # Extract numeric damage percentage
    method = factor(method, levels = strategy_order),  # Ensure the method column has the specified order
    # Add minimum and maximum scores for error bars
    min_F1 = case_when(damage_percent == 2.5 ~ min_F1_2.5,
                       damage_percent == 5 ~ min_F1_5,
                       damage_percent == 10 ~ min_F1_10,
                       damage_percent == 15 ~ min_F1_15,
                       damage_percent == 20 ~ min_F1_20),
    max_F1 = case_when(damage_percent == 2.5 ~ max_F1_2.5,
                       damage_percent == 5 ~ max_F1_5,
                       damage_percent == 10 ~ max_F1_10,
                       damage_percent == 15 ~ max_F1_15,
                       damage_percent == 20 ~ max_F1_20)
  )

# Define plot theme with error bars for DEG_F1
plot_similarity <- function(damage, data, strategy_colours) {
  
  # Extract DropletQC value for the current damage percentage
  dropletqc_value <- data %>%
    filter(method == "DropletQC" & damage_percent == damage) %>%
    pull(median_F1)  # Pull median F1 score for DropletQC
  
  # Create the plot with error bars
  p <- ggplot(data %>% filter(damage_percent == damage), 
              aes(y = method, x = median_F1, fill = method)) +
    geom_bar(stat = "identity") +
    geom_errorbarh(aes(xmin = min_F1, xmax = max_F1), height = 0.3, color = "gray40") +  # Use min/max F1 values
    geom_vline(xintercept = dropletqc_value, color = "black", linetype = "dashed", linewidth = 0.5) +
    scale_fill_manual(values = strategy_colours) +
    coord_cartesian(xlim = c(0.65, 1)) +
    labs(title = "",
         x = "F1 Score",
         y = "") +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 16, vjust = -0.5),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_text(size = 16, vjust = -1.5),
          panel.border = element_rect(fill = NA, color = "black"))
  
  return(p)
}


plots <- lapply(unique(DEG_F1$damage_percent), function(damage) {
  plot_similarity(damage, DEG_F1, strategy_colours)
})

# Name the plots
names(plots) <- paste0("damage_", unique(DEG_F1$damage_percent))

# Compile and save plots
DEG_F1_plot  <- plots$damage_2.5 | plots$damage_5 | plots$damage_10 | plots$damage_15 | plots$damage_20

ggsave(filename = "/home/alicen/Projects/ReviewArticle/damage_perturbation/analysis_results/DEG_correctness_plot.png",
       plot = DEG_F1_plot ,
       width = 17, height = 6.4, units = "in")


#-------------------------------------------------------------------------------
# BASELINE HVG
#-------------------------------------------------------------------------------

# 3. Baseline HVG ----

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
filtered_genes <- names(gene_counts[gene_counts >= 2])

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



### End 

