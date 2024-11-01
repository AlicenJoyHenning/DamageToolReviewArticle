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


#-------------------------------------------------------------------------------
# HVG for filtered cases 
#-------------------------------------------------------------------------------

# 2. HVG similarity -----

HVG_jaccard <- read.csv("/home/alicen/Projects/ReviewArticle/damage_perturbation/analysis_results/hvg_similarity_5000.csv")

# Rank the tools by increasing similarity ----

# Function to rank and assign points for similarity 
rank_and_assign_points <- function(df, metric) {
  df %>%
    arrange(desc(!!sym(metric))) %>%
    mutate(rank = row_number()) %>%
    mutate(points = 11- rank) %>%
    select(method, points)
}

# Perform ranking
results_ordered <- HVG_jaccard[, c("method", "damage_percent", "median_score")]
ranked_results <- results_ordered %>% group_by(damage_percent)
similarity_ranked <- rank_and_assign_points(ranked_results, "median_score") %>% 
  group_by(method) %>% 
  summarise(similarity = sum(points))


# Save and view final results
write.csv(similarity_ranked, 
          file = "/home/alicen/Projects/ReviewArticle/damage_perturbation/analysis_results/HVG_similarity_ranked.csv",
          quote = FALSE, 
          row.names = FALSE
)

View(similarity_ranked)

# Plot the similarity bars ----

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

# Define plot theme
plot_similarity <- function(damage, data, strategy_colours) {
  # Extract DropletQC value for the current damage percentage
  dropletqc_value <- data %>%
    filter(method == "DropletQC" & damage_percent == damage) %>%
    pull(median_score)
  
  # Create the plot
  p <- ggplot(data %>% filter(damage_percent == damage), aes(y = method, x = median_score, fill = method)) +
    geom_bar(stat = "identity") +
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

# Reshape and order for plotting
HVG_jaccard <- HVG_jaccard %>%
  pivot_longer(cols = starts_with("median_"), names_to = "damage_percent", values_to = "median_score") %>%
  mutate(
    damage_percent = as.numeric(sub("median_", "", damage_percent)),
    method = factor(method, levels = strategy_order)  # Ensure the method column has the specified order
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
         width = 17, height = 6, units = "in")


#-------------------------------------------------------------------------------
# DEG correctness of filtered cases
#-------------------------------------------------------------------------------

# 3. DEG correctness -----

negative_controls <- read.csv("/home/alicen/Projects/ReviewArticle/damage_perturbation/analysis_results/DEG_negative_controls.csv")
DEG_F1 <-  read.csv("/home/alicen/Projects/ReviewArticle/damage_perturbation/analysis_results/deg_correctness.csv")

# Function to rank and assign points for similarity
rank_and_assign_points <- function(df, metric) {
  df %>%
    arrange(desc(!!sym(metric))) %>%
    mutate(rank = row_number()) %>%
    mutate(points = 11 - rank) %>%
    select(method, points)
}

# Initialize an empty dataframe to store the points
points_df <- data.frame(method = DEG_F1$method)

# List of damage columns
damage_columns <- c("damaged_2.5", "damaged_5", "damaged_10", "damaged_15", "damaged_20")

# Loop through each damage column, rank, and assign points
for (damage in damage_columns) {
  ranked <- rank_and_assign_points(DEG_F1, damage)
  points_df <- points_df %>%
    left_join(ranked, by = "method", suffix = c("", paste0("_", damage)))
}

# Sum the points across all damage columns
DEG_correctness_ranked <- points_df %>%
  rowwise() %>%
  mutate(total_points = sum(c_across(starts_with("points_")), na.rm = TRUE)) %>%
  select(method, total_points)

# Save and view final results
write.csv(DEG_correctness_ranked, 
          file = "/home/alicen/Projects/ReviewArticle/damage_perturbation/analysis_results/DEG_correctness_ranked.csv",
          quote = FALSE, 
          row.names = FALSE)
 
View(DEG_correctness_ranked)

# Plot the correctness bars -----

# Define plot theme
plot_similarity <- function(damage, data, strategy_colours) {
  # Extract DropletQC value for the current damage percentage
  dropletqc_value <- data %>%
    filter(method == "DropletQC") %>%
    pull(!!sym(damage))
  
  # Create the plot
  p <- ggplot(data, aes(y = method, x = !!sym(damage), fill = method)) +
    geom_bar(stat = "identity") +
    geom_vline(xintercept = dropletqc_value, color = "black", linetype = "dashed", linewidth = 0.5) +
    scale_fill_manual(values = strategy_colours) +
    coord_cartesian(xlim = c(0.4, 1)) +  # Set x-axis range from 0.65 to 1
    labs(title = "",
         x = "F1 score",
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

# Ensure the method column has the specified order
DEG_F1 <- DEG_F1 %>%
  mutate(method = factor(method, levels = strategy_order))

# List of damage columns
damage_columns <- c("damaged_2.5", "damaged_5", "damaged_10", "damaged_15", "damaged_20")

# Create the plots using the function
plots <- lapply(damage_columns, function(damage) {
  plot_similarity(damage, DEG_F1, strategy_colours)
})

# Name the plots
names(plots) <- damage_columns

# Compile and save plots
DEG_correctness_plot <- plots$damaged_2.5 | plots$damaged_5 | plots$damaged_10 | plots$damaged_15 | plots$damaged_20

ggsave(filename = "/home/alicen/Projects/ReviewArticle/damage_perturbation/analysis_results/DEG_correctness_plot.png",
       plot = DEG_correctness_plot,
       width = 17, height = 6, units = "in")

### End 

