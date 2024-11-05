# SCRIPT CONTEXT 
#
# Collect results from performance and comparison metrics to provide overall rankings 
# for the tools.
# 
# A. Performance ranked 
#    Score combined from 4 equally weighted components 
#    1 - Precision ranking from all ground truth datasets
#    2 - PR-AUC ranking from all ground truth datasets 
#    3 - HVG similarity (by Jaccard Index) to undamaged control on all simulated data 
#    4 - DEG correctness (by F1 Score) to undamaged control on all simulated data
# 
# B. Consistency ranked 
#    Similarity score combined from ground truth, non-ground truth, and simulated datasets 
# 
# C. Usability score 
#    Manual scoring of methods based on 4 equally weighted objective metrics 
#    and 2 more subjective scores, the first based on our experience with the tool 
#    and the second based on the stars assigned to the GitHub page (with the intention of reflecting other User's experience with the tool)
#    1 - Documentation present (Yes - 1, No - 0) 
#    2 - More than one GitHub/ deployment Upload (Yes - 1, No - 0)
#    3 - Last adjustment to tool less than two years ago (Yes - 1, No - 0)
#    4 - Did the tool require alterations of the core code to run (Yes - 0, No - 1)
#    5 - Ease of use ( * / 5 ) - Score based on many tool-specific factors contributing to ease of use 
#    6 - GitHub stars ( * / 5 ) - 

#-------------------------------------------------------------------------------
# PREPARATIONS
#-------------------------------------------------------------------------------

# Load necessary libraries
packages <- c("cowplot", "dplyr", "ggrepel", "ggplot2", 'irr', "Seurat", "purrr")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    library(pkg)
  }
}



#-------------------------------------------------------------------------------
# PERFORMANCE
#-------------------------------------------------------------------------------

# Performance -----
performance_rankings <- read.csv("/home/alicen/Projects/ReviewArticle/summarised_results/summarised_rankings.csv")

# Sum the rankings 
performance_rankings$Final_Rank <- (performance_rankings$precision + performance_rankings$PR.AUC + performance_rankings$HVG.Jaccard + performance_rankings$DEG.F1)

# Convert the ranking to a score out of 40 where higher score means better performance 
performance_rankings$Final_Rank <- 40 - performance_rankings$Final_Rank


#-------------------------------------------------------------------------------
# Consistency 
#-------------------------------------------------------------------------------

# Consistency -----
consistecy_rankings <- read.csv("/home/alicen/Projects/ReviewArticle/summarised_results/consistency_rankings.csv")

# Sum the rankings 
consistecy_rankings$Final_Rank <- (consistecy_rankings$groundtruth + consistecy_rankings$non_groundtruth + consistecy_rankings$simulated)
                                     
# Convert the ranking to a score out of 40 where higher score means better performance 
consistecy_rankings$Final_Rank <- 40 - consistecy_rankings$Final_Rank

#-------------------------------------------------------------------------------
# Usability 
#-------------------------------------------------------------------------------

# Usability -----
usability_rankings <- read.csv("/home/alicen/Projects/ReviewArticle/summarised_results/usability_scores.csv")

# Apply min-max normalisation for the number of GitHub repo stars 
min_stars <- min(usability_rankings$Stars)
max_stars <- max(usability_rankings$Stars)

# Apply min-max normalization to scale the Stars to a range [0, 5]
usability_rankings$Stars_score <- 5 * (usability_rankings$Stars - min_stars) / (max_stars - min_stars)

# Add together all 6 scores to get final score for usability out of 15 (closer to 15 ~ best)
usability_rankings$Final_Rank <- (usability_rankings$GitHub_Upload + usability_rankings$GitHub_update + 
                                  usability_rankings$Documentation + usability_rankings$Alterations + 
                                  usability_rankings$Ease_of_use + usability_rankings$Stars_score)


#-------------------------------------------------------------------------------
# Plot 
#-------------------------------------------------------------------------------

# Combine all final ranks into one data frame 
Final_Ranks <- performance_rankings[, c("method", "Final_Rank")]
Final_Ranks$Performance <- Final_Ranks$Final_Rank
Final_Ranks$Final_Rank <- NULL
Final_Ranks$Consistency <- consistecy_rankings[, c("Final_Rank")]
Final_Ranks$Usability <- usability_rankings[, c("Final_Rank")]
View(Final_Ranks)

# For creating 'progress bar' like plot, add complementary column for each metric
Final_Ranks$Performance_inverse <- 40 - Final_Ranks$Performance
Final_Ranks$Consistency_inverse <- 40 - Final_Ranks$Consistency
Final_Ranks$Usability_inverse <- 15 - Final_Ranks$Usability

# Create stacked bar plot 

# Order according to performance 
Final_Ranks$method <- factor(Final_Ranks$method, levels = Final_Ranks$method[order(Final_Ranks$Performance)])


# Performance 

# Reshape & ensure performance is plotted first
data_long <- reshape2::melt(Final_Ranks, id.vars = "method", measure.vars = c("Performance", "Performance_inverse"))
data_long$variable <- factor(data_long$variable, levels = c("Performance_inverse", "Performance"))

# Create the stacked bar graph
performance_plot <- ggplot(data_long, aes(x = value, y = method, fill = variable)) +
  geom_bar(stat = "identity", width = 0.55) +
  scale_fill_manual(values = c("Performance" = "#001E5C", "Performance_inverse" = "#E6E6E6")) +
  theme_classic() + 
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.line = element_blank(),
        legend.position = "none"
        )

ggsave(filename = file.path("/home/alicen/Projects/ReviewArticle/summarised_results/performance_plot.png"), 
       plot = performance_plot, width = 7, height = 12, dpi = 300)


# Consistency 

# Reshape & ensure performance is plotted first
data_long <- reshape2::melt(Final_Ranks, id.vars = "method", measure.vars = c("Consistency", "Consistency_inverse"))
data_long$variable <- factor(data_long$variable, levels = c("Consistency_inverse", "Consistency"))

# Create the stacked bar graph
consistency_plot <- ggplot(data_long, aes(x = value, y = method, fill = variable)) +
  geom_bar(stat = "identity", width = 0.55) +
  scale_fill_manual(values = c("Consistency" = "#8DC5BD", "Consistency_inverse" = "#E6E6E6")) +
  theme_classic() + 
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.line = element_blank(),
        legend.position = "none"
  )

ggsave(filename = file.path("/home/alicen/Projects/ReviewArticle/summarised_results/consistency_plot.png"), 
       plot = consistency_plot, width = 7, height = 12, dpi = 300)

# Usability 

# Reshape & ensure performance is plotted first
data_long <- reshape2::melt(Final_Ranks, id.vars = "method", measure.vars = c("Usability", "Usability_inverse"))
data_long$variable <- factor(data_long$variable, levels = c("Usability_inverse", "Usability"))

# Create the stacked bar graph
usability_plot <- ggplot(data_long, aes(x = value, y = method, fill = variable)) +
  geom_bar(stat = "identity", width = 0.55) +
  scale_fill_manual(values = c("Usability" = "#CF9EBB", "Usability_inverse" = "#E6E6E6")) +
  theme_classic() + 
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.line = element_blank(),
        legend.position = "none"
  )

ggsave(filename = file.path("/home/alicen/Projects/ReviewArticle/summarised_results/usability_plot.png"), 
       plot = usability_plot, width = 7, height = 12, dpi = 300)

# Combine all 

# Remove labels for arranging
consistency_plot <- consistency_plot + theme(axis.text = element_blank())
usability_plot <- usability_plot + theme(axis.text = element_blank())

combined_ranking <- performance_plot + consistency_plot + usability_plot 

ggsave(filename = file.path("/home/alicen/Projects/ReviewArticle/summarised_results/usability_plot.png"), 
       plot = combined_ranking, width = 14, height = 12, dpi = 300)


### End 











