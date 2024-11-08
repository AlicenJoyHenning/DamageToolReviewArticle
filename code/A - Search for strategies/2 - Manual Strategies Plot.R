# CONTEXT 
# 
# Following manual review, the protocols were summarized according to:
# 1. Explanation 
#    Whether damaged filtering protocol was explained 
#
# 2. Strategy used 
#    If explained, what strategy was used. First strategy was matched to one of the following 11  
#    but if no appropriate match was found the method was added to 'Other' 
# 
#    Tools: ddqc, DropletQC, ensembleKQC, miQC, scater, valiDrops
#    Manaul (outlier detection of any kind): 
#     i. All
#        Any outlier based filtering based on 5 or more criteria such as feature and UMI counts,
#        mitochondrial, ribosomal and MALAT1 expression/ percentages.
#    ii. Mito -ribo 
#        Any outlier based filtering focusing mitochondrial and ribosomal expression/ 
#        percentages, may include feature and UMI counts.
#   iii. Mito 
#        Any outlier based filtering focusing mitochondrial expression/ 
#        percentages, may include feature and UMI counts.
#    iv. Mito -isolated
#        Any outlier based filtering focusing on only mitochondrial expression/ 
#        percentages and nothing else
#     v. MALAT1
#        Any outlier based filtering focusing on only MALAT1 expression/ 
#        percentages and nothing else


#-------------------------------------------------------------------------------
# Preparations
#-------------------------------------------------------------------------------

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Set the working directory to the Zenodo directory 
# setwd("/Users/name/Zenodo")

#-------------------------------------------------------------------------------
# ALL STRATEGIES 
#-------------------------------------------------------------------------------

# Read the CSV file into a data frame
df <- read.csv("./A_Search_Strategies/data/PubMed_search_papers_reviewed.csv")

# Edit for consistency 
df$Strategy <- ifelse(df$Strategy == "Not specified ", "Not specified", df$Strategy)
df$Strategy <- ifelse(df$Strategy == "Manual_mito_isolated" | df$Strategy == "manual_mito", "Manual_mito", df$Strategy)
df$Strategy <- ifelse(df$Strategy == "SampleQC" | df$Strategy == "scater", "Tool", df$Strategy)


# Colours 
unique(df$Strategy)
Colours <- c("Other" = "#CF9EBB",
             "Not specified" = "#E6E6E6", 
             "Manual_feature_UMI" = "#E6E6E6", 
             "Tool" =  "#8DC5BD",
             "Manual_mito_ribo" = "#8294AC",
             "Manual_mito" = "#001E5C",
             "Yes" = "#001E5C",
             "No" = "#E6E6E6"
             )


# Prepare data frame for plotting
included_counts <- table(df$Included)
included_df <- as.data.frame(included_counts)
colnames(included_df) <- c("Included", "Count")

# Create a pie chart for the strategies used
methods_included_plot <- ggplot(included_df, aes(x = "", y = Count, fill = Included)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = Colours) + 
  coord_polar("y") +
  theme_void() 

# Save (with transparent background)
ggsave("./A_Search_Strategies/img/included_pie_chart.png", 
       plot = methods_included_plot, width = 5, height = 5, units = "in", dpi = 300, bg = "transparent")



#-------------------------------------------------------------------------------
# OF ALL PRESENT STRATEGIES 
#-------------------------------------------------------------------------------

# Of those where the method is specified, what are predominant methods? 
# Filter the data frame to include only rows where Included is "Yes"
df_yes <- df %>% filter(Included == "Yes")

# Prepare data frame for plotting
strategy_counts <- table(df_yes$Strategy)
strategy_df <- as.data.frame(strategy_counts )
colnames(strategy_df) <- c("Strategy", "Count")

# As before, create a pie chart for the Strategy column using the filtered data frame
strategy_counts <- table(df_yes$Strategy)
strategy_df <- as.data.frame(strategy_counts)
colnames(strategy_df) <- c("Strategy", "Count")

methods_prevalence_plot <- ggplot(strategy_df, aes(x = "", y = Count, fill = Strategy)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = Colours) +
  coord_polar("y") +
  theme_void() 

# Save (with transparent background)
ggsave("./A_Search_Strategies/img/methods_pie_chart.png", 
       plot = methods_prevalence_plot, width = 5, height = 5, units = "in", dpi = 300, bg = "transparent")



### End 

#-------------------------------------------------------------------------------
# OF MITO-CENTERED STRATEGIES
#-------------------------------------------------------------------------------

# Of those where the method is specified, what are predominant methods? 
# Filter the data frame to include only rows where Included is "Yes"
df_mito <- df %>% filter(Threshold != "-")

# Prepare data frame for plotting
mito_counts <- table(df_mito$Threshold)
mito_df <- as.data.frame(mito_counts )
colnames(mito_df) <- c("Threshold", "Count")

Colours <- c(
  "3" = "#4C9A92",         
  "5" = "#8DC5BD",        
  "6" = "#A8D8D3",       
  "7.5" = "#D5E7E3",       
  "10" = "#E6E6E6",        
  "15" = "#D1D1D1",        
  "20" = "#B4B4B4",        
  "25" = "#8C8C8C",        
  "30" = "#4D5C7B",        
  "40" = "#3B4A65",        
  "50" = "#2A3A54",        
  "80" = "#1A2A3F",        
  "Not specified" = "#F1D0D6",  
  "Sample-specific" = "#CF9EBB" 
)

# Get the unique Threshold values and order them based on the colour palette
ordered_thresholds <- names(Colours)
mito_df$Threshold <- factor(mito_df$Threshold, levels = ordered_thresholds)
mito_df <- mito_df[order(mito_df$Threshold), ]

# PLOT
mito_plot <- ggplot(mito_df, aes(x = "", y = Count, fill = Threshold)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = Colours) +
  coord_polar("y") +
  theme_void() 

# Save (with transparent background)
ggsave("./A_Search_Strategies/img/mito_pie_chart.png", 
       plot = mito_plot, width = 5, height = 5, units = "in", dpi = 300, bg = "transparent")



### End 
