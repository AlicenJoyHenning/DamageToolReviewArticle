# CONTEXT 
# 
# Following manual review, the protocols were summarized according to:
# 1. Explanation 
#    Whether damaged filtering protocol was explained 
#
# 2. Strategy used 
#    If explained, what strategy was used. First strategy was matched to of the following 11  
#    but if no appropriate match was found the method was added to other 
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



#-------------------------------------------------------------------------------
# Plotting 
#-------------------------------------------------------------------------------

# Read the CSV file into a data frame
df <- read.csv("/home/alicen/Projects/ReviewArticle/article_search/PubMed_search_outcome.csv")

# Colours 
Colours <- c("Yes" = "#DCDBDD",
             "No" = "#011E5C",
             "Manual_feature_UMI" = "#00907F",
             "Manual_mito" = "#DDECE2",
             "Manual_mito_isolated" = "#A6BEAF",
             "Manual_mito_ribo" = "#CFDFB7",
             "Other" = "#011E5C",
             "SampleQC" = "#AFD0F7",
             "scater" = "#88A1BD", 
             "Manual" = "#F4F4F4",
             "Tool-based" = "#011E5C"
             )



# Create a pie c/home/alicen/Projects/ReviewArticle/article_search/PubMed_search_outcome.csv# Create a pie chart for the Included column
included_counts <- table(df$Included)
included_df <- as.data.frame(included_counts)
colnames(included_df) <- c("Included", "Count")


methods_included_plot <- ggplot(included_df, aes(x = "", y = Count, fill = Included)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = Colours) + 
  coord_polar("y") +
  theme_void() 

# Filter the data frame to include only rows where Included is "Yes"
df_yes <- df %>% filter(Included == "Yes")

# Create a pie chart for the Strategy column using the filtered data frame
strategy_counts <- table(df_yes$Strategy)
strategy_df <- as.data.frame(strategy_counts)
colnames(strategy_df) <- c("Strategy", "Count")

ggplot(strategy_df, aes(x = "", y = Count, fill = Strategy)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = Colours) +
  coord_polar("y") +
  theme_void() 

# Summarize the Strategy column to categorize methods as "Manual" or "Tool-based"
df_yes <- df_yes %>%
  mutate(Method = case_when(
    grepl("^Manual_", Strategy, ignore.case = TRUE) ~ "Manual",
    Strategy == "Not_specified" ~ "Not_specified",
    TRUE ~ "Tool-based"
  ))

method_counts <- table(df_yes$Method)
method_df <- as.data.frame(method_counts)
colnames(method_df) <- c("Method", "Count")

ggplot(method_df, aes(x = "", y = Count, fill = Method)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  theme_void() +
  ggtitle("Summarized Strategies Pie Chart")

