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
library(lubridate)

# Set the working directory to the Zenodo directory 
# setwd("/Users/name/Zenodo")
setwd("/Users/alicen/Projects/ReviewArticle/Zenodo/")

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
df_yes <- df %>% dplyr::filter(Included == "Yes")

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


# Stacked bar plot with same results
strategy_df$Count_inverse <- 63 - strategy_df$Count

# Order according to increasing count
strategy_df$Strategy <- factor(strategy_df$Strategy, levels = strategy_df$Strategy[order(strategy_df$Count)])


# Performance 
# Reshape & ensure performance is plotted first
data_long <- reshape2::melt(strategy_df, id.vars = "Strategy", measure.vars = c("Count", "Count_inverse"))
data_long$variable <- factor(data_long$variable, levels = c("Count_inverse", "Count"))

# Create the stacked bar graph
strategy_ranks <- ggplot(data_long, aes(x = value, y = Strategy, fill = variable)) +
  geom_bar(stat = "identity", width = 0.55) +
  scale_fill_manual(values = c("Count" = "#001E5C", "Count_inverse" = "#E6E6E6")) +
  theme_classic() + 
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.line = element_blank(),
        legend.position = "none"
  )

ggsave("./A_Search_Strategies/img/methods_bar_chart.png", 
       plot = strategy_ranks, width = 4, height = 3, units = "in", dpi = 300, bg = "transparent")



#-------------------------------------------------------------------------------
# OF MITO-CENTERED STRATEGIES
#-------------------------------------------------------------------------------

# Of those where the method is specified, what are predominant methods? 
# Filter the data frame to include only rows where Included is "Yes"
df_mito <- df %>% dplyr::filter(Threshold != "-")


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


# Create the stacked bar graph
threshold_ranks <- ggplot(mito_df, aes(x = Threshold, y = Count)) +
  geom_bar(stat = "identity", width = 0.55, fill = "#001E5C") +
  scale_y_continuous(breaks = seq(0, 12, 1), limits = c(0, 12)) + 
  theme_classic() + 
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_blank(),
        axis.line = element_blank(),
        legend.position = "none"
  )

ggsave("./A_Search_Strategies/img/thresholds_bar_chart.png", 
       plot = threshold_ranks, width = 14, height = 3.7, units = "in", dpi = 300, bg = "transparent")






#-------------------------------------------------------------------------------
# Values for mitochondrial threshold
#-------------------------------------------------------------------------------

# Test correlation of mitochondrial thresholds with covariates ------

# General stats of threshold values 

df_mito_numbers <- df_mito %>%
  filter(!is.na(as.numeric(as.character(Threshold))))
df_mito_numbers$Explanation <- NULL

# Convert Year and Month into single date 
df_mito_numbers$Month <- match(df_mito_numbers$Month, month.abb)
df_mito_numbers$Date <- make_date(year = df_mito_numbers$Year, month = df_mito_numbers$Month, day = 1)


# Convert the Threshold column to numeric
df_mito_numbers$Threshold <- as.numeric(as.character(df_mito_numbers$Threshold))
median(as.numeric(df_mito_numbers$Threshold))
mean(as.numeric(df_mito_numbers$Threshold))
sd(as.numeric(df_mito_numbers$Threshold))


# Question 1: relationship between time & thresholds ----

# Linear regression model
model_year <- lm(Threshold ~ Date, data = df_mito_numbers)

# Summary of the model
model_summary <- summary(model_year)
r_squared <- model_summary$r.squared
p_value <- model_summary$coefficients[2, 4] 

# Visualization
ggplot(df_mito_numbers, aes(x = Date, y = Threshold)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  labs(x = "", y = "Mitochondrial Threshold (%)") +
  theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12)) +
  annotate("text", 
           x = min(df_mito_numbers$Date), 
           y = max(df_mito_numbers$Threshold), 
           label = paste0("R² = ", round(r_squared, 4), "\n", 
                          "p = ", signif(p_value, 3)), 
           hjust = 0, 
           size = 5, 
           color = "black")

ggsave("./A_Search_Strategies/img/correlation_mito_date.png", 
      width = 7, height = 6.5, units = "in", dpi = 300, bg = "transparent")


# Question 2: Do the thresholds used differ across the species under investigation -----

# Only including species with at least three data points
df_mito_numbers_species <- subset(df_mito_numbers, subset = (Species %in% c("Human","Mouse")))

# Compare humans vs other species 
kruskal_test <- kruskal.test(Threshold ~ Species, data = df_mito_numbers_species)
p_value <- kruskal_test$p.value

# Visualization with p-value
ggplot(df_mito_numbers_species, aes(x = Species, y = Threshold)) +
  geom_boxplot() +
  labs(title = "",
       x = "",
       y = "Mitochondrial Threshold (%)") +
  theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12), 
        legend.position = "none") +
  annotate("text", 
           x = 1.5,  # Adjust this based on your data to position the p-value
           y = max(df_mito_numbers_species$Threshold, na.rm = TRUE), 
           label = paste0("p = ", signif(p_value, 3)), 
           size = 5, 
           color = "black")


ggsave("./A_Search_Strategies/img/correlation_mito_species.png", 
       width = 4, height = 6.6, units = "in", dpi = 300, bg = "transparent")



# Question 3: Do the thresholds used differ across the (broad) study application -----

# Only including application with at least three data points
df_mito_numbers_application <- subset(df_mito_numbers, subset = Application != "Development")

# Normality testing 
shapiro_test <- by(df_mito_numbers_application$Threshold, df_mito_numbers_application$Application, shapiro.test)
# Cancer        p-value 8.337e - 05 (no)
# Cell line     p-value 0.6369      (yes)
# Tissue        p-value 0.0164      (no)
# (ANOVA not appropriate)

# Compare humans vs other species
kruskal_test <- .test(Threshold ~ Application, data = df_mito_numbers_application)
p_value <- kruskal_test$p.value

# Perform pairwise Kruskal-Wallis tests (post-hoc)
kruskal_pairwise <- pairwise.wilcox.test(df_mito_numbers_application$Threshold, 
                                         df_mito_numbers_application$Application, 
                                         p.adjust.method = "BH")  

# Extract adjusted p-values
p_adjusted_values <- kruskal_pairwise$p.value
p_value_cancer_tissue <- p_adjusted_values["Cancer", "Tissue"]


# Visualization with p-value
ggplot(df_mito_numbers_application, aes(x = Application, y = Threshold)) +
  geom_boxplot() +
  labs(title = "",
       x = "",
       y = "Mitochondrial Threshold (%)") +
  theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12), 
        legend.position = "none") +
  annotate("text", 
           x = 1.5,  # Adjust x to position p-value between groups
           y = max(df_mito_numbers_application$Threshold, na.rm = TRUE), 
           label = paste0("p_adj = ", signif(p_value_cancer_tissue, 3)), 
           size = 5, 
           color = "black")


ggsave("./A_Search_Strategies/img/correlation_mito_application.png", 
       width = 5, height = 6.6, units = "in", dpi = 300, bg = "transparent")



### End 
