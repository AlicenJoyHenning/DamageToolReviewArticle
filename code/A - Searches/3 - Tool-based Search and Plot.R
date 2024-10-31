# Visualise tool statistics 

# Load packages -----

install.packages("rcrossref")
library(rcrossref)          
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(cowplot)

# Prepare for visualisations -----

# Read in the data frame downloaded from https://www.scrna-tools.org/table
scRNAseqtools_df <- read.csv("/home/alicen/Projects/ReviewArticle/article_search/scRNA-tools-dates.csv")

# For those tools with no publications, add 0 
scRNAseqtools_df$Citations <- ifelse(scRNAseqtools_df$Citations == "'-", 0, scRNAseqtools_df$Citations)
scRNAseqtools_df$Citations <- as.numeric(scRNAseqtools_df$Citations)

# Fix the Pub dates

# Function to process the Pub Dates column
process_pub_dates <- function(pub_dates) {
  # Split the dates by comma
  dates <- unlist(strsplit(pub_dates, ", "))
  # Remove NA entries
  dates <- dates[dates != "NA"]
  # Keep only the first date if there are multiple dates
  if (length(dates) > 0) {
    return(dates[1])
  } else {
    return(NA)
  }
}

# Apply the function to the Pub Dates column and create a new column with processed dates
scRNAseqtools_df <- scRNAseqtools_df %>%
  mutate(First_Pub_Date = sapply(Pub.Dates, process_pub_dates))

# Remove older dates and replace NA with unknown 
scRNAseqtools_df$Pub.Dates <- NULL
scRNAseqtools_df$First_Pub_Date[is.na(scRNAseqtools_df$First_Pub_Date)] <- "unknown"

# Remove tools with no citations and no DOI 
scRNAseqtools_df$DOIs[is.na(scRNAseqtools_df$DOIs) | scRNAseqtools_df$DOIs == ""] <- "unknown"
scRNAseqtools_df <- subset(scRNAseqtools_df, scRNAseqtools_df$DOIs != "unknown")


# Manually add unknown dates -----

unknown_dates <- subset(scRNAseqtools_df, scRNAseqtools_df$First_Pub_Date == "unknown")

# Address the multiple entries for a column problem 
keep_first_entry <- function(entry) {
  
  split_entry <- strsplit(entry, ",")[[1]]  # Keep the first DOI in the column 
  return(trimws(split_entry[1]))
  
}
 
unknown_dates$DOIs <- sapply(unknown_dates$DOIs, keep_first_entry)


# Use crossref to extract dates from DOI column
extract_date <- function(doi) {
  
  doi_info <- tryCatch({cr_works(dois = doi)}, error = function(e) {return(NULL)})
  
  if (is.null(doi_info)) {
    
    publication_date <- "unknown"
    
  } 
  
  else 
    
  {
    
    publication_date <- doi_info$data$created
    
  }
  
  # Show user progress
  cat(paste0("Date from ", doi, " obtained\n"))
  
  return(publication_date)
  
}


# Apply the function to each entry in the DOIs column and create a new Date column
unknown_dates$Date <- sapply(unknown_dates$DOIs, extract_date)
unknown_dates <- subset(unknown_dates, unknown_dates$Date != "unknown")
unknown_dates$First_Pub_Date <- unknown_dates$Date

# Merge the data frames based on the 'Name' column
merged_df <- merge(scRNAseqtools_df, unknown_dates, by = "Name", all.x = TRUE)

# Replace the First_Pub_Date in scRNAseqtools_df with Date from unknown_dates where First_Pub_Date is "unknown"
merged_df$Final_Date <- ifelse(merged_df$First_Pub_Date.x == "unknown", merged_df$First_Pub_Date.y, merged_df$First_Pub_Date.x)

# Create the final data frame
scRNAseqtools_dff <- data.frame(
  Name = merged_df$Name
 # Date = merged_df$Final_Date
 # Citations = merged_df$Citations.x,
 # Categories = merged_df$Categories.x
)

scRNAseqtools_dff$Date <- merged_df$Final_Date
scRNAseqtools_dff$DOI <- merged_df$DOIs.x
scRNAseqtools_dff$Citations <- merged_df$Citations.x
scRNAseqtools_dff$Categories <- merged_df$Categories.x


# Convert list columns to character columns
scRNAseqtools_dff <- data.frame(lapply(scRNAseqtools_dff, function(x) {
  if (is.list(x)) {
    return(sapply(x, toString))
  } else {
    return(x)
  }
}), stringsAsFactors = FALSE)


scRNAseqtools_dff$Date <- ifelse(scRNAseqtools_dff$Date == "", "unknown", scRNAseqtools_dff$Date)

# Write out to manually enter in unknown 
write.csv2(scRNAseqtools_dff, "C:/Users/alice/OneDrive/Documents/University/Masters/limiric_article/scRNAseq-tools-unknown-to-be-filled-in.csv")
scRNAseqtools_dff <- read.csv2("C:/Users/alice/OneDrive/Documents/University/Masters/limiric_article/scRNAseq-tools-unknown-to-be-filled-in.csv", sep = ";")

scRNAseqtools_dff$DOI <- NULL
scRNAseqtools_dff$X <- NULL

# Add -01 to all entries in the Date column that don't have a day
scRNAseqtools_dff$Date <- ifelse(grepl("^\\d{4}-\\d{2}$", scRNAseqtools_dff$Date), 
                                 paste0(scRNAseqtools_dff$Date, "-01"), 
                                 scRNAseqtools_dff$Date)


# Mark categories as either QC related, QC unrelated, QC focus ------

scRNAseqtools_dff$Colour <- ifelse(
  
  # QC focus if QC is present and there are two or fewer other categories 
  grepl("Quality Control", scRNAseqtools_dff$Categories) & sapply(strsplit(scRNAseqtools_dff$Categories, ", "), length) <= 3, 
  "QC-focus", 
  
  # QC related if QC is present but there are three or more other categories 
  ifelse(grepl("Quality Control", scRNAseqtools_dff$Categories) & sapply(strsplit(scRNAseqtools_dff$Categories, ", "), length) > 3, 
         "QC-related", 
  
  # QC unrelated if QC not present 
         "QC-unrelated")
)


# Damaged cell detection -----

# Manually edit to see which of the QC-focus tools can detect damaged cells 

test <- subset(scRNAseqtools_dff, scRNAseqtools_dff$Colour == "QC-focus")
test$Categories <- NULL
# test$Colour <- NULL

write.csv(test,
          quote = FALSE,
          row.names = FALSE,
          "C:/Users/alice/OneDrive/Documents/University/Masters/limiric_article/complete_details.csv"
)

damaged_focus <- read.csv("C:/Users/alice/OneDrive/Documents/University/Masters/limiric_article/complete_details.csv", sep = ",", header = TRUE)
damaged_focus  <- read.csv("C:/Users/alice/OneDrive/Documents/University/Masters/limiric_article/complete_details.csv", sep = ",", header = FALSE)
colnames(damaged_focus) <- damaged_focus[1,]
damaged_focus <- damaged_focus[-1,]
damaged_focus$Date <- NULL
damaged_focus$Citations <- NULL
damaged_focus$Description <- NULL
damaged_focus <- damaged_focus[, -3]
damaged_focus <- damaged_focus[rowSums(is.na(damaged_focus)) != ncol(damaged_focus), ]


write.csv(damaged_focus,
          quote = FALSE,
          row.names = FALSE,
          "C:/Users/alice/OneDrive/Documents/University/Masters/limiric_article/damage_focus.csv"
)

damaged_focus  <- read.csv("C:/Users/alice/OneDrive/Documents/University/Masters/limiric_article/damage_focus.csv", sep = ",", header = FALSE)
colnames(damaged_focus)[colnames(damaged_focus) == "V1"] <- "Name"
colnames(damaged_focus)[colnames(damaged_focus) == "V2"] <- "Damaged"


# Problem with naming, spacesnullfile()# Problem with naming, spaces 
# look <- as.data.frame(unique(damaged_focus$Damaged_relatedness))

# FIXING
# Trim whitespace from the Damaged_relatedness column
# damaged_focus$Damaged_relatedness <- trimws(damaged_focus$Damaged_relatedness)
# 
# # Standardize the values in the Damaged_relatedness column
# damaged_focus$Damaged_relatedness <- case_when(
#   damaged_focus$Damaged_relatedness %in% c("damage-unrelated", "damage-unrelated ") ~ "damage-unrelated",
#   damaged_focus$Damaged_relatedness %in% c("damage-related", "damage-related ") ~ "damage-related",
#   damaged_focus$Damaged_relatedness %in% c("damage-focus", "damage-focus ") ~ "damage-focus",
#   TRUE ~ damaged_focus$Damaged_relatedness
# )
# 
# # Verify the unique values in the Damaged_relatedness column
# unique(damaged_focus$Damaged_relatedness)
# 
# Save corrected
# write.csv(damaged_focus,
#           quote = FALSE,
#           row.names = FALSE,
#           "C:/Users/alice/OneDrive/Documents/University/Masters/limiric_article/complete_details.csv"
# )

# Plot! Pies ------

# One pie
big_pie_df <- merge(scRNAseqtools_dff, damaged_focus, by = "Name", all.x = TRUE)
big_pie_df$Colour_BP <- ifelse(is.na(big_pie_df$Damaged_relatedness), big_pie_df$Colour, big_pie_df$Damaged_relatedness)
unique(big_pie_df$Colour_BP)

big_pie <- big_pie_df %>%
  group_by(Colour_BP) %>%
  summarize(Count = n())

big_pie$Colour_BP <- ifelse(big_pie$Colour_BP == "damage-unrelated", "QC-focus_damaged-unrelated", big_pie$Colour_BP)
big_pie$Colour_BP <- ifelse(big_pie$Colour_BP == "damage-related", "QC-focus_damaged-related", big_pie$Colour_BP)
big_pie$Colour_BP <- ifelse(big_pie$Colour_BP == "damage-focus", "QC-focus_damaged-focus", big_pie$Colour_BP)



# Define the colors for each category
unique(big_pie$Colour_BP)
color_mapping <- c( 
                   "QC-unrelated" = "#CBCBCB",
                   "QC-related"   = "#AFADF5",
                   "QC-focus_damaged-unrelated" = "#201BE5", 
                   "QC-focus_damaged-related" = "#92DAB3",
                   "QC-focus_damaged-focus" = "#3CB371"
                   )

# QC Pie Chart
QC_pie <- ggplot(big_pie, aes(x = "", y = Count, fill = Colour_BP)) +
  geom_bar(width = 2, stat = "identity", color = "black") +  # Adjust the width to decrease the diameter
  coord_polar(theta = "y", start = 0, direction = 0.7) +
  scale_fill_manual(values = color_mapping) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Category") +
  theme_void() + 
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    legend.position = "none",
    #legend.position = c(1.2, 0.95),  # Position the legend in the top right
    legend.justification = c(1, 1),  # Anchor the legend to the top right
    plot.margin = margin(10, 10, 10, 10)  # Adjust margins to control the plot area
  )

QC_pie

# Pie chart to show that very few QC tools 

broadQC_pie_df <- scRNAseqtools_dff %>%
  group_by(Colour) %>%
  summarize(Count = n())

# Define the colors for each category
color_mapping <- c("QC-related"   = "#AFADF5", 
                   "QC-focus"     = "#201BE5", 
                   "QC-unrelated" = "#CBCBCB")

# QC Pie Chart
QC_pie <- ggplot(broadQC_pie_df, aes(x = "", y = Count, fill = Colour)) +
  geom_bar(width = 2, stat = "identity") +  # Adjust the width to decrease the diameter
  coord_polar(theta = "y", start = 0, direction = 0.7) +
  scale_fill_manual(values = color_mapping) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Category") +
  theme_void() + 
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
     legend.position = "none",
   #legend.position = c(1.2, 0.95),  # Position the legend in the top right
    legend.justification = c(1, 1),  # Anchor the legend to the top right
    plot.margin = margin(10, 10, 10, 10)  # Adjust margins to control the plot area
  )

QC_pie

# Define the colors for each category
damaged_focus <- damaged_focus[-1,]

narrowQC_pie_df <- damaged_focus %>% 
  group_by(Damaged) %>%
  summarize(Count = n())


color_mapping <- c("damage-related"   = "#92DAB3", 
                   "damage-focus"     = "#3CB371",
                   "damage-unrelated" = "#E1E1E1")

# Damage Pie Chart
damage_pie <- ggplot(narrowQC_pie_df, aes(x = "", y = Count, fill = Damaged)) +
  geom_bar(width = 2, stat = "identity") +  
  coord_polar(theta = "y", start = 0, direction = 0.7) +
  scale_fill_manual(values = color_mapping) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Category") +
  theme_void() + 
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    legend.position = "none",
    #legend.position = c(1.2, 0.95),  # Position the legend in the top right
    legend.justification = c(1, 1),  # Anchor the legend to the top right
    plot.margin = margin(10, 10, 10, 10)  # Adjust margins to control the plot area
  )

damage_pie

plot <- QC_pie / damage_pie

# Add a white background to the combined plot
plot <- ggdraw(plot) + theme(plot.background = element_rect(fill = "white", color = NA))


ggsave("C:/Users/alice/OneDrive/Documents/University/Masters/limiric_article/pie_legend.png", plot = plot, width = 10, height = 12, units = "in", dpi = 300)



# Plot! Scatter -----

# Ensure the Date column is of Date type
scatter_df <- merge(scRNAseqtools_dff, damaged_focus, by = "Name", all.x = TRUE)

# Rename & clean our columns
scatter_df$Colour_Final <- ifelse(is.na(scatter_df$Damaged), scatter_df$Colour, scatter_df$Damaged)

scatter_df$Date <- as.Date(scatter_df$Date)
unique(scatter_df$Colour_Final)

# Subset 
scatter_df_subset <- subset(scatter_df, scatter_df$Colour_Final != "QC-unrelated")

# Define the colors for each category
color_mapping <- c("QC-related"       = "#E1E1E1",  
                   "damage-unrelated" = "#AFADF5", 
                   "damage-related"   = "#201BE5", 
                   "damage-focus"     = "#3CB371")

# Define custom breaks for the x-axis
year_2010 <- as.Date("2010-01-01")
year_2012 <- as.Date("2012-01-01")
year_2014 <- as.Date("2014-01-01")
year_2016 <- as.Date("2016-01-01")
year_2018 <- as.Date("2018-01-01")
year_2020 <- as.Date("2020-01-01")
year_2022 <- as.Date("2022-01-01")
year_2024 <- as.Date("2024-01-01")

# Add a small constant to the Citations column to avoid log(0) issues
scatter_df_subset$Citations <- scatter_df_subset$Citations + 1

colnames(scatter_df_subset)

scatter_df$popular <- ifelse(scatter_df$Citations >= 30000, "popular", "average")
scatter_df$popular <- ifelse(scatter_df$Citations <= 30000 & scatter_df$Citations >= 5000, "high", scatter_df$popular)
scatter_df$popular <- ifelse(scatter_df$Name  %in% c("alevin-fry", "CellChat", "scvi-tools", "CellRanger", "Harmony"), "high", scatter_df$popular)
scatter_df$popular <- ifelse(scatter_df$Name  %in% c("cell2location","BEANIE", "SingleR", "CellChat"), "average", scatter_df$popular)
scatter_df <- subset(scatter_df, scatter_df$Name != "scQA") 
scatter_df_subset$Colour_Final <- ifelse(scatter_df_subset$Name  %in% c("valiDrops", "ddqc"), "damage-focus", scatter_df_subset$Colour_Final)
scatter_df_subset$Colour_Final <- ifelse(scatter_df_subset$Name  %in% c("EnsembleKQC"), "damage-related", scatter_df_subset$Colour_Final)
scatter_df_subset$Colour_Final <- ifelse(scatter_df_subset$Name  %in% c("SoCube", "vaeda"), "damage-unrelated", scatter_df_subset$Colour_Final)


scatter <- ggplot(scatter_df, aes(x = Date, y = Citations, colour = popular)) +
  geom_point(data = subset(scatter_df, Name != "alevin-fry"), size = 3, alpha = 0.7) +
  geom_point(data = subset(scatter_df, Name == "alevin-fry"), size = 3) +  # Bring "alevin-fry" to the front
  scale_color_manual(values = c("popular" = "#201BE5", 
                                "high"    = "#AFADF5",
                                "average" = "#E2E2E2")) +
  scale_y_continuous(breaks = c(0, 10000, 20000, 30000)) +
  scale_x_date(breaks = c(year_2010, year_2012, year_2014, year_2016, year_2018, year_2020, year_2022, year_2024), 
               labels = c("2010", "2012", "2014", "2016", "2018", "2020", "2022", "2024")) +
  labs(y = "Citations") +
  theme_classic() + 
  theme(
    axis.title.y = element_text(vjust = 3, face = "bold", size = 10),
    axis.text.x = element_text(hjust = 0),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 10),
    legend.position = "none",  # Position the legend inside the plot on the top right
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.title = element_blank()
  ) + 
  geom_text_repel(data = subset(scatter_df, popular %in%  c("popular", "high")), 
                  aes(label = Name), 
                  size = 3, 
                  color = "black",
                  nudge_y = 0.2,  # Adjust the vertical position of the labels
                  nudge_x = 0.2,
                  segment.color = "grey",  # Color of the line pointing to the point
                  segment.size = 0.5)  # Thickness of the line pointing to the point

scatter

ggsave("C:/Users/alice/OneDrive/Documents/University/Masters/limiric_article/scRNAseq-tools_overall.png", plot = scatter, width = 8, height = 5, units = "in", dpi = 300)



# Create the scatter plot
scatter <- ggplot(scatter_df_subset, aes(x = Date, y = Citations, color = Colour_Final)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = color_mapping) +
  scale_y_log10() +
  scale_x_date(breaks = c(year_2016, year_2018, year_2020, year_2022, year_2024), labels = c("2016", "2018", "2020", "2022", "2024")) +
  labs(y = "Citations") +
  theme_classic() + 
  theme(
    axis.title.y = element_text(vjust = 3, face = "bold", size = 10),
    axis.text.x = element_text(hjust = 0),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.position = c(0.95, 0.95),  # Position the legend inside the plot on the top right
    legend.justification = c(0.8, 1),  # Anchor the legend to the top right corner
    legend.background = element_rect(color = "black", fill = alpha('white', 0.6), size = 0.5),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.title = element_blank()
  ) + 
  geom_text_repel(data = subset(scatter_df_subset, Colour_Final %in% c("damage-related","damage-focus")), 
                  aes(label = Name), 
                  size = 3, 
                  color = "black",
                  nudge_y = 0.2,  # Adjust the vertical position of the labels
                  nudge_x = 0.2,
                  segment.color = "grey",  # Color of the line pointing to the point
                  segment.size = 0.5)  # Thickness of the line pointing to the point


scatter

ggsave("C:/Users/alice/OneDrive/Documents/University/Masters/limiric_article/scRNAseq-tools_subset.png", plot = scatter, width = 9, height = 6, units = "in", dpi = 300)



# Plot! Combine -----

# Combine the plots with the desired layout and labels
combined_plot <- plot_grid(
  plot_grid(QC_pie, damage_pie, ncol = 1, labels = c("a", "b"), rel_heights = c(1, 1)),
  scatter, 
  ncol = 2, 
  rel_widths = c(1, 2), 
  labels = c("", "c")
)

# Add a white background to the combined plot
combined_plot <- ggdraw(combined_plot) + theme(plot.background = element_rect(fill = "white", color = NA))

# Display the combined plot
print(combined_plot)

ggsave("C:/Users/alice/OneDrive/Documents/University/Masters/limiric_article/scRNAseq-tools_plot_version_two_pies.png", plot = combined_plot, width = 15, height = 5, units = "in", dpi = 300)
