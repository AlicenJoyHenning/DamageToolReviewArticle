# CONTEXT 
# 
# Using scrna-tools table, https://www.scrna-tools.org/tools, explore the tools 
# available that are related to cell quality control, specifically damaged cell detection. 

#-------------------------------------------------------------------------------
# Libraries 
#-------------------------------------------------------------------------------

# Load packages -----

# install.packages("rcrossref")
library(rcrossref)          
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(cowplot)

# Set the working directory to the Zenodo directory 
# setwd("/Users/name/Zenodo")

#-------------------------------------------------------------------------------
# Data cleaning and manual editing 
#-------------------------------------------------------------------------------

# Prepare for visualisations -----

# Read in the data frame downloaded from https://www.scrna-tools.org/table
scRNAseqtools_df <- read.csv("./A_Search_Strategies/data/scRNA-tools-dates.csv")

# For those tools with no publications, add 0 
scRNAseqtools_df$Citations <- ifelse(scRNAseqtools_df$Citations == "'-", 0, scRNAseqtools_df$Citations)
scRNAseqtools_df$Citations <- as.numeric(scRNAseqtools_df$Citations)

# Function to process the PubDates column
process_pub_dates <- function(pub_dates) {
  
  # Remove NA
  dates <- unlist(strsplit(pub_dates, ", "))
  dates <- dates[dates != "NA"]
  
  # Keep only the first date if there are multiple dates
  if (length(dates) > 0) {
    return(dates[1])
  } else {
    return(NA)
  }
  
}

# Apply the function to the PubDates column and create a new column with processed dates
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
write.csv2(scRNAseqtools_dff, "./A_Search_Strategies/data/scRNAseq-tools-unknown-to-be-filled-in.csv")

#-------------------------------------------------------------------------------
# Plot
#-------------------------------------------------------------------------------

scRNAseqtools_dff <- read.csv("./A_Search_Strategies/data/scRNAseq-tools-unknown-to-be-filled-in.csv", sep = ";", header = TRUE)
scRNAseqtools_dff$X <- NULL
damaged_focus  <- read.csv("./A_Search_Strategies/data/scRNA-tools-damage-focus.csv", sep = ",", header = TRUE)

# Plot Scatter -----

# Define custom breaks for the x-axis
year_2010 <- as.Date("2010-01-01")
year_2012 <- as.Date("2012-01-01")
year_2014 <- as.Date("2014-01-01")
year_2016 <- as.Date("2016-01-01")
year_2018 <- as.Date("2018-01-01")
year_2020 <- as.Date("2020-01-01")
year_2022 <- as.Date("2022-01-01")
year_2024 <- as.Date("2024-01-01")

# Merge dataframes
scatter_df <- merge(scRNAseqtools_dff, damaged_focus, by = "Name", all.x = TRUE)

# Ensure the Date column is of Date type (1 of month given if no date specified)
scatter_df$Date <- as.Date(ifelse(nchar(scatter_df$Date) == 7, paste0(scatter_df$Date, "-01"), scatter_df$Date), format = "%Y-%m-%d")

# Filter data for quick explorations 
scatter_df_sub <- subset(scatter_df, subset = Damaged_relatedness != "Cell QC unrelated")
scatter_df_sub$Categories <- NULL
scatter_df_sub$Date <- NULL
scatter_df_sub$DOI <- NULL

# Statistics of tools only capable of damaged cell removal 
scatter_df_sub <- subset(scatter_df, subset = Damaged_relatedness == "Damage focus")
median(scatter_df_sub$Citations)
sd(scatter_df_sub$Citations)

# Statistics of popular low quality cell removal tools 
scatter_df_sub <- subset(scatter_df, subset = Name %in% c("DoubletFinder", "scrublet", "DropletUtils"))
median(scatter_df_sub$Citations)
sd(scatter_df_sub$Citations)

# Label non QC-related tools 
scatter_df$Damaged_relatedness <- ifelse(is.na(scatter_df$Damaged_relatedness), "Cell QC unrelated", scatter_df$Damaged_relatedness)

# Add a small constant to the Citations column to avoid log(0) issues
scatter_df$Citations <- scatter_df$Citations + 1

# Correct the color mapping key
scatter_df$Damaged_relatedness <- ifelse(scatter_df$Damaged_relatedness == "General cell QC ", "General cell QC", scatter_df$Damaged_relatedness)

# Define the colors for each category
color_mapping <- c("QC unrelated"          = "#EDEEF1",
                   "Cell QC unrelated"     = "#EDEEF1",
                   "Doublet focus"         = "#8DC5BD", 
                   "Damage focus"          = "#011E5C",  
                   "Empty droplet focus"   = "#CF9EBB", 
                   "General cell QC"       = "#B3B5BC")

# Filter data for plotting 
scatter_df$viewing <- ifelse(scatter_df$Citations <= 500 & scatter_df$Damaged_relatedness == "Cell QC unrelated" & scatter_df$Date >= as.Date("2018-01-01"), "Filter", "Keep")
table(scatter_df$viewing)
scatter_df_sub <- subset(scatter_df, scatter_df$viewing == "Keep")

# Create the scatter plot
scatter <- ggplot(scatter_df_sub, aes(x = Date, y = Citations, color = Damaged_relatedness)) +
  geom_point(aes(fill = Damaged_relatedness), size = 4, shape = 21, stroke = 0.5, color = "grey") +
  scale_fill_manual(values = color_mapping) +
  scale_y_log10() +
  scale_x_date(breaks = c(as.Date("2010-01-01"), as.Date("2012-01-01"), as.Date("2014-01-01"), as.Date("2016-01-01"), as.Date("2018-01-01"), as.Date("2020-01-01"), as.Date("2022-01-01"), as.Date("2024-01-01")), labels = c("2010", "2012", "2014", "2016", "2018", "2020", "2022", "2024")) +
  labs(y = "log(Citations)") +
  theme_classic() + 
  theme(
    axis.title.y = element_text(vjust = 2, face = "bold", size = 12),
    axis.text.x = element_text(hjust = 0, size = 12),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    panel.background = element_blank(),      
    panel.border = element_blank(),
    plot.background = element_blank()       
  ) +
  geom_text_repel(data = subset(scatter_df_sub, (Damaged_relatedness %in% c("Damage focus") | Name %in% c("Scater", "DoubletFinder", "Scrublet", "DropletUtils", "valiDrops"))),
                  aes(label = Name),
                  size = 5,
                  color = "black",
                  nudge_y = 0.2,
                  nudge_x = -0.3,
                  segment.color = "grey",
                  segment.size = 0.5)

# Save (with transparent background)
ggsave("./A_Search_Strategies/img/scRNAseq-tools-subset.png", 
       plot = scatter, width = 12, height = 6, units = "in", dpi = 300, bg = "transparent")


scatter_full <- ggplot(scatter_df, aes(x = Date, y = Citations, color = Damaged_relatedness)) +
  geom_point(aes(fill = Damaged_relatedness), size = 4, shape = 21, stroke = 0.5, color = "grey") +
  scale_fill_manual(values = color_mapping) +
  scale_y_log10() +
  scale_x_date(breaks = c(as.Date("2010-01-01"), as.Date("2012-01-01"), as.Date("2014-01-01"), as.Date("2016-01-01"), as.Date("2018-01-01"), as.Date("2020-01-01"), as.Date("2022-01-01"), as.Date("2024-01-01")), labels = c("2010", "2012", "2014", "2016", "2018", "2020", "2022", "2024")) +
  labs(y = "log(Citations)") +
  theme_classic() + 
  theme(
    axis.title.y = element_text(vjust = 2, face = "bold", size = 12),
    axis.text.x = element_text(hjust = 0, size = 12),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    panel.background = element_blank(),      
    panel.border = element_blank(),
    plot.background = element_blank()       
  ) +
  geom_text_repel(data = subset(scatter_df_sub, (Damaged_relatedness %in% c("Damage focus") | Name %in% c("Scater", "DoubletFinder", "Scrublet", "DropletUtils", "valiDrops"))),
                  aes(label = Name),
                  size = 5,
                  color = "black",
                  nudge_y = 0.2,
                  nudge_x = -0.3,
                  segment.color = "grey",
                  segment.size = 0.5)

# Save (with transparent background)
ggsave("./A_Search_Strategies/img/scRNAseq-tools.png", 
       plot = scatter_full, width = 12, height = 6, units = "in", dpi = 300, bg = "transparent")



scatter_full_no_labels <- ggplot(scatter_df, aes(x = Date, y = Citations, color = Damaged_relatedness)) +
  geom_point(aes(fill = Damaged_relatedness), size = 4, shape = 21, stroke = 0.5, color = "grey") +
  scale_fill_manual(values = color_mapping) +
  scale_y_log10() +
  scale_x_date(breaks = c(as.Date("2010-01-01"), as.Date("2012-01-01"), as.Date("2014-01-01"), as.Date("2016-01-01"), as.Date("2018-01-01"), as.Date("2020-01-01"), as.Date("2022-01-01"), as.Date("2024-01-01")), labels = c("2010", "2012", "2014", "2016", "2018", "2020", "2022", "2024")) +
  labs(y = "log(Citations)") +
  theme_classic() + 
  theme(
    axis.title.y = element_text(vjust = 2, face = "bold", size = 12),
    axis.text.x = element_text(hjust = 0, size = 12),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    panel.background = element_blank(),      
    panel.border = element_blank(),
    plot.background = element_blank()       
  ) 

# Save (with transparent background)
ggsave("./A_Search_Strategies/img/scRNAseq-tools-no-labels.png", 
       plot = scatter_full_no_labels, width = 12, height = 6, units = "in", dpi = 300, bg = "transparent")



### End 
