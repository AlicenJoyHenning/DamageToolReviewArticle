# CONTEXT 
# 
# Using scrna-tools table, https://www.scrna-tools.org/tools, explore the tools 
# available that are related to cell quality control, specifically damaged cell detection. 

#-------------------------------------------------------------------------------
# Libraries 
#-------------------------------------------------------------------------------

# Load packages -----

install.packages("rcrossref")
library(rcrossref)          
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(cowplot)

#-------------------------------------------------------------------------------
# Data cleaning and manual editing 
#-------------------------------------------------------------------------------

# Prepare for visualisations -----

# Read in the data frame downloaded from https://www.scrna-tools.org/table
scRNAseqtools_df <- read.csv("/home/alicen/Projects/ReviewArticle/article_search/scRNA-tools-dates.csv")

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
write.csv2(scRNAseqtools_dff, "Projects/ReviewArticle/article_search/scRNAseq-tools-unknown-to-be-filled-in.csv")

# Continue after manual entry
scRNAseqtools_dff <- read.csv("Projects/ReviewArticle/article_search/scRNAseq-tools-unknown-to-be-filled-in.csv")
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

# Manually search to see which of the QC-focus tools can detect damaged cells 
test <- subset(scRNAseqtools_dff, scRNAseqtools_dff$Colour == "QC-focus")
test$Categories <- NULL


write.csv(test,
          quote = FALSE,
          row.names = FALSE,
          "C:/Users/alice/OneDrive/Documents/University/Masters/limiric_article/complete_details.csv"
)

damaged_focus <- read.csv("Projects/ReviewArticle/article_search/complete_details.csv", sep = ",", header = TRUE)
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

damaged_focus  <- read.csv("/home/alicen/Projects/ReviewArticle/article_search/scRNA-tools-damage-focus.csv", sep = ",", header = FALSE)
colnames(damaged_focus)[colnames(damaged_focus) == "V1"] <- "Name"
colnames(damaged_focus)[colnames(damaged_focus) == "V2"] <- "Damaged"


#-------------------------------------------------------------------------------
# Plot
#-------------------------------------------------------------------------------

damaged_focus  <- read.csv("/home/alicen/Projects/ReviewArticle/article_search/scRNA-tools-damage-focus.csv", sep = ",", header = FALSE)

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


# Ensure the Date column is of Date type
scatter_df <- merge(scRNAseqtools_dff, damaged_focus, by = "Name", all.x = TRUE)

# Rename & clean our columns
scatter_df$Colour_Final <- ifelse(is.na(scatter_df$Damaged), scatter_df$Colour, scatter_df$Damaged)
scatter_df$Colour_Final <- ifelse(scatter_df$Colour_Final %in% c("QC-unrelated", "QC-related"), "QC unrelated", scatter_df$Colour_Final)
scatter_df$Date <- as.Date(scatter_df$Date)
scatter_df$Colour_Final <- ifelse(scatter_df$Name  %in% c("ddqc", "DropletQC", "EnsembleKQC", "miQC", "scater", "scuttle", "valiDrops"), "Damage focus", scatter_df$Colour_Final)

# Subset 
scatter_df <- subset(scatter_df, !is.null(scatter_df$Colour_Final))
scatter_df <- subset(scatter_df, !is.na(scatter_df$Colour_Final))

# Add a small constant to the Citations column to avoid log(0) issues
scatter_df$Citations <- scatter_df$Citations + 1

# Define the colors for each category
unique(scatter_df$Colour_Final)
color_mapping <- c("QC unrelated"          = "#EDEEF1",
                   "Cell QC unrelated"     = "#EDEEF1",
                   "Doublet focus"         = "#8DC5BD", 
                   "Damage focus"          = "#011E5C",  
                   "Empty droplet focus"   = "#CF9EBB", 
                   "General cell QC "      = "#B3B5BC")

scatter_df$viewing <- ifelse(scatter_df$Citations <= 500 & scatter_df$Colour_Final == "QC unrelated" & scatter_df$Date >= year_2016, "Filter", "Keep")
scatter_df_sub <- subset(scatter_df, scatter_df$viewing == "Keep")

# Create the scatter plot
scatter <- ggplot(scatter_df_sub, aes(x = Date, y = Citations, color = Colour_Final)) +
  geom_point(aes(fill = Colour_Final), size = 4, shape = 21, stroke = 0.5, color = "grey") +
  scale_fill_manual(values = color_mapping) +
  scale_y_log10() +
  scale_x_date(breaks = c(year_2010, year_2012, year_2014, year_2016, year_2018, year_2020, year_2022, year_2024), labels = c("2010", "2012", "2014", "2016", "2018", "2020", "2022", "2024")) +
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
  +
  geom_text_repel(data = subset(scatter_df, Colour_Final %in% c("Damage focus")),
  aes(label = Name),
  size = 5,
  color = "black",
  nudge_y = 0.2,
  nudge_x = -0.3,
  segment.color = "grey",
  segment.size = 0.5)

# Save (with transparent background)
ggsave("Projects/ReviewArticle/article_search/scRNAseq-tools-no-labels.png", 
       plot = scatter, width = 12, height = 6, units = "in", dpi = 300, bg = "transparent")



### End 