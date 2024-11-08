# CONTEXT 
#
# To get a snapshot of the damaged protocols currently used in the field, papers from the 
# past decade that have analysed single cell data will be collected and protocols therein will be reviewd.
# This script allows for that collection from PubMed. 


#-------------------------------------------------------------------------------
# Libraries 
#-------------------------------------------------------------------------------

install.packages("rentrez")
library(rentrez)

# Set the working directory to the Zenodo directory 
# setwd("/Users/name/Zenodo")

#-------------------------------------------------------------------------------
# Perform the search 
#-------------------------------------------------------------------------------

# Define years of interest 
years <- c(2024, 2023, 2022, 2020, 2019, 2018, 2017, 2016, 2015)
sample_size <- 30
sampled_papers <- list()

for (year in years) {
  # Define search query with date range for each year
  query <- sprintf('"single cell RNA sequencing" OR "scRNA-seq" OR "single-cell transcriptomics" AND (dataset OR data OR protocol OR experiment) AND ("%d"[PDAT] : "%d"[PDAT])', year, year)
  
  # Perform search on PubMed
  search_results <- entrez_search(db="pubmed", term=query, retmax=150, use_history=TRUE)
  
  # Use search IDs to retrieve paper information 
  search_ids <- search_results$ids
  papers <- entrez_summary(db="pubmed", id=search_ids)
  papers_df <- data.frame(title = sapply(papers, function(x) x$title),
                          pubdate = sapply(papers, function(x) x$pubdate))
  
  # Filter papers without complete dates (day, month, and year)
  papers_df <- papers_df[grepl("^\\d{4} [A-Za-z]{3} \\d{2}$", papers_df$pubdate), ]
  
  # Extract year and month
  papers_df$year <- as.numeric(sub("^(\\d{4}).*", "\\1", papers_df$pubdate))
  papers_df$month <- sub("^\\d{4} ([A-Za-z]{3}).*", "\\1", papers_df$pubdate)
  
  year_sample <- data.frame()
  
  for (month in month.abb) {
    monthly_papers <- papers_df[papers_df$month == month, ]
    
    # If there are more papers than needed, sample approx 2.5 (30/12) from each month
    if (nrow(monthly_papers) > 0) {
      month_sample <- monthly_papers[sample(seq_len(nrow(monthly_papers)), min(3, nrow(monthly_papers))), ]
      year_sample <- rbind(year_sample, month_sample)
    }
  }
  
  # If there are fewer than sample_size, sample additional papers
  if (nrow(year_sample) < sample_size) {
    remaining_needed <- sample_size - nrow(year_sample)
    additional_papers <- papers_df[!papers_df$title %in% year_sample$title, ]
    
    if (nrow(additional_papers) > 0) {
      additional_sample <- additional_papers[sample(seq_len(nrow(additional_papers)), min(remaining_needed, nrow(additional_papers))), ]
      year_sample <- rbind(year_sample, additional_sample)
    }
  }
  
  sampled_papers[[as.character(year)]] <- year_sample
}

# Combine all samples 
final_sampled_papers <- do.call(rbind, sampled_papers)

# Save the search result to a file for reproducibility
write.csv(final_sampled_papers, "./A_Search_Strategies/data/PubMed_search_papers.csv", quote = FALSE, row.names = FALSE)


### End 
