# CONTEXT 
#
# To get a snapshot of the damaged protocols currently used in the field, papers from the 
# past decade that have analysed single cell data will be collected and protocols therein will be reviewd.
# This script allows for that collection from PubMed. 


#-------------------------------------------------------------------------------
# Libraries 
#-------------------------------------------------------------------------------

# Install and load rentrez package
install.packages("rentrez")
library(rentrez)

#-------------------------------------------------------------------------------
# Perform the search
#-------------------------------------------------------------------------------

# Define search query (PubMed)
query <- '"single cell RNA sequencing" OR "scRNA-seq" OR "single-cell transcriptomics" AND (dataset OR data OR protocol OR experiment)'

# Perform search on PubMed
search_results <- entrez_search(db="pubmed", term=query, retmax=150, use_history=TRUE)

# Use search IDs to retrieve paper information 
search_ids <- search_results$ids
papers     <- entrez_summary(db="pubmed", id=search_ids)
papers_df  <- data.frame(title = sapply(papers_list, function(x) x$title),
                         pubdate = sapply(papers_list, function(x) x$pubdate)
                         )

# Filter papers without complete dates (day, month, and year)
papers_df <- papers_df[grepl("^\\d{4} [A-Za-z]{3} \\d{2}$", papers_df$pubdate), ] # 9999 > 3864


# Select subset of papers from each year for review

# Add columns for filtering
papers_df$year <- as.numeric(sub("^(\\d{4}).*", "\\1", papers_df$pubdate))
papers_df$month <- sub("^\\d{4} ([A-Za-z]{3}).*", "\\1", papers_df$pubdate)

# Define years of interest 
years <- c(2024, 2023, 2022, 2020, 2019)
sample_size <- 20


# Sample
sampled_papers <- list()

for (year in years) {
  
  yearly_papers <- papers_df[papers_df$year == year, ]
  
  year_sample <- data.frame()
  
  for (month in month.abb) {
    monthly_papers <- yearly_papers[yearly_papers$month == month, ]
    
    # If there are more papers than needed, sample approx 2 (20/12) from each month
    if (nrow(monthly_papers) > 0) {
      month_sample <- monthly_papers[sample(seq_len(nrow(monthly_papers)), min(2, nrow(monthly_papers))), ]
      year_sample <- rbind(year_sample, month_sample)
    }
  }
  
  # If there are fewer than sample_size, sample additional papers
  if (nrow(year_sample) < sample_size) {
    remaining_needed <- sample_size - nrow(year_sample)
    additional_papers <- yearly_papers[!yearly_papers$ids %in% year_sample$ids, ]
    additional_sample <- additional_papers[sample(seq_len(nrow(additional_papers)), remaining_needed), ]
    year_sample <- rbind(year_sample, additional_sample)
  }
  
  sampled_papers[[as.character(year)]] <- year_sample
  
}

# Combine all samples 
final_sampled_papers <- do.call(rbind, sampled_papers) # 116 papers 

# Save the search result to a file for reproducibility
write.csv(final_sampled_papers, "/home/alicen/Projects/ReviewArticle/article_search/PubMed_search_papers.csv", quote = FALSE, row.names = FALSE)


### End 