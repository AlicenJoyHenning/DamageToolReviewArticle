#


# Install and load rentrez package
install.packages("rentrez")
library(rentrez)

#-------------------------------------------------------------------------------
# Perform the search
#-------------------------------------------------------------------------------

# Define search query
query <- '"single cell RNA sequencing" AND protocol AND dataset'

# Perform search on PubMed
search_results <- entrez_search(db="pubmed", term=query, retmax=1000, use_history=TRUE)

# Get citation counts for each paper (you may need to combine this with a citation database like CrossRef)
search_ids <- search_results$ids

# Extract details of top-cited papers
papers <- entrez_summary(db="pubmed", id=search_ids)

# Convert to a data frame
papers_df <- data.frame(
  title = sapply(papers, function(x) x$title),
  citation_count = sapply(papers, function(x) x$citation_count), # Use another source for citation counts if needed
  pubdate = sapply(papers, function(x) x$pubdate)
)

# Sort by citation count (if citation data is available)
top_100_papers <- papers_df[order(-papers_df$citation_count), ][1:100,]

# Save the search result to a file for reproducibility
write.csv(top_100_papers, "scRNA_seq_top_100_papers.csv")
