# SCRIPT CONTEXT 
#
# Gather the following objective usability criterion, 
#    1 - Documentation present (Yes - 1, No - 0) 
#    2 - More than one GitHub/ deployment Upload (Yes - 1, No - 0)
#    3 - Last adjustment to tool less than two years ago (Yes - 1, No - 0)
#    4 - GitHub stars rescaled ( * / 5 ) 
#
# Using GitHub actions, the criterion will be automatically updated everyday. 

# Install required libraries
if (!requireNamespace("httr")) install.packages("httr")
if (!requireNamespace("tidyverse")) install.packages("tidyverse")
library(httr)
library(tidyverse)

# Define repository list
repositories <- data.frame(
  method = c("ddqc", "DropletQC", "ensembleKQC", "miQC", "SampleQC", 
             "scater", "valiDrops", "manual_fixed_mito", "manual_adaptive_mito"),
  repo = c("ayshwaryas/ddqc", 
           "powellgenomicslab/DropletQC", 
           "mzhq/ensembleKQC", 
           "greenelab/miQC", 
           "wmacnair/SampleQC", 
           "alanocallaghan/scater", 
           "madsen-lab/valiDrops", 
           NA, NA) 
)


# Function to fetch stars
get_stars <- function(repo) {
  if (is.na(repo)) {
    return(0)  # Return 0 for non-GitHub repos
  }
  url <- paste0("https://api.github.com/repos/", repo)
  response <- GET(url)
  if (status_code(response) == 200) {
    content(response)$stargazers_count
  } else {
    warning(paste("Failed to fetch:", repo))
    return(0)
  }
}

# Fetch stars and save as CSV
repositories <- repositories %>%
  mutate(Stars = map_int(repo, get_stars)) %>%
  select(method, Stars)

# Write to CSV
write.csv(repositories, "data/D_Summary_Results/usability-flexible.csv", row.names = FALSE)
