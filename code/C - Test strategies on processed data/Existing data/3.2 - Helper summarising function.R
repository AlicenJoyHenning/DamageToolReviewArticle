# SCRIPT CONTEXT 
#
# Helper function for the benchmarking pipeline to summarise the output of benchmark function (3 - Run remaining detection strategies)
# Calculates the percentage of damaged cells detected for each tool. 


#-------------------------------------------------------------------------------
# HELPER FUNCTION DEFINED 
#-------------------------------------------------------------------------------

# Helper function for benchmarking 

summarise_results <- function(seurat, 
                              methods = list("ddqc", "DropletQC", "ensembleKQC", "miQC", "scater", "valiDrops", 
                                             "manual_fixed_mito", "manual_adaptive_mito", "manual_mito_ribo",  "manual_mito_ribo_library",  "manual_malat1", "manual_malat1_mito_ribo")
                              ){
  
  # Function to automate live and dead cell calculations for each method
  calculate_percentage <- function(data, method) {
    
    data %>%
      mutate(Category = data[[method]]) %>%
      group_by(Category) %>%
      summarise(Count = n(), .groups = 'drop') %>%
      mutate(Percentage = (Count / sum(Count)) * 100)
    }
  
  # Calculate percentages for each method and store in a list
  percentages_list <- lapply(methods, function(method) {
    calculate_percentage(seurat@meta.data, method) %>%
      mutate(Method = method)
  })

  # Combine all percentages into a single data frame
  combined_percentages <- bind_rows(percentages_list)

  # Reshape the data frame
  final_df <- combined_percentages %>%
    pivot_wider(names_from = Category, values_from = Percentage, values_fill = list(Percentage = 0))

  final_df$Count <- NULL

  # Ensure the columns are named correctly
  final_df <- final_df %>%
    rename_with(~ c("Method", "cell", "damaged"), everything())

  # Collapse repeated rows by retaining the maximum value for each column
  final_df <- final_df %>%
    group_by(Method) %>%
    summarise(across(everything(), \(x) max(x, na.rm = TRUE)))

  return(final_df)
  
}

### End
