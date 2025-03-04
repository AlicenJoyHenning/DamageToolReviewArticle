


# Load the matrix
PH_CT6_matrix <- readRDS("/Users/alicen/Projects/Damage_analsyis/damage_perturbation/scDesign2/modelling/unfiltered_genes/models/PBMC_high_CT6_reference_matrix.rds")


PH_CT6_matrix_fix <- PH_CT6_matrix
colnames(PH_CT6_matrix_fix) <- names(colnames(PH_CT6_matrix))

mito_idx <- grep("^MT-", rownames(PH_CT6_matrix_fix), ignore.case = FALSE)
mito <- colSums(PH_CT6_matrix_fix[mito_idx, ]) / colSums(PH_CT6_matrix_fix)
filtered_cells <- names(mito[mito < 0.1])
PH_CT6_filtered_matrix <- PH_CT6_matrix_fix[, filtered_cells]
mito <- colSums(PH_CT6_filtered_matrix[mito_idx, ]) / colSums(PH_CT6_filtered_matrix)

df <- data.frame(cells = colnames(PH_CT6_filtered_matrix))
df$features <- colSums(PH_CT6_filtered_matrix != 0)
df$mito <- colSums(PH_CT6_filtered_matrix[mito_idx, ]) / colSums(PH_CT6_filtered_matrix)

SimulateDamage <- function(
    counts,            # Single cell count matrix
    damage_proportion, # What is the target proportion of cells to be damaged? 
    annotated = FALSE, # Whether or not the cells have labels associated with cell type
    target_damage = c(0.1, 0.8), # Range of damage
    damage_steepness = "moderate",
    damage_distribution = "right_skewed",
    penalty_factor = 0.01,  # Reduce the probability for ribosomal genes by 50%
    organism = "Hsap"
) {
  # Calculate the number of damaged cells given target damage proportion 
  total_cells <- ncol(counts)
  damaged_cell_number <- round(total_cells * damage_proportion)
  
  if (annotated){
    
    # Damage cell selections must be distributed across cell types evenly
    cell_types <- as.factor(sub("_.*", "", colnames(counts))) 
    cell_type_counts <- table(cell_types)
    damage_per_type <- round(cell_type_counts * (damage_proportion))
    
    # Adjust to ensure the total matches the target
    total_damaged <- sum(damage_per_type)
    
    # Adjust if necessary
    if (total_damaged != damaged_cell_number) {
      diff <- damaged_cell_number - total_damaged
      adj_indices <- sample(seq_along(damage_per_type), abs(diff), replace = TRUE)
      damage_per_type[adj_indices] <- damage_per_type[adj_indices] + sign(diff)
    }
    
    # Select damaged cells
    damaged_cell_selections <- unlist(lapply(names(damage_per_type), function(ct) {
      cells_of_type <- which(cell_types == ct)
      sample(cells_of_type, size = damage_per_type[ct], replace = FALSE)
    }))
    
    # Make the column names unique (no longer need the celltypes in isolation)
    colnames(counts) <- paste0(colnames(counts), "_", seq_len(dim(counts)[2])) # Keep cell types but ensure uniqueness
    
  } else {
    
    # If no cell types specified, cells are sampled randomly 
    damaged_cell_selections <- sample(seq_len(total_cells), damaged_cell_number)
    
  }
  
  # Storage of damage levels for all cell barcodes for plotting later
  damage_label <- data.frame(barcode = colnames(counts)[damaged_cell_selections], status = rep("damaged", length(damaged_cell_selections)))
  undamaged_cell_number_cells <- setdiff(seq_len(total_cells), damaged_cell_selections)
  undamaged_cell_number <- data.frame(barcode = colnames(counts)[undamaged_cell_number_cells], status = rep("control", length(undamaged_cell_number_cells)))
  damage_label <- rbind(damage_label, undamaged_cell_number)
  
  # Assign damage levels to the selected cells based on beta distribution 
  
  # Assign the steepness (peak height) to meet target
  steepness_levels <- list(
    shallow = 4,
    moderate = 7,
    steep = 12
  )
  
  # Retrieve user specified steepness, default is "moderate"
  steepness_value <- steepness_levels[[damage_steepness]]
  
  # Define shape parameters to achieve target distribution, default "right_skewed"
  if (damage_distribution == "right_skewed") {
    a <- steepness_value * 0.3
    b <- steepness_value * 0.7
  } else if (damage_distribution == "left_skewed") {
    a <- steepness_value * 0.7
    b <- steepness_value * 0.3
  } else if (damage_distribution == "symmetric") {
    a <- steepness_value * 0.5
    b <- steepness_value * 0.5
  }
  
  # Consolidate into target beta distribution 
  damage_levels <- rbeta(length(damaged_cell_selections), shape1 = a, shape2 = b)
  
  # Scale values lie within target range, default c(0.1, 0.8)
  damage_levels <- target_damage[1] + (target_damage[2] - target_damage[1]) * damage_levels 
  
  # Store the assigned damage levels in a dataframe
  damage_label$damage_level <- 0  # Default is 0 for undamaged cells
  damage_label$damage_level[match(colnames(counts)[damaged_cell_selections], damage_label$barcode)] <- damage_levels
  
  # Retrieve genes corresponding to the organism of interest
  if (organism == "Hsap") {
    mito_pattern <- "^MT-"
    ribo_pattern <- "^(RPS|RPL)"
    nuclear <- c("FIRRE", "NEAT1","XIST", "MALAT1", "MIAT", "MEG3", "KCNQ1OT1", "HOXA11-AS", "FTX") 
  }
  
  if (organism == "Mmus") {
    mito_pattern <- "^mt-"
    ribo_pattern <- "^(rps|rpl)"
    nuclear <- c("Firre", "Neat1","Xist", "Malat1", "Miat", "Meg3", "Kcnq1ot1", "Hoxa11-as", "Ftx") 
  }
  
  # Isolate gene set indices (consistent across cells, not subsetting the matrix)
  mito_idx <- grep(mito_pattern, rownames(counts), ignore.case = FALSE)
  nucl_idx <- which(rownames(counts) %in% nuclear)
  mito_idx <- c(mito_idx, nucl_idx)
  non_mito_idx <- setdiff(seq_len(nrow(counts)), mito_idx)
  ribo_idx <- grep(ribo_pattern, rownames(counts), ignore.case = FALSE)
  
  # Initialize for storing modified counts
  damaged_counts <- counts
  
  for (cell in damaged_cell_selections) {
    # Return unchanged if no loss is applied
    if (damage_label$damage_level[match(colnames(counts)[cell], damage_label$barcode)] == 0) {
      next
    } else {
      
      # Find total non-mito transcripts
      total_counts <- sum(damaged_counts[non_mito_idx, cell])
      
      # Determine number of transcripts to lose
      total_loss <- round(damage_label$damage_level[match(colnames(counts)[cell], damage_label$barcode)] * total_counts)
      
      # Expand genes into individual transcript-level representations
      transcript_df <- rep(non_mito_idx, times = counts[non_mito_idx, cell])
      
      # Assign probability of loss based on gene abundance
      gene_totals <- damaged_counts[non_mito_idx, cell]
      probabilities <- gene_totals / total_counts
      
      # Apply penalty to ribosomal genes
      probabilities[ribo_idx] <- probabilities[ribo_idx] * penalty_factor
      
      # Normalize probabilities so that the total sums to 1
      probabilities <- probabilities / sum(probabilities)
      
      # Sample transcripts to be lost
      lost_transcripts <- sample(
        transcript_df,  
        size = total_loss, 
        replace = FALSE, 
        prob = rep(probabilities, times = gene_totals)
      )
      
      # Sum remaining transcripts per gene
      remaining_counts <- table(factor(transcript_df[!transcript_df %in% lost_transcripts], levels = non_mito_idx))
      
      # Assign updated counts to non-mito genes
      damaged_counts[non_mito_idx, cell] <- as.integer(remaining_counts)
    }
  }
  
  # isolated indices for correct ordering of damage status 
  matched_indices <- match(colnames(damaged_counts), damage_label$barcode)
  
  # Generate qc_summary with numeric values
  qc_summary <- data.frame(
    Cell = colnames(counts),
    Damaged_Level = as.numeric(damage_label$damage_level[matched_indices]),
    Original_Features = as.numeric(colSums(counts != 0)),
    New_Features = as.numeric(colSums(damaged_counts != 0)),
    Original_MitoProp = as.numeric(colSums(counts[mito_idx, , drop = FALSE]) / colSums(counts)),
    New_MitoProp = as.numeric(colSums(damaged_counts[mito_idx, , drop = FALSE]) / colSums(damaged_counts)),
    Original_RiboProp = as.numeric(colSums(counts[ribo_idx, , drop = FALSE]) / colSums(counts)),
    New_RiboProp = as.numeric(colSums(damaged_counts[ribo_idx, , drop = FALSE]) / colSums(damaged_counts))
  )
  
  # Function to generate QC plots
  generate_plot <- function(data, x, y, altered = FALSE, legend = FALSE, damage_column = "Damaged_Level") {
    
    if (altered) {
      p <- ggplot(data, aes_string(x = x, y = y, colour = damage_column)) + 
        scale_color_gradientn(
          colours = c("grey", "#7023FD", "#E60006"), 
          values = scales::rescale(c(0, 0.3, 1)), 
          limits = c(0, 1)
        ) + 
        scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + 
        geom_point(size = 0.2) + 
        theme_classic() + 
        theme(
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(fill = NA, color = "black"),
          legend.position = "bottom",
          legend.justification = "center",
          legend.title = element_text(face = "bold", hjust = 0.5) # Center and bold the legend title
        ) + 
        guides(color = guide_colorbar(title = "Damage Level", title.position = "top"))
      
      return(p)
      
    } else {
      p <- ggplot(data, aes_string(x = x, y = y)) + 
        scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + 
        geom_point(size = 0.2, color = "gray") + 
        theme_classic() + 
        theme(
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(fill = NA, color = "black"),
          legend.position = "none"
        )
      
      return(p)
    }
  }
  
  # Generate individual plots 
  mito_ribo_old <- generate_plot(qc_summary, x = "Original_RiboProp", y = "Original_MitoProp") + 
    scale_x_continuous(limits = c(0, 1), labels = scales::number_format(accuracy = 0.1)) + 
    labs(x = "Ribosomal proportion", y = "Mitochondrial proportion")
  
  mito_ribo_new <- generate_plot(qc_summary, x = "New_RiboProp", y = "New_MitoProp", altered = TRUE) + 
    scale_x_continuous(limits = c(0, 1), labels = scales::number_format(accuracy = 0.1))  + 
    labs(x = "Ribosomal proportion", y = "Mitochondrial proportion")
  
  mito_features_old <- generate_plot(qc_summary, x = "Original_Features", y = "Original_MitoProp") + 
    labs(x = "Features expressed", y = "Mitochondrial proportion")
  
  mito_features_new <- generate_plot(qc_summary, x = "New_Features", y = "New_MitoProp", altered = TRUE) + 
    labs(x = "Features expressed", y = "Mitochondrial proportion")
  
  # Extract the legend from mito_ribo_new (ensure one legend is present)
  legend <- ggpubr::get_legend(mito_ribo_new)
  
  # Create titles for the plots
  title_original <- ggdraw() + 
    draw_label("Original counts", fontface = 'bold', hjust = 0.5)
  
  title_altered <- ggdraw() + 
    draw_label("Altered counts", fontface = 'bold',  hjust = 0.5)
  
  # Arrange original plots in a single row
  original_plots <- plot_grid(mito_features_old, mito_ribo_old, ncol = 2)
  
  # Arrange altered plots in a single row 
  mito_ribo_new_no_legend <- mito_ribo_new + theme(legend.position = "none")
  mito_features_new_no_legend <- mito_features_new + theme(legend.position = "none")
  altered_plots <- plot_grid(mito_features_new_no_legend, mito_ribo_new_no_legend, ncol = 2)
  
  # Combine the original and altered rows with their titles
  original_with_title <- plot_grid(title_original, original_plots, ncol = 1, rel_heights = c(0.2, 1))
  altered_with_title <- plot_grid(title_altered, altered_plots, ncol = 1, rel_heights = c(0.2, 1))
  
  # Combine the original and altered rows, and position the legend in its own row
  final_plot <- plot_grid(
    original_with_title, 
    altered_with_title, 
    legend,  # Legend in its own row
    ncol = 1, 
    rel_heights = c(1, 1, 0.25)  # Adjust relative height to fit legend
  )
  
  # Increase margins around the total plot area
  final_plot <- final_plot + 
    theme(
      plot.margin = margin(10, 20, 20, 20)  # Increase margins (top, right, bottom, left)
    )
  
  # Display the final plot
  print(final_plot)
  
  
  # Return a list containing items of interest
  return(list(matrix = damaged_counts, qc_summary = qc_summary, plot = final_plot))
}

output <- SimulateDamage(counts = PH_CT6_filtered_matrix,
                         damage_proportion = 0.5,
                         annotated = FALSE)  


