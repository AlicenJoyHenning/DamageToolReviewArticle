# Source script for valiDrops label_dead function
# Usually comes as last step in process of uploading raw data, identifying empty droplets etc. but since 
# I'm only interested in comparing dead/damaged cells, I'm only going to benchmark this function 

label_dead <- function(
    counts, 
    metrics) {
  
  # Default values
  cor.threshold <- NULL
  train <- TRUE
  rep <- 10
  n.min <- 8
  n.relabel <- 1
  feature.try <- 3
  verbose <- FALSE
  label.thrs <- NULL
  label.frac <- 0.1
  nfeats <- 2000
  alpha <- 0
  npcs <- 100
  weight <- TRUE
  epochs <- 20
  nfolds <- 5
  nrep <- 10
  fail.weight <- 0.2
  cor.min <- 0.0001
  cor.max <- 0.005
  cor.steps <- 50
  nrep.cor <- 10
  min.dead <- 100
  max.live <- 500
  plot <- TRUE
  
  
  
  # Soft label the dataset
  metrics$logUMIs <- scale(metrics$logUMIs, scale = FALSE)
  metrics$logFeatures <- scale(metrics$logFeatures, scale = FALSE)
  metrics$ribosomal_fraction <- asin(sqrt((metrics$ribosomal_fraction)))/(pi/2)
  metrics$coding_fraction <- asin(sqrt((metrics$coding_fraction)))/(pi/2)
  metrics$mitochondrial_fraction <- asin(sqrt(metrics$mitochondrial_fraction))/(pi/2)
  metrics$score <- metrics$logUMIs * -11.82 + metrics$logFeatures * 2.08 + metrics$ribosomal_fraction * 158.98 + metrics$logFeatures * metrics$coding_fraction * 18.87 + metrics$ribosomal_fraction * metrics$coding_fraction * -125.9
  
  # Define a flag for a successful run
  flag <- "Success" 
  
  # Determine the optimal threshold
  if (is.null(label.thrs)) {
    if (verbose) { message("Optimizing labeling threshold") }
    # Try default quantile cutoff
    max.quantile <- 0.1
    break.data <- as.data.frame(matrix(ncol=2, nrow = 200))
    cnt <- 1
    for (brk in seq(0.0001, max.quantile, 0.0001)) {
      break.data[cnt,1] <- brk
      break.data[cnt,2] <- quantile(metrics$score, brk)
      cnt <- cnt + 1
    }
    last.thrs <- quantile(metrics$score, inflection::uik(break.data$V1, break.data$V2))
    metrics$label <- "live"
    metrics[ metrics$score <= last.thrs,"label"] <- "dead"
    last.distribution <- table(metrics$label)
    
    # Increase quantile cutoff and evaluate
    stop <- 0
    while (stop == 0) {
      # Label using a new cutoff
      max.quantile <- max.quantile + 0.1
      break.data <- as.data.frame(matrix(ncol=2, nrow = 200))
      cnt <- 1
      for (brk in seq(0.0001, max.quantile, 0.0001)) {
        break.data[cnt,1] <- brk
        break.data[cnt,2] <- quantile(metrics$score, brk)
        cnt <- cnt + 1
      }
      new.thrs <- quantile(metrics$score, inflection::uik(break.data$V1, break.data$V2))
      metrics$label <- "live"
      metrics[ metrics$score <= new.thrs,"label"] <- "dead"
      new.distribution <- table(metrics$label)
      
      # Evaluate whether or not to stop
      if (min(new.distribution) > 0) {
        if (min(last.distribution) > 0) {
          if (last.distribution[1] == new.distribution[1]) {
            last.distribution <- new.distribution
            last.thrs <- new.thrs
          } else if (last.distribution[1] < new.distribution[1]) {
            label.thrs <- last.thrs
            stop <- 1
          }
        } else {
          label.thrs <- new.thrs
          stop <- 1
        }
      } else {
        if (max.quantile >= 0.95) {
          label.thrs <- new.thrs
          flag <- "Failed"
          stop <- 1
        } else {
          last.distribution <- new.distribution
          last.thrs <- new.thrs
        }
      }
    }
    
    # Label using final cutoff
    metrics$label <- "live"
    metrics[ metrics$score <= label.thrs,"label"] <- "dead"
  } else {
    # Label the barcodes
    metrics$label <- "live"
    metrics[ metrics$score <= label.thrs,"label"] <- "dead"
  }
  
  # Break if there are not enough observations or if there are too many
  if (nrow(metrics[ metrics$label == "dead",]) < 3) {
    if (verbose) { message("Soft-labeling identified less than 3 barcodes associated with dead cells. Returning results from soft-thresholding") }
    return(list(labels = metrics$label, metrics = metrics, flag = "Failed"))
  }
  
  # Uncertain fraction calculation
  uncertain.fraction <- nrow(metrics[ metrics$label == "uncertain",]) / nrow(metrics)
  if (uncertain.fraction >= 0.0125 & uncertain.fraction < 0.025) {
    if (verbose) { message("There is a medium number of barcodes with uncertain classification. Interpret results with caution!") }
    flag <- "Caution"
  } else if (uncertain.fraction >= 0.025) {
    if (verbose) { message("There is a high number of barcodes with uncertain classification. The model did not converge. Interpret with extreme caution.") }
    flag <- "Failed"
  }
  
  # Return
  return(list(labels = metrics$label, metrics = metrics, flag = flag))
}
