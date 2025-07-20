#!/usr/bin/env Rscript

# Script to evaluate cell counts from epiAneufinder CNV calls
# Reads cnv_calls.rds from each epianeufinder_* folder and counts cells

require(data.table)
options(stringsAsFactors = FALSE)

cat("=== epiAneufinder Cell Count Evaluation ===\n")
cat("Scanning for epianeufinder sample directories...\n")

# Find all directories starting with "epianeufinder_"
sample_dirs <- list.dirs(".", recursive = FALSE, full.names = TRUE)
sample_dirs <- sample_dirs[grepl("epianeufinder_", basename(sample_dirs))]

if(length(sample_dirs) == 0) {
  stop("No sample directories found with pattern 'epianeufinder_*'")
}

cat("Found", length(sample_dirs), "sample directories\n\n")

# Initialize results list
results <- list()

# Process each sample directory
for(sample_dir in sample_dirs) {
  sample_name <- basename(sample_dir)
  cat("Processing:", sample_name, "\n")
  
  # Path to CNV calls file
  cnv_calls_file <- file.path(sample_dir, "epiAneufinder_results", "cnv_calls.rds")
  
  if(file.exists(cnv_calls_file)) {
    tryCatch({
      # Read CNV calls
      cnv_calls <- readRDS(cnv_calls_file)
      
      # Count cells (length of the list)
      n_cells <- length(cnv_calls)
      
      cat("  Found", n_cells, "cells\n")
      
      # Store result
      results[[sample_name]] <- data.frame(
        sample = sample_name,
        n_cells = n_cells,
        file_path = cnv_calls_file,
        stringsAsFactors = FALSE
      )
      
    }, error = function(e) {
      cat("  Error reading CNV calls:", e$message, "\n")
      results[[sample_name]] <- data.frame(
        sample = sample_name,
        n_cells = NA,
        file_path = cnv_calls_file,
        error = e$message,
        stringsAsFactors = FALSE
      )
    })
  } else {
    cat("  CNV calls file not found:", cnv_calls_file, "\n")
    results[[sample_name]] <- data.frame(
      sample = sample_name,
      n_cells = NA,
      file_path = cnv_calls_file,
      error = "File not found",
      stringsAsFactors = FALSE
    )
  }
}

# Combine results
if(length(results) > 0) {
  results_df <- do.call(rbind, results)
  rownames(results_df) <- NULL
  
  # Write to output file
  output_file <- "evals_cb.txt"
  write.table(results_df, file = output_file, sep = "\t", quote = FALSE, 
              row.names = FALSE, col.names = TRUE)
  
  cat("\n=== Summary ===\n")
  cat("Results saved to:", output_file, "\n")
  
  # Print summary statistics
  valid_results <- results_df[!is.na(results_df$n_cells), ]
  if(nrow(valid_results) > 0) {
    cat("Successfully processed samples:", nrow(valid_results), "\n")
    cat("Total cells across all samples:", sum(valid_results$n_cells), "\n")
    cat("Average cells per sample:", round(mean(valid_results$n_cells), 1), "\n")
    cat("Range of cells per sample:", min(valid_results$n_cells), "-", max(valid_results$n_cells), "\n")
  }
  
  failed_samples <- nrow(results_df) - nrow(valid_results)
  if(failed_samples > 0) {
    cat("Failed to process samples:", failed_samples, "\n")
  }
  
  # Print the results table
  cat("\nDetailed results:\n")
  print(results_df)
  
} else {
  cat("No results to write.\n")
}

cat("\nEvaluation complete!\n")
