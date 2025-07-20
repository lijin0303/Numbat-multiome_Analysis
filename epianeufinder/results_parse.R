#!/usr/bin/env Rscript

# epiAneufinder Results Parser
# Analyzes CNV calls to identify bins with >50 cells having CNV != 1
# and generates summary tables of deletions (0) and amplifications (2)

require(optparse)
require(data.table)
options(stringsAsFactors = FALSE)

# Command line options
option_list <- list(
    make_option("--sample_dir", type="character", default=NULL,
              help="Path to sample directory containing epiAneufinder_results/",
               metavar="character"),
  make_option(c("--output_dir"), type="character", default=NULL,
              help="Output directory for results [default= %default]", metavar="character"),
  make_option(c("--min_cells"), type="integer", default=50,
              help="Minimum number of cells with CNV != 1 to flag a bin [default= %default]", metavar="integer"),
  make_option(c("--all_samples"), action="store_true", default=FALSE,
              help="Process all sample directories in current folder")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Debug: Print parsed arguments
cat("Parsed arguments:\n")
cat("  sample_dir:", opt$sample_dir, "\n")
cat("  all_samples:", opt$all_samples, "\n")
cat("  Is sample_dir NULL?", is.null(opt$sample_dir), "\n")
if(is.null(opt$output_dir)) {
  opt$output_dir <- opt$sample_dir
}
if(is.null(opt$sample_dir) && !opt$all_samples){
  print_help(opt_parser)
  stop("At least one of --sample_dir or --all_samples must be specified.", call.=FALSE)
}

# Function to process a single sample using the correct epiAneufinder structure
process_sample <- function(sample_path, min_cells = 50, output_dir = ".") {
  
  cat("Processing sample:", sample_path, "\n")
  
  # Extract sample name
  sample_name <- basename(sample_path)
  
  # Check if output file already exists
  cnv_calls_file <- file.path(output_dir, paste0(sample_name, "_cnv_calls.tsv"))
  if(file.exists(cnv_calls_file)) {
    cat("  CNV calls file already exists, skipping:", cnv_calls_file, "\n")
    
    # Return basic summary stats without reprocessing
    summary_stats <- list(
      sample = sample_name,
      total_bins = NA,
      total_cells = NA,
      flagged_bins = NA,
      total_deletions_in_flagged = NA,
      total_amplifications_in_flagged = NA,
      cnv_calls_deletions = NA,
      cnv_calls_amplifications = NA,
      min_cells_threshold = min_cells,
      skipped = TRUE
    )
    return(summary_stats)
  }
  
  # Check for required files
  results_dir <- file.path(sample_path, "epiAneufinder_results")
  cnv_file <- file.path(results_dir, "cnv_calls.rds")
  results_file <- file.path(results_dir, "results_table.tsv")
  
  if(!dir.exists(results_dir)) {
    cat("Warning: Results directory not found for", sample_name, "\n")
    return(NULL)
  }
  
  if(!file.exists(cnv_file) || !file.exists(results_file)) {
    cat("  Warning: Required files not found for", sample_name, "\n")
    cat("    CNV calls file exists:", file.exists(cnv_file), "\n")
    cat("    Results table exists:", file.exists(results_file), "\n")
    return(NULL)
  }
  
  # Load bin information from results_table.tsv
  cat("  Loading bin information from results table...\n")
  tryCatch({
    # Read the results table, drop first column (likely row numbers), take first 3 columns
    binInfo <- data.table::fread(results_file, header = F, stringsAsFactors = FALSE)
    # Drop first column and take columns 1:3 as seqnames, start, end
    binInfo <- binInfo[, 2:4]  # Skip first column, take next 3
    colnames(binInfo) <- c("seqnames", "start", "end")
    binInfo$bin_id <- paste0("bin", 1:nrow(binInfo))
    cat("    Loaded", nrow(binInfo), "bins\n")
  }, error = function(e) {
    cat("  Error reading results table:", e$message, "\n")
    return(NULL)
  })
  
  # Load CNV calls
  cat("  Loading CNV calls from RDS file...\n")
  tryCatch({
    CNVcall <- readRDS(cnv_file)
    cat("    Found CNV calls for", length(CNVcall), "cells\n")
    
    # Convert list to matrix (cells x bins)
    CNVmat <- do.call(rbind, CNVcall)
    colnames(CNVmat) <- binInfo$bin_id
    rownames(CNVmat) <- gsub("cell-", "", names(CNVcall))
    
    # Transpose to get bins x cells matrix for easier analysis
    cnv_matrix <- t(CNVmat)
    cat("    CNV matrix dimensions:", nrow(cnv_matrix), "bins x", ncol(cnv_matrix), "cells\n")
    
  }, error = function(e) {
    cat("  Error processing CNV calls:", e$message, "\n")
    return(NULL)
  })
  
  if(!exists("cnv_matrix")) {
    cat("  Warning: Could not create CNV matrix for", sample_name, "\n")
    return(NULL)
  }
  
  # Analyze CNVs per bin
  cnv_summary <- data.frame(
    bin_id = binInfo$bin_id,
    seqnames = binInfo$seqnames,
    start = binInfo$start,
    end = binInfo$end,
    total_cells = ncol(cnv_matrix),
    cells_with_cnv_neq_1 = 0,
    deletions_0 = 0,
    normal_1 = 0,
    amplifications_2 = 0,
    other_cnv = 0,
    stringsAsFactors = FALSE
  )
  
  # Count CNVs for each bin
# Use apply for efficient row-wise computation with \(x) syntax
cnv_summary$cells_with_cnv_neq_1 <- apply(cnv_matrix, 1, \(x) sum(!is.na(x) & x != 1))
cnv_summary$deletions_0         <- apply(cnv_matrix, 1, \(x) sum(x == 0, na.rm = TRUE))
cnv_summary$normal_1            <- apply(cnv_matrix, 1, \(x) sum(x == 1, na.rm = TRUE))
cnv_summary$amplifications_2    <- apply(cnv_matrix, 1, \(x) sum(x == 2, na.rm = TRUE))
cnv_summary$other_cnv           <- apply(cnv_matrix, 1, \(x) sum(!is.na(x) & !(x %in% c(0, 1, 2))))
  
  # Filter bins with more than min_cells having CNV != 1
  flagged_bins <- cnv_summary[cnv_summary$cells_with_cnv_neq_1 >= min_cells, ]
  
  cat("  Total bins:", nrow(cnv_summary), "\n")
  cat("  Flagged bins (>= ", min_cells, " cells with CNV != 1):", nrow(flagged_bins), "\n")
  
  # Generate summary statistics
  if(nrow(flagged_bins) > 0) {
    cat("  Summary of flagged bins:\n")
    cat("    Total deletions (CNV=0):", sum(flagged_bins$deletions_0), "\n")
    cat("    Total amplifications (CNV=2):", sum(flagged_bins$amplifications_2), "\n")
    cat("    Total other CNVs:", sum(flagged_bins$other_cnv), "\n")
  }
  
  # Create CNV calls file with amp/del format
  cnv_calls_list <- list()
  
  # Process deletions (CNV=0)
  del_bins <- cnv_summary[cnv_summary$deletions_0 >= min_cells, ]
  if(nrow(del_bins) > 0) {
    del_df <- data.frame(
      seqnames = del_bins$seqnames,
      start = del_bins$start,
      end = del_bins$end,
      cnv = "del",
      n_cells = del_bins$deletions_0,
      stringsAsFactors = FALSE
    )
    cnv_calls_list[["deletions"]] <- del_df
  }
  
  # Process amplifications (CNV=2)
  amp_bins <- cnv_summary[cnv_summary$amplifications_2 >= min_cells, ]
  if(nrow(amp_bins) > 0) {
    amp_df <- data.frame(
      seqnames = amp_bins$seqnames,
      start = amp_bins$start,
      end = amp_bins$end,
      cnv = "amp",
      n_cells = amp_bins$amplifications_2,
      stringsAsFactors = FALSE
    )
    cnv_calls_list[["amplifications"]] <- amp_df
  }
  
  # Combine all CNV calls
  if(length(cnv_calls_list) > 0) {
    cnv_calls_combined <- do.call(rbind, cnv_calls_list)
    rownames(cnv_calls_combined) <- NULL
    
    # Sort by chromosome and position
    cnv_calls_combined <- cnv_calls_combined[order(cnv_calls_combined$seqnames, 
                                                   cnv_calls_combined$start), ]
  } else {
    cnv_calls_combined <- data.frame(
      seqnames = character(0),
      start = integer(0),
      end = integer(0),
      cnv = character(0),
      n_cells = integer(0),
      stringsAsFactors = FALSE
    )
  }
  
  # Save results
  output_file <- file.path(output_dir, paste0(sample_name, "_cnv_analysis.tsv"))
  flagged_output_file <- file.path(output_dir, paste0(sample_name, "_flagged_bins.tsv"))
  cnv_calls_file <- file.path(output_dir, paste0(sample_name, "_cnv_calls.tsv"))
  
  # Write all bins summary
  write.table(cnv_summary, file = output_file, sep = "\t", quote = FALSE, 
              row.names = FALSE, col.names = TRUE)
  cat("  All bins summary saved to:", output_file, "\n")
  
  # Write flagged bins only
  if(nrow(flagged_bins) > 0) {
    write.table(flagged_bins, file = flagged_output_file, sep = "\t", quote = FALSE, 
                row.names = FALSE, col.names = TRUE)
    cat("  Flagged bins saved to:", flagged_output_file, "\n")
  }
  
  # Write CNV calls in requested format
  write.table(cnv_calls_combined, file = cnv_calls_file, sep = "\t", quote = FALSE, 
              row.names = FALSE, col.names = TRUE)
  cat("  CNV calls (amp/del format) saved to:", cnv_calls_file, "\n")
  cat("    Found", nrow(cnv_calls_combined), "CNV calls (", 
      sum(cnv_calls_combined$cnv == "del"), "deletions,", 
      sum(cnv_calls_combined$cnv == "amp"), "amplifications)\n")
  
  # Create a summary report
  summary_stats <- list(
    sample = sample_name,
    total_bins = nrow(cnv_summary),
    total_cells = ncol(cnv_matrix),
    flagged_bins = nrow(flagged_bins),
    total_deletions_in_flagged = sum(flagged_bins$deletions_0),
    total_amplifications_in_flagged = sum(flagged_bins$amplifications_2),
    cnv_calls_deletions = sum(cnv_calls_combined$cnv == "del"),
    cnv_calls_amplifications = sum(cnv_calls_combined$cnv == "amp"),
    min_cells_threshold = min_cells,
    skipped = FALSE
  )
  
  return(summary_stats)
}

# Main execution
cat("=== epiAneufinder Results Parser ===\n")
cat("Minimum cells threshold:", opt$min_cells, "\n")
cat("Output directory:", opt$output_dir, "\n\n")

# Create output directory if it doesn't exist
if(!is.null(opt$output_dir)){
 if(!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE)
}
}


all_summaries <- list()

if(opt$all_samples) {
  # Process all sample directories
  cat("Processing all sample directories...\n")
  
  # Find all directories starting with "epianeufinder_"
  sample_dirs <- list.dirs(".", recursive = FALSE, full.names = TRUE)
  sample_dirs <- sample_dirs[grepl("epianeufinder_", basename(sample_dirs))]
  
  if(length(sample_dirs) == 0) {
    stop("No sample directories found with pattern 'epianeufinder_*'")
  }
  
  cat("Found", length(sample_dirs), "sample directories\n")
  
  for(sample_dir in sample_dirs) {
    summary_stats <- process_sample(sample_dir, opt$min_cells, sample_dir)
    if(!is.null(summary_stats)) {
      all_summaries[[length(all_summaries) + 1]] <- summary_stats
    }
  }
  
} else {
  # Process single sample
  summary_stats <- process_sample(opt$sample_dir, opt$min_cells, opt$output_dir)
  if(!is.null(summary_stats)) {
    all_summaries[[1]] <- summary_stats
  }
}

# Generate overall summary
if(length(all_summaries) > 0) {
  cat("\n=== Overall Summary ===\n")
  
  # Convert to data frame
  summary_df <- do.call(rbind, lapply(all_summaries, function(x) {
    data.frame(
      sample = x$sample,
      total_bins = ifelse(x$skipped, NA, x$total_bins),
      total_cells = ifelse(x$skipped, NA, x$total_cells),
      flagged_bins = ifelse(x$skipped, NA, x$flagged_bins),
      total_deletions = ifelse(x$skipped, NA, x$total_deletions_in_flagged),
      total_amplifications = ifelse(x$skipped, NA, x$total_amplifications_in_flagged),
      cnv_calls_deletions = ifelse(x$skipped, NA, x$cnv_calls_deletions),
      cnv_calls_amplifications = ifelse(x$skipped, NA, x$cnv_calls_amplifications),
      min_cells_threshold = x$min_cells_threshold,
      skipped = x$skipped,
      stringsAsFactors = FALSE
    )
  }))
  
  # Save overall summary
  overall_summary_file <- file.path(opt$output_dir, "overall_cnv_summary.tsv")
  write.table(summary_df, file = overall_summary_file, sep = "\t", quote = FALSE, 
              row.names = FALSE, col.names = TRUE)
  
  cat("Overall summary saved to:", overall_summary_file, "\n")
  
  # Print summary to console
  print(summary_df)
  
  # Count processed vs skipped samples
  processed_samples <- sum(!summary_df$skipped, na.rm = TRUE)
  skipped_samples <- sum(summary_df$skipped, na.rm = TRUE)
  
  cat("\nSample processing summary:\n")
  cat("  Total samples:", nrow(summary_df), "\n")
  cat("  Processed samples:", processed_samples, "\n")
  cat("  Skipped samples (already processed):", skipped_samples, "\n")
  
  if(processed_samples > 0) {
    # Only calculate totals for processed samples
    processed_df <- summary_df[!summary_df$skipped, ]
    cat("\nTotal across processed samples:\n")
    cat("  Total flagged bins:", sum(processed_df$flagged_bins, na.rm = TRUE), "\n")
    cat("  Total deletions in flagged bins:", sum(processed_df$total_deletions, na.rm = TRUE), "\n")
    cat("  Total amplifications in flagged bins:", sum(processed_df$total_amplifications, na.rm = TRUE), "\n")
    cat("  Total CNV deletion calls:", sum(processed_df$cnv_calls_deletions, na.rm = TRUE), "\n")
    cat("  Total CNV amplification calls:", sum(processed_df$cnv_calls_amplifications, na.rm = TRUE), "\n")
  }
}

cat("\nAnalysis complete!\n")
