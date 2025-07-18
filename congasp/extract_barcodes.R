# Extract barcodes from bin_counts_for_python.csv in each subfolder of congasp/inputData
# For patA: exclude barcodes containing 'NIH'
# For others: exclude barcodes containing 'GEX'
# Output: congasp/tmp_CB/<subfolder>_CB.txt (one barcode per line, no header)

input_dir <- "congasp/inputData"
out_dir <- "congasp/tmp_CB"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

subfolders <- list.dirs(input_dir, full.names = FALSE, recursive = FALSE)

for (sub in subfolders) {
  csv_path <- file.path(input_dir, sub, "bin_counts_for_python.csv")
  if (!file.exists(csv_path)) next
  # Read only the header (first line)
  header <- readLines(csv_path, n = 1)
  # Remove any leading/trailing whitespace
  header <- trimws(header)
  # Split by comma
  barcodes <- unlist(strsplit(header, ","))
  # Remove empty strings
  barcodes <- barcodes[barcodes != ""]
  # Filtering
  if (sub == "patA") {
    keep <- !grepl("NIH", barcodes)
  } else {
    keep <- !grepl("GEX", barcodes)
  }
  filtered_barcodes <- barcodes[keep]
  # Write to file
  out_file <- file.path(out_dir, paste0(sub, "_CB.txt"))
  writeLines(filtered_barcodes, out_file)
}
cat("Done. Output written to", out_dir, "\n") 