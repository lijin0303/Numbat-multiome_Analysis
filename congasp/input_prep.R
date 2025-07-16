library(Matrix)
library(readr)
library(optparse)

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input RDS file path", metavar="file")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input)) {
  stop("Please provide an input RDS file with -i <file>")
}

input_rds <- opt$input
basename <- tools::file_path_sans_ext(basename(input_rds))
outdir <- file.path("congasp/inputData", basename)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# 1. Read the bin count matrix from RDS
bin_counts <- readRDS(input_rds)

# 2. Convert sparse matrix to dense matrix, then to data frame
bin_counts_dense <- as.matrix(bin_counts)
bin_counts_df <- as.data.frame(bin_counts_dense)

# 3. Write the count matrix as CSV (rows = bins, columns = cells)
write_csv(bin_counts_df, file.path(outdir, "bin_counts_for_python.csv"))

# 4. Write the bin/segment names (if rownames are bins)
if (!is.null(rownames(bin_counts))) {
  bins <- rownames(bin_counts)
} else {
  bins <- paste0("bin", seq_len(nrow(bin_counts)))
}
write_csv(data.frame(bin=bins), file.path(outdir, "bin_segments.csv"))

# 5. Write the ploidy/prior file (all 2s, one per bin)
ploidy <- rep(2, length(bins))
write_csv(data.frame(bin=bins, ploidy=ploidy), file.path(outdir, "bin_ploidy.csv"))