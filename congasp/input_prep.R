library(Matrix)
library(readr)

# 1. Read the bin count matrix from RDS
bin_counts <- readRDS("benchmark/patA_bin200kb_inputs/bin200kb_comb_bincnt.rds")  

# 2. Convert sparse matrix to dense matrix, then to data frame
bin_counts_dense <- as.matrix(bin_counts)
bin_counts_df <- as.data.frame(bin_counts_dense)

# 3. Write the count matrix as CSV (rows = bins, columns = cells)
write_csv(bin_counts_df, "congasp/bin_counts_for_python.csv")

# 4. Write the bin/segment names (if rownames are bins)
if (!is.null(rownames(bin_counts))) {
  bins <- rownames(bin_counts)
} else {
  bins <- paste0("bin", seq_len(nrow(bin_counts)))
}
write_csv(data.frame(bin=bins), "congasp/bin_segments.csv")

# 5. Write the ploidy/prior file (all 2s, one per bin)
ploidy <- rep(2, length(bins))
write_csv(data.frame(bin=bins, ploidy=ploidy), "congasp/bin_ploidy.csv")