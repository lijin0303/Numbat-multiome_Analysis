Rscript -e 'library(GenomicRanges)
# get chr7 length
chr7len <- as.integer(numbat::chrom_sizes_hg38[7,"size"])
# set up a named seqlengths vector
seql   <- c(chr7=chr7len);
# define your bin sizes
# bsizes <- c(100e3, 200e3, 300e3, 500e3);
bsizes <- c(20e3, 50e3, 80e3);
# tile and save
for (b in bsizes) {
  gr <- tileGenome(seqlengths=seql, tilewidth=b, cut.last.tile.in.chrom=TRUE);
  saveRDS(gr, paste0("binsize/bin", b/1e3, "kb_chr7.rds"));
}'