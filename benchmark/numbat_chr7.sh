#!/usr/bin/env bash
set -euo pipefail

# 0. Define your custom binGR and GTF
binGR="my_bins.rds"        # ← your custom GRanges RDS
gtfF="hg38"                # or path/to/gtf_hg38.gtf

# 1. Generate gene-to-bin map
Rscript get_gene_binned_intersections.R \
  --numbatGTFname ${gtfF} \
  --binGR          ${binGR} \
  --outfile        gene2bin_map.csv

# 2. Prepare binned inputs for your sample(s)
sample="Sample1"           # replace with your sample name

# 2a. RNA bins
Rscript get_binned_rna.R \
  --rnaCountsFile     ${sample}.seu.rds \
  --outFile           ${sample}/${sample}_rna_bin.rds \
  --barcodesKeep      ${sample}_barcodes.tsv \
  --geneBinMapCSVFile gene2bin_map.csv

# 2b. ATAC bins
Rscript get_binned_atac.R \
  --CB      ${sample}_atac_barcodes.tsv \
  --frag    ${sample}_fragments.tsv.gz \
  --binGR   ${binGR} \
  --outFile ${sample}/${sample}_atac_bin.rds

# # 2c. Combine RNA & ATAC bins into one matrix
# Rscript -e "
#   library(glue);
#   source('input_prep.R');
#   saveRDS(
#     binCnt(
#       c(
#         glue('${sample}/${sample}_rna_bin.rds'),
#         glue('${sample}/${sample}_atac_bin.rds')
#       ),
#       seed  = 123,
#       maxCB = 10000
#     ),
#     glue('${sample}/${sample}_comb_bincnt.rds')
#   )
# "

# # 3. Build references from your normal samples
# binGR="my_bins.rds"  # same as above
# refsamples=(normal1 normal2)

# # 3a. Per-sample binned references
# for ref in "${refsamples[@]}"; do
#   Rscript get_binned_rna.R \
#     --rnaCountsFile     ${ref}.seu.rds \
#     --outFile           Reference/lambdas_RNA_bincnt_${ref}.rds \
#     --barcodesKeep      ${ref}_barcodes.tsv \
#     --geneBinMapCSVFile gene2bin_map.csv \
#     --generateAggRef

#   Rscript get_binned_atac.R \
#     --CB      ${ref}_atac_barcodes.tsv \
#     --frag    ${ref}_fragments.tsv.gz \
#     --binGR   ${binGR} \
#     --outFile Reference/lambdas_ATAC_bincnt_${ref}.rds \
#     --generateAggRef
# done

# # 3b. Merge per-sample refs into combined Reference
# Rscript -e "
#   library(dplyr);
#   # read in each per-sample ref
#   rna_lst  <- list.files('Reference','lambdas_RNA_bincnt_.*\\.rds$', full.names=TRUE)
#   atac_lst <- list.files('Reference','lambdas_ATAC_bincnt_.*\\.rds$', full.names=TRUE)
#   ref_rna  <- do.call(cbind, lapply(rna_lst, readRDS))
#   ref_atac <- do.call(cbind, lapply(atac_lst, readRDS))
#   shared   <- intersect(rownames(ref_rna), rownames(ref_atac))
#   ref_comb <- cbind(ref_rna[shared,], ref_atac[shared,])
#   saveRDS(ref_comb, 'Reference/lambdas_comb_bincnt.rds')
# "