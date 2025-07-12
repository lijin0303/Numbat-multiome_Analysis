#!/usr/bin/env bash
set -euo pipefail

prefix="bin100kb"
binGR="bin100kb_chr7.rds"        # ← your custom GRanges RDS
gtfF="hg38"                # or path/to/gtf_hg38.gtf

# # 1. Generate gene-to-bin map
# Rscript benchmark/get_gene_binned_intersections.R \
#   --numbatGTFname ${gtfF} \
#   --binGR          benchmark/${binGR} \
#   --outfile        benchmark/${prefix}_gene2bin_map.csv

# 2. Prepare binned inputs for your sample(s)
# sample="patA" 

# # 2a. RNA bins
# Rscript benchmark/get_binned_rna.R \
#   --rnaCountsFile     benchmark/${sample}.rds \
#   --outFile           benchmark/${prefix}_rna_bin.rds \
#   --barcodesKeep      benchmark/test_RNA.txt \
#   --geneBinMapCSVFile benchmark/${prefix}_gene2bin_map.csv

# # 2b. ATAC bins
# Rscript benchmark/get_binned_atac.R \
#   --CB      benchmark/test_ATAC.txt \
#   --frag    benchmark/chr7.fragments.tsv.gz \
#   --binGR   benchmark/${binGR} \
#   --outFile benchmark/${prefix}_atac_bin.rds

# # 2c. Combine RNA & ATAC bins into one matrix
# Rscript -e "
#   library(glue);
#   source('benchmark/input_prep.R');
#   saveRDS(
#     binCnt_union(
#       c(
#         glue('benchmark/${prefix}_rna_bin.rds'),
#         glue('benchmark/${prefix}_atac_bin.rds')
#       ),
#       seed  = 123,
#       maxCB = 10000
#     ),
#     glue('benchmark/${prefix}_comb_bincnt.rds')
#   )
# "

# # 3. Build references from your normal samples
# # 3a. Per-sample binned references
# Rscript benchmark/get_binned_rna.R \
#     --rnaCountsFile     benchmark/${sample}.rds \
#     --outFile           benchmark/lambdas_RNA_bincnt_${prefix}.rds \
#     --barcodesKeep      benchmark/normal_RNA_annot.tsv \
#     --geneBinMapCSVFile benchmark/${prefix}_gene2bin_map.csv \
#     --generateAggRef

# Rscript benchmark/get_binned_atac.R \
#     --frag    benchmark/chr7.fragments.tsv.gz \
#     --CB      benchmark/normal_ATAC_annot.tsv \
#     --binGR   benchmark/${binGR} \
#     --outFile benchmark/lambdas_ATAC_bincnt_${prefix}.rds \
#     --generateAggRef

# 3b. Merge per-sample refs into combined Reference
# Rscript -e "
#   library(glue);
#   source('benchmark/input_prep.R');
#   saveRDS(
#     binCnt_union(
#       c(
#         glue('benchmark/lambdas_RNA_bincnt_${prefix}.rds'),
#         glue('benchmark/lambdas_ATAC_bincnt_${prefix}.rds')
#       ),
#       seed  = 123,
#       maxCB = 10000
#     ),
#     glue('benchmark/lambdas_comb_bincnt_${prefix}.rds')
#   )
# "


# parL="par_numbatm.rds" # a list of any run_numbat parameters you would like to optimize
Rscript benchmark/run_numbat_multiome.R  \
        --countmat benchmark/benchmark/${prefix}_comb_bincnt.rds \
        --alleledf benchmark/${sample}_allele_counts_chr7.tsv.gz \
        --out_dir benchmark/${sample}_${prefix}/ \
        --ref  benchmark/lambdas_comb_bincnt_${prefix}.rds \
        --gtf  benchmark/${binGR}\
        --parL benchmark/par_comb_bincnt.rds
