#!/usr/bin/env bash
set -euo pipefail

prefix="bin100kb"
binGR="bin100kb_chr7.rds"        # ← your custom GRanges RDS
gtfF="hg38"                # or path/to/gtf_hg38.gtf

# 1. Generate gene-to-bin map
Rscript binsize/get_gene_binned_intersections.R \
  --numbatGTFname ${gtfF} \
  --binGR          binsize/${binGR} \
  --outfile        binsize/${prefix}_gene2bin_map.csv

# 2. Prepare binned inputs for your sample(s)
sample="patA" 

# 2a. RNA bins
Rscript binsize/get_binned_rna.R \
  --rnaCountsFile     binsize/${sample}.rds \
  --outFile           binsize/${prefix}_rna_bin.rds \
  --barcodesKeep      binsize/test_RNA.txt \
  --geneBinMapCSVFile binsize/${prefix}_gene2bin_map.csv

# 2b. ATAC bins
Rscript binsize/get_binned_atac.R \
    --CB      binsize/test_ATAC.txt \
    --frag    binsize/chr7.fragments.tsv.gz \
    --binGR   binsize/${binGR} \
    --outFile binsize/${prefix}_atac_bin.rds

# 2c. Combine RNA & ATAC bins into one matrix
Rscript -e "
  library(glue);
  source('binsize/input_prep.R');
  saveRDS(
    binCnt_union(
      c(
        glue('binsize/${prefix}_rna_bin.rds'),
        glue('binsize/${prefix}_atac_bin.rds')
      ),
      seed  = 123,
      maxCB = 10000
    ),
    glue('binsize/${prefix}_comb_bincnt.rds')
  )
"

# 3. Build references from your normal samples
# 3a. Per-sample binned references
Rscript binsize/get_binned_rna.R \
    --rnaCountsFile     binsize/${sample}.rds \
    --outFile           binsize/lambdas_RNA_bincnt_${prefix}.rds \
    --barcodesKeep      binsize/normal_RNA_annot.tsv \
    --geneBinMapCSVFile binsize/${prefix}_gene2bin_map.csv \
    --generateAggRef

Rscript binsize/get_binned_atac.R \
    --frag    binsize/chr7.fragments.tsv.gz \
    --CB      binsize/normal_ATAC_annot.tsv \
    --binGR   binsize/${binGR} \
    --outFile binsize/lambdas_ATAC_bincnt_${prefix}.rds \
    --generateAggRef

3b. Merge per-sample refs into combined Reference
Rscript -e "
  library(glue);
  source('binsize/input_prep.R');
  saveRDS(
    binCnt_union(
      c(
        glue('binsize/lambdas_RNA_bincnt_${prefix}.rds'),
        glue('binsize/lambdas_ATAC_bincnt_${prefix}.rds')
      ),
      seed  = 123,
      maxCB = 10000
    ),
    glue('binsize/lambdas_comb_bincnt_${prefix}.rds')
  )
"


Rscript binsize/run_numbat_multiome.R  \
        --countmat binsize/${prefix}_comb_bincnt.rds \
        --alleledf binsize/${sample}_allele_counts_chr7.tsv.gz \
        --out_dir binsize/${sample}_${prefix}_outputs/ \
        --ref  binsize/lambdas_comb_bincnt_${prefix}.rds \
        --gtf  binsize/${binGR}\
        --parL binsize/par_comb_bincnt.rds
