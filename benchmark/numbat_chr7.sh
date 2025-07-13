#!/usr/bin/env bash
set -euo pipefail

# 0. Config

workdir="benchmark"
rsrc="benchmark/inst"
prefix="bin100kb"
binGR="${prefix}_chr7.rds"    # custom GRanges RDS under $rsrc
gtfF="hg38"                   # or path/to/gtf_hg38.gtf
sample="patA"

# 1. Make dirs
inputdir="${workdir}/${sample}_${prefix}_inputs"
outdir="${workdir}/${sample}_${prefix}_outputs"
mkdir -p "$inputdir" "$outdir"

# 3. Gene→bin map (prefix-specific → $inputdir)
Rscript "${rsrc}/get_gene_binned_intersections.R" \
  --numbatGTFname "${gtfF}" \
  --binGR          "${rsrc}/${binGR}" \
  --outfile        "${inputdir}/${prefix}_gene2bin_map.csv"

# 4. Binned inputs for your sample
## 4a. RNA bins → $inputdir
Rscript "${rsrc}/get_binned_rna.R" \
  --rnaCountsFile     "${workdir}/${sample}.rds" \
  --outFile           "${inputdir}/${prefix}_rna_bin.rds" \
  --barcodesKeep      "${workdir}/test_RNA.txt" \
  --geneBinMapCSVFile "${inputdir}/${prefix}_gene2bin_map.csv"

## 4b. ATAC bins → $inputdir
Rscript "${rsrc}/get_binned_atac.R" \
  --CB      "${workdir}/test_ATAC.txt" \
  --frag    "${workdir}/chr7.fragments.tsv.gz" \
  --binGR   "${rsrc}/${binGR}" \
  --outFile "${inputdir}/${prefix}_atac_bin.rds"

## 4c. Combine RNA & ATAC → $inputdir
Rscript -e "
  library(glue);
  source('${rsrc}/input_prep.R');
  saveRDS(
    binCnt_union(
      c(
        glue('${inputdir}/${prefix}_rna_bin.rds'),
        glue('${inputdir}/${prefix}_atac_bin.rds')
      ),
      seed  = 123,
      maxCB = 10000
    ),
    glue('${inputdir}/${prefix}_comb_bincnt.rds')
  )
"

# 5. Build refs from normal cells
## 5a. RNA reference → $inputdir
Rscript "${rsrc}/get_binned_rna.R" \
  --rnaCountsFile     "${workdir}/${sample}.rds" \
  --outFile           "${inputdir}/lambdas_RNA_bincnt_${prefix}.rds" \
  --barcodesKeep      "${workdir}/normal_RNA_annot.tsv" \
  --geneBinMapCSVFile "${inputdir}/${prefix}_gene2bin_map.csv" \
  --generateAggRef

## 5b. ATAC reference → $inputdir
Rscript "${rsrc}/get_binned_atac.R" \
  --frag    "${workdir}/chr7.fragments.tsv.gz" \
  --CB      "${workdir}/normal_ATAC_annot.tsv" \
  --binGR   "${rsrc}/${binGR}" \
  --outFile "${inputdir}/lambdas_ATAC_bincnt_${prefix}.rds" \
  --generateAggRef

## 5c. Merge RNA+ATAC refs → $inputdir
Rscript -e "
  library(glue);
  source('${rsrc}/input_prep.R');
  saveRDS(
    binCnt_union(
      c(
        glue('${inputdir}/lambdas_RNA_bincnt_${prefix}.rds'),
        glue('${inputdir}/lambdas_ATAC_bincnt_${prefix}.rds')
      ),
      seed  = 123,
      maxCB = 10000
    ),
    glue('${inputdir}/lambdas_comb_bincnt_${prefix}.rds')
  )
"

# 6. Run Numbat → final outputs in $outdir
Rscript "${rsrc}/run_numbat_multiome.R" \
  --countmat "${inputdir}/${prefix}_comb_bincnt.rds" \
  --alleledf "${workdir}/${sample}_allele_counts_chr7.tsv.gz" \
  --out_dir  "${outdir}" \
  --ref      "${inputdir}/lambdas_comb_bincnt_${prefix}.rds" \
  --gtf      "${rsrc}/${binGR}" \
  --parL     "${workdir}/par_comb_bincnt.rds"

echo "Done!
- Non-prefix inputs (test_RNA.txt, test_ATAC.txt, sample.rds, fragment & allele files, parL) → ${workdir}
- Prefix-specific intermediates (gene2bin map, binned counts, refs) → ${inputdir}
- Final NumBat outputs → ${outdir}"