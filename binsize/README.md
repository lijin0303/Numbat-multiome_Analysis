# Bin Size Screening Analysis

This directory contains a systematic evaluation of different genomic bin sizes for Numbat-multiome CNV inference, focusing on chromosome 7 data from patient A (patA).

## Overview

The goal of this analysis is to determine the optimal genomic bin size for single-cell multi-omics CNV detection by comparing performance across different bin sizes. This systematic screening helps identify the trade-off between resolution and computational efficiency.

## Bin Sizes Tested

The analysis evaluates the following bin sizes:
- **20kb** - High resolution, fine-grained analysis
- **50kb** - Medium-high resolution
- **80kb** - Medium resolution
- **100kb** - Standard resolution
- **200kb** - Medium-low resolution
- **300kb** - Low resolution
- **500kb** - Very low resolution, coarse-grained analysis

## Directory Structure

```
binsize/
├── inst/                           # Core analysis scripts
│   ├── get_gene_binned_intersections.R  # Gene-to-bin mapping
│   ├── get_binned_rna.R            # RNA count binning
│   ├── get_binned_atac.R           # ATAC count binning
│   ├── input_prep.R                # Data preparation utilities
│   └── run_numbat_multiome.R       # Main Numbat execution
├── patA_bin{size}kb_inputs/        # Input data for each bin size
│   ├── bin{size}kb_chr7.rds        # Genomic bin definitions
│   ├── bin{size}kb_gene2bin_map.csv # Gene-to-bin mapping
│   ├── bin{size}kb_rna_bin.rds     # Binned RNA counts
│   ├── bin{size}kb_atac_bin.rds    # Binned ATAC counts
│   ├── bin{size}kb_comb_bincnt.rds # Combined RNA+ATAC counts
│   ├── lambdas_RNA_bincnt_bin{size}kb.rds    # RNA reference
│   ├── lambdas_ATAC_bincnt_bin{size}kb.rds   # ATAC reference
│   └── lambdas_comb_bincnt_bin{size}kb.rds   # Combined reference
├── patA_bin{size}kb_outputs/       # Numbat results for each bin size
│   ├── log.txt                     # Execution log
│   ├── clones_1.rds, clones_2.rds  # Clone assignments
│   ├── segs_consensus_1.tsv, segs_consensus_2.tsv # CNV segments
│   ├── bulk_clones_*.png           # Clone visualization
│   ├── panel_*.png                 # CNV panel plots
│   └── tree_*.rds                  # Phylogenetic trees
├── bin_generate.sh                 # Script to generate bin definitions
├── bin_runtime.sh                  # Script to extract runtime data
├── runtime.tsv                     # Runtime comparison results
└── [data files]                    # Raw data and annotations
```

## Analysis Pipeline

### 1. Bin Generation (`bin_generate.sh`)
```bash
# Generate genomic bins for different sizes
Rscript -e '
library(GenomicRanges)
chr7len <- as.integer(numbat::chrom_sizes_hg38[7,"size"])
seql <- c(chr7=chr7len)
bsizes <- c(20e3, 50e3, 80e3, 100e3, 200e3, 300e3, 500e3)
for (b in bsizes) {
  gr <- tileGenome(seqlengths=seql, tilewidth=b, cut.last.tile.in.chrom=TRUE)
  saveRDS(gr, paste0("binsize/bin", b/1e3, "kb_chr7.rds"))
}'
```

### 2. Data Processing for Each Bin Size
For each bin size, the pipeline:
1. **Gene-to-bin mapping**: Maps genes to genomic bins
2. **RNA binning**: Aggregates RNA counts into bins
3. **ATAC binning**: Aggregates ATAC fragment counts into bins
4. **Combination**: Merges RNA and ATAC data
5. **Reference generation**: Creates normal cell references

### 3. Numbat Execution (`numbat_chr7.sh`)
```bash
# Run Numbat for each bin size
for prefix in bin20kb bin50kb bin80kb bin100kb bin200kb bin300kb bin500kb; do
  bash binsize/numbat_chr7.sh $prefix
done
```

### 4. Runtime Analysis (`bin_runtime.sh`)
Extracts execution times from log files to compare computational efficiency across bin sizes.

## Results

### Runtime Performance
Based on `runtime.tsv`, execution times decrease with larger bin sizes:

| Bin Size | Runtime (minutes) |
|----------|------------------|
| 50kb     | 21               |
| 80kb     | 17               |
| 100kb    | 15               |
| 200kb    | 11               |
| 300kb    | 9                |
| 500kb    | 6                |

**Key Observation**: Larger bin sizes (300-500kb) are ~2-3x faster than smaller bins (50-100kb).

### Output Files for Each Bin Size

Each `patA_bin{size}kb_outputs/` directory contains:

#### **Core Results:**
- `clones_1.rds`, `clones_2.rds` - Cell-to-clone assignments
- `segs_consensus_1.tsv`, `segs_consensus_2.tsv` - CNV segment calls
- `geno_1.tsv`, `geno_2.tsv` - Genotype matrices

#### **Visualizations:**
- `bulk_clones_*.png` - Clone composition plots
- `panel_*.png` - CNV panel visualizations
- `exp_roll_clust.png` - Expression clustering

#### **Phylogenetic Analysis:**
- `tree_final_*.rds` - Final phylogenetic trees
- `treeML_*.rds` - Maximum likelihood trees
- `treeUPGMA_*.rds` - UPGMA trees

## Key Findings

### **Resolution vs. Performance Trade-off:**
- **Small bins (20-50kb)**: Higher resolution, longer runtime
- **Large bins (300-500kb)**: Lower resolution, faster runtime
- **Optimal range**: 100-200kb provides good balance

### **Data Quality Considerations:**
- Smaller bins may capture fine-scale CNVs but require more cells
- Larger bins are more robust to sparse data but may miss small CNVs
- 100kb appears to be a good compromise for most applications

## Usage

### **Running the Complete Screen:**
```bash
# 1. Generate bins
bash binsize/bin_generate.sh

# 2. Run analysis for all bin sizes
bash bin_screen.sh

# 3. Extract runtime data
bash binsize/bin_runtime.sh
```

### **Running Individual Bin Size:**
```bash
# Run specific bin size (e.g., 100kb)
bash binsize/numbat_chr7.sh bin100kb
```

### **Analyzing Results:**
```r
# Load results for comparison
library(dplyr)
results <- list()
for(size in c("20kb", "50kb", "80kb", "100kb", "200kb", "300kb", "500kb")) {
  results[[size]] <- readRDS(paste0("binsize/patA_bin", size, "_outputs/clones_2.rds"))
}
```

## Recommendations

1. **For high-resolution studies**: Use 50-100kb bins
2. **For large-scale screening**: Use 200-300kb bins
3. **For computational efficiency**: Use 300-500kb bins
4. **Default recommendation**: 100kb provides optimal balance

## Dependencies

- **R packages**: GenomicRanges, numbat, dplyr, ggplot2
- **Data**: scRNA-seq and scATAC-seq data for chromosome 7
- **Reference**: Normal cell annotations for reference generation

## Notes

- Analysis focuses on chromosome 7 for computational efficiency
- Results may vary with different datasets and cell types
- Consider cell count and data quality when choosing bin size
- Runtime estimates are based on the specific hardware configuration used 