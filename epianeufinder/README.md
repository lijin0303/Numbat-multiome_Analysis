# epiAneufinder Benchmarking Setup

This folder contains the setup and scripts for benchmarking scATAC-seq data-based CNV inference using epiAneufinder. The benchmarking uses the same raw `fragments.tsv.gz` files that were used to generate our cell-by-bin matrix, ensuring fair comparison between methods.

## Overview

epiAneufinder is a tool for detecting copy number variations (CNVs) from single-cell ATAC-seq data. This benchmarking setup evaluates its performance using:

- **Same input data**: Identical `fragments.tsv.gz` files used in our primary analysis
- **Matched cell barcodes**: Same cell barcode (CB) lists to ensure comparable cell populations
- **Systematic evaluation**: Automated processing with timing and performance metrics

## Key Components

### 1. Docker Image with Latest epiAneufinder

#### Why We Built a Custom Image
- **Latest version**: Uses the most recent epiAneufinder from GitHub commit `6969d6c4926acb707087f7a03dcfabca7edd5cfe`
- **Simplified installation**: Eliminates dependency issues and ensures reproducible environment
- **Consistent environment**: All analyses run in the same containerized environment

#### Image Building Process
The custom Docker image is built using:

```bash
# Build the updated image
./build_and_push.sh
```

**Dockerfile highlights:**
- Base image: `gcr.io/broad-getzlab-workflows/scatac_tool:v41`
- Updated epiAneufinder from: `colomemaria/epiAneufinder@6969d6c`
- Additional dependencies for the latest version
- Optimized for reproducible CNV analysis

### 2. Input Configuration

#### File Structure
```
epianeufinder_inputs.txt    # Sample definitions with GCS paths
tmp_CB/                     # Cell barcode files per sample
├── {sample}_CB.txt         # Cell barcodes for each sample
└── ...
```

#### Input Format (`epianeufinder_inputs.txt`)
Tab-separated file with columns:
- **Sample ID**: Unique identifier
- **Fragment file**: GCS path to `fragments.tsv.gz`
- **Index file**: GCS path to `fragments.tsv.gz.tbi`
- **Cell barcodes**: GCS path to cell barcode file

#### Cell Barcode Matching
- Cell barcode files in `tmp_CB/` contain the exact same cells used in our primary analysis
- Format: One barcode per line, no header
- Ensures direct comparability between epiAneufinder and other methods

### 3. Automated Processing Pipeline

#### Main Execution Script: `epiAneufinder.sh`

**Key Features:**
- Processes samples from `epianeufinder_inputs.txt`
- Automatic file management (downloads, caching)
- Cell barcode filtering
- Performance optimization (uses 2/3 available cores)
- Comprehensive timing logging

**Usage:**
```bash
./epiAneufinder.sh
```

**What it does:**
1. **File Management**: Downloads fragment files from GCS if not present locally
2. **Cell Filtering**: Applies sample-specific cell barcode lists
3. **CNV Detection**: Runs epiAneufinder with optimized parameters
4. **Result Upload**: Uploads CNV calls to GCS bucket
5. **Performance Tracking**: Records detailed timing information

#### Processing Parameters
- **Window size**: 100kb bins
- **Genome**: hg38
- **Excluded chromosomes**: chrX, chrY, chrM
- **Karyogram plotting**: Disabled (for speed)
- **Cell filtering**: Applied via CB files

### 4. Performance Evaluation

#### Timing Analysis (`epiAneufinder_timing.log`)
Automatically generated log with columns:
- **sample**: Sample identifier
- **start_time**: Processing start timestamp
- **end_time**: Processing completion timestamp  
- **duration_seconds**: Total processing time (numeric)
- **duration_formatted**: Human-readable duration (e.g., "5h 43m 49s")

**Current Performance Metrics:**
| Sample  | Duration | 
|---------|----------|
| pM10114 | 5h 50m   |
| pM10975 | 3h 35m   |
| pM11004 | 5h 43m   |
| pM1835  | 4h 16m   |
| pM9916  | 2h 14m   |

#### Results Analysis
The `results_parse.R` script provides:
- CNV call summarization
- Cell-level statistics
- Bin-level analysis
- Output in standardized format for comparison

### 5. Output Files

#### Per Sample Outputs
```
epianeufinder_{sample}/
├── epiAneufinder_results/
│   ├── cnv_calls.rds           # Raw CNV calls (list format)
│   ├── results_table.tsv       # Bin information and calls
│   └── ...
├── {sample}_cnv_calls.tsv      # Processed CNV calls (amp/del format)
├── {sample}_cnv_analysis.tsv   # Detailed bin analysis
└── {sample}_flagged_bins.tsv   # High-confidence CNV regions
```

#### Summary Outputs
- `epiAneufinder_timing.log`: Performance metrics
- `overall_cnv_summary.tsv`: Cross-sample summary
- `evals_cb.txt`: Cell count validation

### 6. Quality Control and Benchmarking

#### Input Validation
- Verifies fragment file accessibility
- Confirms cell barcode file availability
- Checks for required dependencies

#### Performance Metrics
- **Processing time**: Per-sample timing analysis
- **Cell coverage**: Number of cells successfully processed
- **CNV detection**: Bins with confident CNV calls
- **Resource usage**: Memory and CPU utilization

#### Comparison Framework
Results can be directly compared with other CNV detection methods using:
- Same input fragments
- Identical cell populations
- Standardized output formats
- Consistent quality thresholds

## Usage Instructions

### Prerequisites
- Docker installed and configured
- Access to GCS buckets with fragment files
- Cell barcode files in `tmp_CB/` directory

### Running the Analysis

1. **Prepare inputs:**
   ```bash
   # Ensure cell barcode files are present
   ls tmp_CB/*.txt
   ```

2. **Execute pipeline:**
   ```bash
   ./epiAneufinder.sh
   ```

3. **Monitor progress:**
   ```bash
   tail -f epiAneufinder_timing.log
   ```

4. **Analyze results:**
   ```bash
   Rscript results_parse.R --all_samples
   ```

### Customization Options

- **Minimum cell threshold**: Modify `--min_cells` in results parsing
- **Processing cores**: Adjust core usage in `epiAneufinder.R`
- **Output location**: Change GCS bucket in `epiAneufinder.sh`
- **Window size**: Modify `windowSize` parameter in `epiAneufinder.R`

## Files Description

| File | Purpose |
|------|---------|
| `Dockerfile` | Custom image with latest epiAneufinder |
| `build_and_push.sh` | Image building and deployment script |
| `epiAneufinder.sh` | Main processing pipeline |
| `epiAneufinder.R` | Core epiAneufinder execution script |
| `results_parse.R` | Results analysis and summarization |
| `epiAneufinder_inputs.txt` | Sample and file specifications |
| `epiAneufinder_timing.log` | Performance tracking log |
| `image_load.sh` | Helper for loading Docker environment |

## Benchmarking Considerations

This setup enables fair comparison with other CNV detection methods by:

1. **Identical inputs**: Same fragment files and cell populations
2. **Standardized processing**: Consistent parameters and thresholds  
3. **Comprehensive metrics**: Timing, accuracy, and resource usage
4. **Reproducible environment**: Containerized execution
5. **Automated workflow**: Reduces manual intervention and errors

The timing logs and standardized outputs facilitate systematic evaluation of epiAneufinder's performance characteristics compared to other approaches in our CNV detection benchmark.
