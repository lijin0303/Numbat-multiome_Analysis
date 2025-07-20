# CONGASp Benchmark

This directory contains the implementation and benchmarking of CONGASp (Copy Number inference from Gene expression and ATAC-Seq Profiles) for multi-omics single-cell copy number variant (CNV) detection.

## Overview

CONGASp is a probabilistic model that integrates RNA-seq and ATAC-seq data from single cells to infer copy number variations. This benchmark evaluates CONGASp's performance across multiple patient samples with different data characteristics.

## Directory Structure

```
congasp/
├── inputData/                    # Processed input data for each patient
│   ├── patA/                     # Patient A data
│   ├── DL3267/                   # Patient DL3267 data
│   ├── pM10109/                  # Patient pM10109 data
│   └── ...                       # Additional patient directories
├── cnv_congasp/                  # CNV segment results from CONGASp
│   ├── patA_cnv_segments_formatted.tsv
│   ├── DL3267_cnv_segments_formatted.tsv
│   └── ...                       # CNV results for each patient
├── tmp_CB/                       # Temporary cell barcode files
├── tmp_rds/                      # Temporary RDS files during processing
├── test_case/                    # Test case with example data and results
├── input_prep.R                  # R script for input data preparation
├── batch_input_prep.sh           # Batch script for processing multiple patients
├── run_congasp_training.py       # Main CONGASp training script
├── run_congasp_training.sh       # Batch script for running CONGASp on all patients
└── extract_barcodes.R            # Script for extracting cell barcodes
```

## Input Data Preparation

### 1. Data Sources
The benchmark uses combined bin count matrices from RNA-seq and ATAC-seq data stored in Google Cloud Storage (`gs://numbat_wrapper/combined_rawcnt/`). Each patient has a corresponding `.rds` file containing:
- **Bin count matrix**: Rows = genomic bins, Columns = cells
- **Cell annotations**: RNA-seq vs ATAC-seq cell identification

### 2. Input Processing Pipeline

#### Step 1: Batch Input Preparation
```bash
bash congasp/batch_input_prep.sh
```

This script:
- Downloads all `.rds` files from Google Cloud Storage
- Processes each file using `input_prep.R`
- Creates standardized input format for CONGASp

#### Step 2: Input Format Conversion (`input_prep.R`)
For each patient, the script generates three CSV files:

1. **`bin_counts_for_python.csv`**: Bin count matrix (bins × cells)
2. **`bin_segments.csv`**: Bin/segment identifiers
3. **`bin_ploidy.csv`**: Prior ploidy information (default: 2 for all bins)

#### Step 3: Cell Type Separation
The pipeline separates RNA-seq and ATAC-seq cells using different separators:
- **patA**: Uses "NIH-A" separator
- **Other patients**: Use "GEX" separator

## CONGASp Model Configuration

### Model Parameters
The CONGASp model is configured with the following parameters:

```python
param_dict = {
    'K': 3,                                    # Number of clusters
    'theta_shape_rna': torch.ones(segments) * 20.0,    # RNA shape parameter
    'theta_rate_rna': torch.ones(segments) * 0.5,      # RNA rate parameter
    'theta_shape_atac': torch.ones(segments) * 10.0,   # ATAC shape parameter
    'theta_rate_atac': torch.ones(segments) * 0.5,     # ATAC rate parameter
    'lambda': 0.5,                             # Regularization parameter
    'multiome': False,                         # Multi-omics mode
    'equal_mixture_weights': True,             # Equal cluster weights
    'likelihood_rna': 'NB',                    # RNA likelihood (Negative Binomial)
    'likelihood_atac': 'NB',                   # ATAC likelihood (Negative Binomial)
    'nb_size_init_rna': torch.ones(segments) * 10.0,   # RNA NB size initialization
    'nb_size_init_atac': torch.ones(segments) * 5.0,   # ATAC NB size initialization
    'hidden_dim': 5                            # Hidden dimension
}
```

### Training Configuration
- **Optimizer**: ClippedAdam
- **Loss Function**: TraceGraph_ELBO
- **Training Steps**: 200 iterations
- **Model**: LatentCategorical

## Running CONGASp

### Single Patient Run
```bash
python3 run_congasp_training.py \
    --separator "NIH-A" \
    --output_dir "/workspace/inputData/patA" \
    --input_dir "/workspace/inputData/patA"
```

### Batch Run (All Patients)
```bash
bash run_congasp_training.sh
```

This script automatically:
- Iterates through all patient directories in `/workspace/inputData/`
- Sets appropriate separators for each patient
- Runs CONGASp training with error handling
- Includes memory management and cleanup

## Outputs

### 1. CONGASp Results (`congas_results.npy`)
Each patient directory contains a NumPy file with learned model parameters:
- Cluster assignments
- Copy number estimates
- Model convergence information

### 2. CNV Segments (`cnv_congasp/*_cnv_segments_formatted.tsv`)
Formatted CNV segment files with columns:
- `seqnames`: Chromosome name
- `start`: Segment start position
- `end`: Segment end position
- `cnv`: Copy number variant type (amp/del)
- `cluster`: Assigned cluster

Example:
```
seqnames	start	end	cnv	cluster
chr1	0	1042457	del	2
chr1	1042457	1265484	del	2
chr1	9392992	9634910	amp	2
```

## Benchmark Evaluation

### Test Case
The `test_case/` directory contains a complete example with:
- Sample input data
- CONGASp training results
- CNV segment extraction
- Detailed analysis notebooks

### Performance Metrics
The benchmark evaluates CONGASp performance across:
- **Multiple patients**: 8 different patient samples
- **Data quality**: Various RNA-seq and ATAC-seq data characteristics
- **Computational efficiency**: Training time and memory usage
- **CNV detection accuracy**: Comparison with ground truth (if available)

## Dependencies

### Python Dependencies
- `torch`: PyTorch for deep learning
- `pyro`: Probabilistic programming
- `pandas`: Data manipulation
- `numpy`: Numerical computing
- `congas`: CONGASp package

### R Dependencies
- `Matrix`: Sparse matrix operations
- `readr`: Fast data reading
- `optparse`: Command-line argument parsing

### System Requirements
- Sufficient memory for large matrices (220MB+ per patient)
- Python 3.x and R 4.x
- Access to Google Cloud Storage (for data download)

## Troubleshooting

### Common Issues
1. **Memory errors**: Large matrices may require increased memory allocation
2. **Separator mismatch**: Ensure correct separator for each patient
3. **File permissions**: Check write permissions for output directories

### Error Handling
The batch scripts include error handling and continue processing other patients if one fails.

## References

For more information about CONGASp, see the original publication and documentation in the `CONGASp_github/` directory. 