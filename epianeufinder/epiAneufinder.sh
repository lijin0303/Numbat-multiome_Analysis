#!/usr/bin/env bash
set -euo pipefail

# Path to your inputs file
INPUTS="epiAneufinder_inputs.txt"

# Path to your blacklist BED (adjust as needed)
BLACKLIST="hg38-blacklist.v2.bed"

# R script to run
RSCRIPT="epiAneufinder.R"

# Cell barcode directory
CB_DIR="tmp_CB"

# GCS bucket for uploads
GCS_BUCKET="gs://numbat_wrapper/epianeufinder/"

# Timing log file
TIMING_LOG="epiAneufinder_timing.log"

# Initialize timing log with header if it doesn't exist
if [[ ! -f "$TIMING_LOG" ]]; then
  echo -e "sample\tstart_time\tend_time\tduration_seconds\tduration_formatted" > "$TIMING_LOG"
fi

# Loop over each line: sampleid <tab> frag <tab> frag_tbi
while IFS=$'\t' read -r sampleid frag frag_tbi CB; do
  outdir="epianeufinder_${sampleid}"
  karyo="${outdir}/epiAneufinder_results/results_table.tsv"

  # If the karyogram already exists, skip this sample
  if [[ -f "$karyo" ]]; then
    echo ">>> Skipping $sampleid: '$karyo' already exists."
    continue
  fi

  echo "=== Processing sample: $sampleid ==="
  start_time=$(date +%s)
  start_time_readable=$(date)
  echo "  Start time: $start_time_readable"
  mkdir -p "$outdir"

  # Check if fragment files already exist, copy if needed
  frag_local="${outdir}/$(basename "$frag")"
  frag_tbi_local="${outdir}/$(basename "$frag_tbi")"
  
  if [[ -f "$frag_local" ]]; then
    echo "  Fragment file already exists: $frag_local"
  else
    echo "  Copying fragment file..."
    gsutil -m cp "$frag" "$outdir/"
  fi
  
  if [[ -f "$frag_tbi_local" ]]; then
    echo "  Fragment index already exists: $frag_tbi_local"
  else
    echo "  Copying fragment index..."
    gsutil -m cp "$frag_tbi" "$outdir/"
  fi

  # Check for cell barcode file
  cb_file="${CB_DIR}/${sampleid}_CB.txt"
  if [[ -f "$cb_file" ]]; then
    echo "  Found cell barcode file: $cb_file"
    selected_cells_arg="--selected_cells $cb_file"
  else
    echo "  No cell barcode file found for $sampleid, processing all cells"
    selected_cells_arg=""
  fi

  # Run the R script
  echo "  Running epiAneufinder for $sampleid..."
  Rscript "$RSCRIPT" \
    --frag "$frag_local" \
    --blacklist "$BLACKLIST" \
    --outdir "$outdir" \
    $selected_cells_arg

  echo "  → Output in $outdir/"
  
  # Upload CNV calls file to GCS bucket
  cnv_calls_file="${outdir}/epiAneufinder_results/${sampleid}_cnv_calls.tsv"
  if [[ -f "$cnv_calls_file" ]]; then
    echo "  Uploading CNV calls to GCS bucket..."
    gsutil -m cp "$cnv_calls_file" "$GCS_BUCKET"
    echo "  ✓ Uploaded: $cnv_calls_file"
  else
    echo "  ⚠ Warning: CNV calls file not found: $cnv_calls_file"
  fi
  
  # Calculate and display execution time
  end_time=$(date +%s)
  duration=$((end_time - start_time))
  hours=$((duration / 3600))
  minutes=$(((duration % 3600) / 60))
  seconds=$((duration % 60))
  
  echo "  End time: $(date)"
  echo "  Execution time: ${hours}h ${minutes}m ${seconds}s (${duration} seconds total)"
  
  # Log timing information to file
  echo -e "${sampleid}\t${start_time_readable}\t$(date)\t${duration}\t${hours}h ${minutes}m ${seconds}s" >> "$TIMING_LOG"
  
  echo
done < "$INPUTS"