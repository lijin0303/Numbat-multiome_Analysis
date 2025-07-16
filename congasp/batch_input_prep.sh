#!/bin/bash
# Usage: bash congasp/batch_input_prep.sh
# This script processes all .rds files in gs://numbat_wrapper/combined_rawcnt/:
# 1. Downloads each .rds file to a temp directory
# 2. Runs congasp/input_prep.R on it
# 3. Removes the .rds file after processing

set -euo pipefail

TMPDIR="congasp/tmp_rds"
mkdir -p "$TMPDIR"

# List all .rds files in the bucket
gsutil ls gs://numbat_wrapper/combined_rawcnt/*.rds | while read -r rds_gs_path; do
    rds_file=$(basename "$rds_gs_path")
    local_rds="$TMPDIR/$rds_file"
    echo "Processing $rds_gs_path ..."
    # Download
    gsutil cp "$rds_gs_path" "$local_rds"
    # Run input_prep.R
    Rscript congasp/input_prep.R -i "$local_rds"
    # Remove local copy
    rm -f "$local_rds"
    echo "Done with $rds_file"
done

echo "All files processed." 