#!/bin/bash

# This script iterates through all subfolders in /workspace/inputData/ and runs run_congasp_training.py
# For 'patA' subfolder, use separator 'NIH-A'; for all others, use 'GEX'

INPUT_PARENT="/workspace/inputData"

for SUBDIR in "$INPUT_PARENT"/*/; do
    # Remove trailing slash and get the folder name
    SUBDIR_NAME=$(basename "${SUBDIR%/}")
    # Set separator
    if [ "$SUBDIR_NAME" = "patA" ]; then
        SEPARATOR="NIH-A"
    else
        SEPARATOR="GEX"
    fi
    echo "Processing $SUBDIR_NAME with separator $SEPARATOR"
    python3 /workspace/run_congasp_training.py \
        --separator "$SEPARATOR" \
        --output_dir "$SUBDIR" \
        --input_dir "$SUBDIR"
done

echo "All subfolders processed." 