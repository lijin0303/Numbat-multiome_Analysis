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
    
    # Run with error handling and memory management
    if python3 /workspace/run_congasp_training.py \
        --separator "$SEPARATOR" \
        --output_dir "$SUBDIR" \
        --input_dir "$SUBDIR"; then
        echo "Successfully processed $SUBDIR_NAME"
    else
        echo "Failed to process $SUBDIR_NAME - continuing with next file"
    fi
    
    # Sleep to allow memory cleanup
    sleep 5
done

echo "All subfolders processed." 