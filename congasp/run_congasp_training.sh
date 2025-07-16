#!/bin/bash
# Usage: ./run_congasp_training.sh <separator_string> <output_dir> [input_dir]
# Example: ./run_congasp_training.sh 146p ./output ./

SEPARATOR=$1
OUTPUT_DIR=$2
INPUT_DIR=${3:-.}

python3 congasp/run_congasp_training.py --separator "$SEPARATOR" --output_dir "$OUTPUT_DIR" --input_dir "$INPUT_DIR" 