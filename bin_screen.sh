#!/usr/bin/env bash
set -euo pipefail

# where your .rds files live
DIR="benchmark"

# path to the script you want to call
NUMBAT_SCRIPT="${DIR}/numbat_chr7.sh"

# sanity check
if [[ ! -x "$NUMBAT_SCRIPT" ]]; then
  echo "ERROR: $NUMBAT_SCRIPT not found or not executable" >&2
  exit 1
fi

# loop all files like binxx_chr.rds
for f in "${DIR}"/bin*_chr7.rds; do
    # skip if no matches
    [[ -e "$f" ]] || { echo "No .rds files found."; break; }

    # extract the basename without .rds and _chr7
    base=$(basename "$f" .rds)
    base_no_chr7="${base/_chr7/}"

    echo "→ Processing '${base_no_chr7}' via numbat_chr7.sh"
    # call the other script, passing the pattern as its first argument
    "$NUMBAT_SCRIPT" "$base_no_chr7"
done
