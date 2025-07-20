#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# detect GNU date vs BSD date
if command -v gdate &>/dev/null; then
  DATE_CMD=gdate
  PARSE_OPT=(-d)
else
  DATE_CMD=date
  PARSE_OPT=(-j -f "%Y-%m-%d %H:%M:%S")
fi

# temp file to collect unsorted data
tmp=$(mktemp)
# write header to final TSV
out=runtime.tsv
printf "size_kb\tminutes\n" > "$out"

for log in patA_bin*_outputs/log.txt; do
  # pull numeric size (e.g. “20” from “patA_bin20kb_outputs”)
  folder=$(basename "$(dirname "$log")")
  size_kb=$(printf '%s' "$folder" | sed -E 's/^patA_bin([0-9]+)kb_outputs$/\1/')

  # grab timestamps
  start_line=$(grep 'INFO \[' "$log" | grep -v 'All done!' | head -n1 || true)
  end_line=$(grep -m1 'All done!' "$log" || true)
  [[ -z $start_line || -z $end_line ]] && continue

  start_ts=$(printf '%s\n' "$start_line" | sed -E 's/.*INFO \[([^]]+)\].*/\1/')
  end_ts=$(printf '%s\n' "$end_line"  | sed -E 's/.*INFO \[([^]]+)\].*/\1/')

  # parse into epoch seconds
  if ! start_s=$($DATE_CMD "${PARSE_OPT[@]}" "$start_ts" +%s 2>/dev/null); then continue; fi
  if ! end_s=$($DATE_CMD "${PARSE_OPT[@]}" "$end_ts"  +%s 2>/dev/null); then continue; fi

  # compute whole minutes
  delta=$(( end_s - start_s ))
  (( delta < 0 )) && continue
  minutes=$(( delta / 60 ))

  # collect
  printf "%s\t%s\n" "$size_kb" "$minutes" >> "$tmp"
done

# sort by first column (numeric) and append to out
sort -k1,1n "$tmp" >> "$out"
rm "$tmp"

echo "Wrote sorted runtimes to $out"