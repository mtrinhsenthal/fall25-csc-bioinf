#!/bin/bash
# set -euxo pipefail 
PATH=${PATH}:${HOME}/.codon/bin
echo "week4"
# Run all alignment scripts with both Python and Codon, measure runtime, and format results.

# Define test cases
declare -A tests
tests=(
  ["week4/code/global_alignment.py,week4/test/MT-human.fa,week4/test/MT-orang.fa"]="global-mt_human"
  ["week4/code/local_alignment.py,week4/test/MT-human.fa,week4/test/MT-orang.fa"]="local-mt_human"
  ["week4/code/affine_alignment.py,week4/test/MT-human.fa,week4/test/MT-orang.fa"]="affine-mt_human"
  ["week4/code/fitting_alignment.py,week4/test/MT-human.fa,week4/test/MT-orang.fa"]="fitting-mt_human"
  ["week4/code/global_alignment.py,week4/test/q1.fa,week4/test/t1.fa"]="global-q1"
  ["week4/code/local_alignment.py,week4/test/q1.fa,week4/test/t1.fa"]="local-q1"
  ["week4/code/affine_alignment.py,week4/test/q1.fa,week4/test/t1.fa"]="affine-q1"
  ["week4/code/fitting_alignment.py,week4/test/q1.fa,week4/test/t1.fa"]="fitting-q1"
)

# Output table header
printf "%-18s %-10s %-10s\n" "Method" "Language" "Runtime"
echo "-----------------------------------------------"

# Loop through test cases
for key in "${!tests[@]}"; do
  IFS=',' read -r script qfile tfile <<< "$key"
  label="${tests[$key]}"

  # --- Python run ---
  start=$(date +%s%3N)
  python3 "$script" "$qfile" "$tfile" >/dev/null 2>&1
  end=$(date +%s%3N)
  runtime_py=$((end - start))
  printf "%-18s %-10s %-10s\n" "$label" "python" "${runtime_py}ms"

  # --- Codon run ---
  start=$(date +%s%3N)
  codon run -release "$script" "$qfile" "$tfile" >/dev/null 2>&1
  end=$(date +%s%3N)
  runtime_cod=$((end - start))
  printf "%-18s %-10s %-10s\n" "$label" "codon" "${runtime_cod}ms"
done
