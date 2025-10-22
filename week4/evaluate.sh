#!/bin/bash
# Run all alignment scripts with both Python and Codon, measure runtime, and format results.

# Define test cases
declare -A tests
tests=(
  ["code/global_alignment.py,test/MT-human.fa,test/MT-orang.fa"]="global-mt_human"
  ["code/global_alignment.py,test/q1.fa,test/t1.fa"]="global-q1"
  ["code/local_alignment.py,test/q2.fa,test/t2.fa"]="local-q2"
  ["code/affine_alignment.py,test/q3.fa,test/t3.fa"]="affine-q3"
  ["code/fitting_alignment.py,test/q4.fa,test/t4.fa"]="fitting-q4"
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
  codon run "$script" "$qfile" "$tfile" >/dev/null 2>&1
  end=$(date +%s%3N)
  runtime_cod=$((end - start))
  printf "%-18s %-10s %-10s\n" "$label" "codon" "${runtime_cod}ms"
done
