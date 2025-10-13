#!/bin/bash
# set -euxo pipefail

PATH=${PATH}:${HOME}/.codon/bin
echo "week3"

# Paths to test scripts
PYTHON_TEST="week3/code/python_versions/p_test_phylo.py"
CODON_TEST="week3/code/test_phylo.py"

# Run Python tests and capture runtime
PY_MS=$(python3 "$PYTHON_TEST"
PY_MS=$(printf "%.2f" "$PY_MS"))

# Run Codon tests and capture runtime
CODON_MS=$(codon run -release "$CODON_TEST")
CODON_MS=$(printf "%.2f" "$CODON_MS")

# Print the timing table
echo
echo "Language    Runtime (ms)"
echo "-------------------------"
printf "python     %s\n" "$PY_MS"
printf "codon      %s\n" "$CODON_MS"
