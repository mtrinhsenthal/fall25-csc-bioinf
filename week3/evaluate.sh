#!/bin/bash
set -euxo pipefail

PATH=${PATH}:${HOME}/.codon/bin
echo "week3"

# echo -e "\nRunning tests\n"
# python3 week3/code/compare_runtimes.py

# Paths to test scripts
PYTHON_TEST="week3/code/python_versions/p_test_phylo.py"
CODON_TEST="week3/code/test_phylo.py"

# Run Python tests and capture runtime
echo "Running Python tests..."
PY_MS=$(python3 "$PYTHON_TEST")

# Run Codon tests and capture runtime
echo "Running Codon tests..."
CODON_MS=$(codon run -release "$CODON_TEST")

# Print the timing table
echo
echo "Language    Runtime (ms)"
echo "-------------------------"
printf "python     %s\n" "$PY_MS"
printf "codon      %s\n" "$CODON_MS"
