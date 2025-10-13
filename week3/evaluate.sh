#!/bin/bash
set -euxo pipefail

PATH=${PATH}:${HOME}/.codon/bin
echo "week3"

echo -e "\nRunning tests\n"
python3 code/compare_runtimes.py