#!/bin/bash
# set -euxo pipefail
PATH=${PATH}:${HOME}/.codon/bin
echo "week2"

echo -e "\nRunning codon tests\n"
codon run -release week2/code/tests_motif.py

echo "----------------------------"
echo -e "\nRunning python tests\n"
python3 week2/code/tests_python.py