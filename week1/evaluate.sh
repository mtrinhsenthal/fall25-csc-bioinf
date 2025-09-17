#!/bin/bash

set -euxo pipefail

# Add Codon to PATH
export PATH="${HOME}/.codon/bin:$PATH"


calculate_n50() {
    local fasta_path="$1"
    python3 <<END
contig_lengths = []
with open("$fasta_path") as f:
    length = 0
    for line in f:
        if line.startswith('>'):
            if length > 0:
                contig_lengths.append(length)
                length = 0
        else:
            length += len(line.strip())
    if length > 0:
        contig_lengths.append(length)

contig_lengths.sort(reverse=True)
total_length = sum(contig_lengths)
half_length = total_length / 2

print(half_length)

running_sum = 0
for length in contig_lengths:
    running_sum += length
    if running_sum >= half_length:
        print(length)
        break

print("probably shouldn't be here...  ", running_sum, half_length)
END
}

format_runtime() {
    local raw="$1"
    local total_ms=$(printf "%.0f" "$(echo "$raw * 1000" | bc)")  # convert to milliseconds

    local minutes=$(( total_ms / 60000 ))
    local seconds=$(( (total_ms % 60000) / 1000 ))
    local millis=$(( total_ms % 1000 / 10 ))  # two-digit milliseconds

    printf "%d:%02d:%02d" "$minutes" "$seconds" "$millis"
}

python_script="week1/code/main.py"
codon_script="week1/code/codon_main.py"
args=("week1/data/data1" "week1/data/data2" "week1/data/data3" "week1/data/data4")
datafiles=()
language=()
runtimes=()
n50_results=()

for arg in "${args[@]}"; do
    language+=("python")
    datafiles+=("${arg#data/}")

    start_time=$(date +%s.%N)
    python3 "$python_script" "$arg"
    # python3 "$python_script" "$arg" > /dev/null # hide print statements
    end_time=$(date +%s.%N)

    runtime=$(echo "$end_time - $start_time" | bc)
    formatted_runtime=$(format_runtime "$runtime")
    runtimes+=("$formatted_runtime")

    n50=$(calculate_n50 "${arg}/contig.fasta")
    n50_results+=("$n50")


    language+=("codon")
    datafiles+=("${arg#data/}")
    start_time=$(date +%s.%N)
    codon run -release "$codon_script" "$arg" > /dev/null # hide print statements
    end_time=$(date +%s.%N)

    runtime=$(echo "$end_time - $start_time" | bc)
    formatted_runtime=$(format_runtime "$runtime")
    runtimes+=("$formatted_runtime")

    n50=$(calculate_n50 "${arg}/contig.fasta")
    n50_results+=("$n50")
done



# Print results in a table
printf "\n%-15s | %-15s | %-15s | %-10s\n" "Dataset" "Language" "Time (s)" "N50"
printf "%s\n" "-----------------------------------------------------------------"

# for i in "${!args[@]}"; do
#     args[$i]="${args[$i]#data/}" # strip "data/" from start of string
# done

for i in "${!language[@]}"; do
    printf "%-15s | %-15s | %-15s | %-10s\n" "${datafiles[$i]}" "${language[$i]}" "${runtimes[$i]}" "${n50_results[$i]}"
done
