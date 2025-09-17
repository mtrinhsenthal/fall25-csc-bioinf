# python code/main.py data/data1
# python code/main.py data/data2
# python code/main.py data/data3

python_script="code/main.py"
codon_script="code/codon_main.py"
args=("data/data1" "data/data2" "data/data3")
language=()
runtimes=()
n50_results=()

for arg in "${args[@]}"; do
    echo "Executing $arg"
    echo "Running Python script..." 

    start_time=$(date +%s.%N)
    python "$python_script" "$arg" > /dev/null # hide print statements
    end_time=$(date +%s.%N)

    runtime=$(echo "$end_time - $start_time" | bc)
    language+=("python")
    runtimes+=("$runtime")

    echo "Calculating n50 for $arg with python"

    n50=$(python <<END
fasta_path = "${arg}/contig.fasta"
contig_lengths = []
with open(fasta_path) as f:
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

running_sum = 0
for length in contig_lengths:
    running_sum += length
    if running_sum >= half_length:
        print(length)
        break
END
)
    n50_results+=("$n50")


    echo "Running codon script..."
    start_time=$(date +%s.%N)
    codon run -release "$codon_script" "$arg" > /dev/null # hide print statements
    end_time=$(date +%s.%N)

    runtime=$(echo "$end_time - $start_time" | bc)
    language+=("codon")
    runtimes+=("$runtime")

    echo "Calculating n50 for $arg with codon"

    n50=$(python <<END
fasta_path = "${arg}/contig.fasta"
contig_lengths = []
with open(fasta_path) as f:
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

running_sum = 0
for length in contig_lengths:
    running_sum += length
    if running_sum >= half_length:
        print(length)
        break
END
)
    n50_results+=("$n50")
done



# Print results in a table
printf "\n%-15s | %-15s | %-15s | %-10s\n" "Dataset" "Language" "Time (s)" "N50"
printf "%s\n" "-----------------------------------------------------------------"

for i in "${!language[@]}"; do
    printf "%-15s | %-15s | %-15s | %-10s\n" "${args[$((i % 3))]}" "${language[$i]}" "${runtimes[$i]}" "${n50_results[$i]}"
done
