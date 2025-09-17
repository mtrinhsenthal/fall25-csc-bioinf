# python code/main.py data/data1
# python code/main.py data/data2
# python code/main.py data/data3

python_script="code/main.py"
n50_script = "n50.py"
args=("data/data1" "data/data2" "data/data3")
results=()
runtimes=()

echo "Running Python script..."

for arg in "${args[@]}"; do
    echo "Executing $arg"

    start_time=$(date +%s.%N)
    # output=$(python "$python_script" "$arg")
    output="hello"
    python "$python_script" "$arg" > /dev/null # hide print statements
    end_time=$(date +%s.%N)

    runtime=$(echo "$end_time - $start_time" | bc)
    results+=("$output")
    runtimes+=("$runtime")
done

# Print results in a table
printf "\n%-15s | %-30s | %-10s\n" "Data" "Result" "Time (s)"
printf "%s\n" "---------------------------------------------------------------"

for i in "${!args[@]}"; do
    printf "%-15s | %-30s | %-10s\n" "${args[$i]}" "${results[$i]}" "${runtimes[$i]}"
done