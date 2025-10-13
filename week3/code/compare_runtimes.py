# import subprocess
# import time

# results = {}

# # Run Python tests
# start = time.perf_counter()
# subprocess.run(["python3", "python_versions/p_test_phylo.py"], check=True)
# end = time.perf_counter()
# results["python"] = (end - start) * 1000

# # Run Codon tests
# start = time.perf_counter()
# subprocess.run(["codon", "run", "test_phylo.py"], check=True)
# end = time.perf_counter()
# results["codon"] = (end - start) * 1000

# # Print table
# print("Language    Runtime")
# print("-------------------")
# for lang, ms in results.items():
#     print(f"{lang:<10} {ms:.0f}ms")

import subprocess
import time

def run_command(cmd):
    start = time.perf_counter()
    result = subprocess.run(cmd, capture_output=True, text=True)
    end = time.perf_counter()
    runtime_ms = (end - start) * 1000
    return runtime_ms, result.stdout.strip(), result.stderr.strip()

# Commands for Python and Codon
commands = {
    "python": ["python3", "week3/code/python_versions/p_test_phylo.py"],
    "codon": ["codon", "run", "-release", "week3/code/test_phylo.py"],
}

results = {}

for lang, cmd in commands.items():
    runtime, stdout, stderr = run_command(cmd)
    results[lang] = runtime
    if stdout:
        print(stdout)
    if stderr:
        print(stderr)

# Print a nice timing table
print("\nLanguage    Runtime (ms)")
print("-------------------------")
for lang, runtime in results.items():
    print(f"{lang:<10} {runtime:.2f}")
