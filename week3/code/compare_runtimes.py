# import subprocess
# import time

# def run_command(cmd):
#     start = time.perf_counter()
#     result = subprocess.run(cmd, capture_output=True, text=True)
#     end = time.perf_counter()
#     runtime_ms = (end - start) * 1000
#     return runtime_ms, result.stdout.strip(), result.stderr.strip()

# # Commands for Python and Codon
# commands = {
#     "python": ["python3", "week3/code/python_versions/p_test_phylo.py"],
#     "codon": ["codon", "run", "-release", "week3/code/test_phylo.py"],
# }

# results = {}

# for lang, cmd in commands.items():
#     runtime, stdout, stderr = run_command(cmd)
#     results[lang] = runtime
#     if stdout:
#         print(stdout)
#     if stderr:
#         print(stderr)

# # Print a nice timing table
# print("\nLanguage    Runtime (ms)")
# print("-------------------------")
# for lang, runtime in results.items():
#     print(f"{lang:<10} {runtime:.2f}")

import subprocess

def get_runtime(cmd):
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    # Expect the file to print only a number (runtime in ms)
    try:
        return float(result.stdout.strip())
    except ValueError:
        print(f"Unexpected output from {cmd}:")
        print(result.stdout)
        print(result.stderr)
        return float('nan')

commands = {
    "python": ["python3", "week3/code/python_versions/p_test_phylo.py"],
    "codon": ["codon", "run", "-release", "week3/code/test_phylo.py"],
}

results = {lang: get_runtime(cmd) for lang, cmd in commands.items()}

print("\nLanguage    Runtime (ms)")
print("-------------------------")
for lang, runtime in results.items():
    print(f"{lang:<10} {runtime:.2f}")
