import sys

def read_file(filename):
    sequences = []
    current_id = None
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_id = line[1:].strip()
                sequences.append((current_id, ""))  # tuple: (id, sequence)
            else:
                seq_id, seq = sequences[-1]
                sequences[-1] = (seq_id, seq + line.strip())
    return sequences

def fitting_alignment(seq1, seq2, match=3, mismatch=-3, gap=-2):
    seq1, seq2 = seq1.upper(), seq2.upper()
    m, n = len(seq1), len(seq2)

    # initialize score and pointer matrices
    matrix = [[0]*(n+1) for _ in range(m+1)]
    pointers = [[""]*(n+1) for _ in range(m+1)]

    # fill matrix
    for i in range(1, m+1):
        for j in range(1, n+1):
            score = match if seq1[i-1] == seq2[j-1] else mismatch
            diag = matrix[i-1][j-1] + score
            up = matrix[i-1][j] + gap
            left = matrix[i][j-1] + gap

            matrix[i][j] = max(diag, up, left)

            if matrix[i][j] == diag:
                pointers[i][j] = "diag"
            elif matrix[i][j] == up:
                pointers[i][j] = "up"
            else:
                pointers[i][j] = "left"

    # find the best score in last row
    max_score = max(matrix[m])
    j = matrix[m].index(max_score)

    # backtrack
    align1, align2 = "", ""
    i = m
    while i > 0 and j > 0:
        if pointers[i][j] == "diag":
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif pointers[i][j] == "up":
            align1 = seq1[i-1] + align1
            align2 = "-" + align2
            i -= 1
        elif pointers[i][j] == "left":
            align1 = "-" + align1
            align2 = seq2[j-1] + align2
            j -= 1
        else:
            break

    return matrix, (align1, align2, max_score)

def compare_pairs(query_file, target_file):
    queries = read_file(query_file)
    targets = read_file(target_file)

    # check for mismatch in pairs
    if len(queries) != len(targets):
        print("Number of queries and targets differ")
        print(f"{len(queries)} queries, {len(targets)} targets")
        return

    for (q_id, q_seq), (t_id, t_seq) in zip(queries, targets):
        print(f"\n=== {q_id} vs {t_id} ===")
        matrix, (a1, a2, max_score) = fitting_alignment(q_seq, t_seq)

        print("Score:", max_score)
        print(a1)
        print(a2)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Please provide filepaths to 2 sequences.")
        sys.exit(0)

    results = compare_pairs(sys.argv[1], sys.argv[2])