
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

def global_alignment(seq1, seq2, match=3, mismatch=-3, gap=-2):
    seq1 = seq1.upper()
    seq2 = seq2.upper()

    m = len(seq1)
    n = len(seq2)

    # initialize matrix
    matrix = [[0] * (n+1) for _ in range(m+1)]
    pointers = [[""] * (n+1) for _ in range(m+1)] # matrix for storing pointers

    # initialize first row and column
    # since there is only one way to reach them (gaps)
    for i in range(1, m+1):
        matrix[i][0] = gap * i
        pointers[i][0] = "up"
    for j in range(1, n+1):
        matrix[0][j] = gap * j
        pointers[0][j] = "left"

    for i in range(1, m+1):
        for j in range(1, n+1):
            score = match if seq1[i-1] == seq2[j-1] else mismatch

            diag = matrix[i-1][j-1] + score
            up = matrix[i-1][j] + gap
            left = matrix[i][j-1] + gap

            best = max(diag, up, left)
            matrix[i][j] = best

            # track direction for backtracking
            if (best == diag):
                pointers[i][j] = "diag"
            elif (best == up):
                pointers[i][j] = "up"
            else:
                pointers[i][j] = "left"

    # backtrack
    align1, align2 = "", ""
    i, j = m, n
    while i > 0 or j > 0:
        if pointers[i][j] == "diag":
            align1 = seq1[i - 1] + align1
            align2 = seq2[j - 1] + align2
            i -= 1
            j -= 1
        elif pointers[i][j] == "up":
            align1 = seq1[i - 1] + align1
            align2 = "-" + align2
            i -= 1
        elif pointers[i][j] == "left":
            align1 = "-" + align1
            align2 = seq2[j - 1] + align2
            j -= 1
        else:
            break

    return matrix, (align1, align2, matrix[m][n])

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
        alignment_matrix, (a1, a2, score) = global_alignment(q_seq, t_seq)

        print("Score:", score)
        print(a1)
        print(a2)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Please provide filepaths to 2 sequences.")
        sys.exit(0)

    # with open(sys.argv[1], 'r') as f1:
    #     seq1 = ''.join(line.strip() for line in f1 if not line.startswith('>'))

    # with open(sys.argv[2], 'r') as f2:
    #     seq2 = ''.join(line.strip() for line in f2 if not line.startswith('>'))

    # matrix, (align1, align2, weight) = global_alignment(seq1, seq1)
    # print(weight)

    results = compare_pairs(sys.argv[1], sys.argv[2])