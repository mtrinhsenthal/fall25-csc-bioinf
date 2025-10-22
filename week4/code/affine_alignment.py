import math
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

def affine_alignment(seq1, seq2, match=3, mismatch=-3, gap_open=-5, gap_extend=-1):
    seq1, seq2 = seq1.upper(), seq2.upper()
    m, n = len(seq1), len(seq2)
    NEG_INF = -math.inf

    # score matrices
    M = [[NEG_INF]*(n+1) for _ in range(m+1)]
    X = [[NEG_INF]*(n+1) for _ in range(m+1)]
    Y = [[NEG_INF]*(n+1) for _ in range(m+1)]

    # pointer matrices
    ptr_M = [["stop"]*(n+1) for _ in range(m+1)]
    ptr_X = [["stop"]*(n+1) for _ in range(m+1)]
    ptr_Y = [["stop"]*(n+1) for _ in range(m+1)]

    # initialization
    M[0][0] = 0
    for i in range(1, m+1):
        X[i][0] = gap_open + (i-1)*gap_extend
        M[i][0] = X[i][0]
        ptr_X[i][0] = "X"
    for j in range(1, n+1):
        Y[0][j] = gap_open + (j-1)*gap_extend
        M[0][j] = Y[0][j]
        ptr_Y[0][j] = "Y"

    # fill matrices
    for i in range(1, m+1):
        for j in range(1, n+1):
            score = match if seq1[i-1] == seq2[j-1] else mismatch

            # match matrix M
            candidates = [
                (M[i-1][j-1], "M"),
                (X[i-1][j-1], "X"),
                (Y[i-1][j-1], "Y")
            ]
            best_val, best_ptr = max(candidates, key=lambda x: x[0])
            M[i][j] = best_val + score
            ptr_M[i][j] = best_ptr

            # gap in seq2 (X)
            candidates = [
                (M[i-1][j] + gap_open + gap_extend, "M"),  # open new gap
                (X[i-1][j] + gap_extend, "X")              # extend
            ]
            best_val, best_ptr = max(candidates, key=lambda x: x[0])
            X[i][j] = best_val
            ptr_X[i][j] = best_ptr

            # gap in seq1 (Y)
            candidates = [
                (M[i][j-1] + gap_open + gap_extend, "M"),
                (Y[i][j-1] + gap_extend, "Y")
            ]
            best_val, best_ptr = max(candidates, key=lambda x: x[0])
            Y[i][j] = best_val
            ptr_Y[i][j] = best_ptr

    # find best final score
    final_scores = [(M[m][n], "M"), (X[m][n], "X"), (Y[m][n], "Y")]
    score, matrix = max(final_scores, key=lambda x: x[0])

    # backtrack
    aligned1, aligned2 = [], []
    i, j = m, n
    curr_matrix = matrix

    while i > 0 or j > 0:
        if curr_matrix == "stop":
            break

        if curr_matrix == "M":
            prev = ptr_M[i][j]
            aligned1.append(seq1[i-1])
            aligned2.append(seq2[j-1])
            i, j = i-1, j-1
            curr_matrix = prev

        elif curr_matrix == "X":
            prev = ptr_X[i][j]
            aligned1.append(seq1[i-1])
            aligned2.append("-")
            i -= 1
            curr_matrix = prev

        elif curr_matrix == "Y":
            prev = ptr_Y[i][j]
            aligned1.append("-")
            aligned2.append(seq2[j-1])
            j -= 1
            curr_matrix = prev

        else:
            break

    # reverse to get final alignment
    aligned1 = "".join(reversed(aligned1))
    aligned2 = "".join(reversed(aligned2))

    return score, aligned1, aligned2

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
        score, a1, a2 = affine_alignment(q_seq, t_seq)

        print("Score:", score)
        print(a1)
        print(a2)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Please provide filepaths to 2 sequences.")
        sys.exit(0)

    results = compare_pairs(sys.argv[1], sys.argv[2])