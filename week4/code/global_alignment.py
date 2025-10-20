
# 1. form a 2d array using the recurrence relation for dynamic programming
# 2. create array containing "backtracking pointers"
# 3. after reaching sink, backtrack to source to produce a max-weight path
# 4. infer the alignment corresponding to this path

def global_alignment(seq1, seq2, match=3, mismatch=-3, gap=-2):
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
    
    # return matrix


# Example usage:
A = "ABCB"
B = "ABB"

matrix, (align1, align2, weight) = global_alignment(A, B)
print("\nScore Matrix:")
for row in matrix:
    print(row)

print(align1)
print(align2)
print(weight)