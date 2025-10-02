from typing import List
from numpy import ndarray, float32

def calculate(sequence: str, matrix: ndarray[float,2], scores: ndarray[float32,1]):
    m = matrix.shape[0]
    n = scores.shape[0]

    for i in range(n):
        score = 0.0
        ok = 1

        for j in range(m):
            c = sequence[i+j]

            if c in 'Aa': score += matrix[j][0]
            elif c in 'Cc': score += matrix[j][1]
            elif c in 'Gg': score += matrix[j][2]
            elif c in 'Tt': score += matrix[j][3]
            else: ok = 0
        
        if ok: scores[i] = float32(score)
        else: scores[i] = float32('nan')