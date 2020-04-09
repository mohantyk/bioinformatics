import numpy as np

from week2 import BLOSUM62

def align_with_affine_gap_penalty(v, w, σ, ε, score_table=BLOSUM62):
    n, m = len(v), len(w)
    lower = np.zeros((n+1, m+1), dtype=int)
    middle = np.zeros((n+1, m+1), dtype=int)
    upper = np.zeros((n+1, m+1), dtype=int)

    for i in range(n+1):
        for j in range(m+1):
            if i>0:
                lower[i, j] = max(lower[i-1, j]-ε, middle[i-1, j]-σ)
            if j>0:
                upper[i, j] = max(upper[i, j-1]-ε, middle[i, j-1]-σ)
            if i>0 and j>0:
                middle[i, j] = max(lower[i,j], upper[i,j],
                                    middle[i-1,j-1] + score_table.loc[v[i-1], w[j-1]])

    return middle[n, m]


