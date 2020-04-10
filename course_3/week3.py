import numpy as np

from week2 import BLOSUM62

def align_with_affine_gap_penalty(v, w, σ, ε, score_table=BLOSUM62):
    n, m = len(v), len(w)

    lower = np.zeros((n+1, m+1), dtype=float) # Float otherwise np.inf doesn't work
    lower[0,:] = -np.inf
    middle = np.zeros((n+1, m+1), dtype=float)
    upper = np.zeros((n+1, m+1), dtype=float)
    upper[:,0] = -np.inf

    backtrack_l = np.empty((n+1, m+1), dtype=str)
    backtrack_m = np.empty((n+1, m+1), dtype=str)
    backtrack_u = np.empty((n+1, m+1), dtype=str)

    for i in range(n+1):
        for j in range(m+1):
            # Lower
            if i>0:
                lower[i, j] = max(lower[i-1, j]-ε, middle[i-1, j]-σ)
                if lower[i, j] == lower[i-1, j]-ε:
                    backtrack_l[i, j] = 'N' # North
                else:
                    backtrack_l[i, j] = 'M' # Middle
            # Upper
            if j>0:
                upper[i, j] = max(upper[i, j-1]-ε, middle[i, j-1]-σ)
                if upper[i, j] == upper[i, j-1]-ε:
                    backtrack_u[i, j] = 'W'
                else:
                    backtrack_u[i, j] = 'M'
            # Middle
            if (i,j) != (0, 0):
                if i > 0 and j == 0:
                    middle[i, j] = lower[i, j]
                    backtrack_m[i, j] = 'L'
                elif i==0 and j>0:
                    middle[i, j] = upper[i, j]
                    backtrack_m[i, j] = 'U'
                elif i>0 and j>0:
                    middle[i, j] = max(lower[i,j], upper[i,j],
                                    middle[i-1,j-1] + score_table.loc[v[i-1], w[j-1]])
                    if middle[i, j] == lower[i, j]:
                        backtrack_m[i, j] = 'L'
                    elif middle[i, j] == upper[i, j]:
                        backtrack_m[i, j] = 'U'
                    else:
                        backtrack_m[i, j] = 'NW'

    final_score = int(middle[n, m])
    curr = backtrack_m
    i, j = n, m
    rev_w = []
    rev_v = []
    while (i,j) != (0, 0):
        direction = curr[i, j]
        if direction == 'N':
            rev_v.append(v[i-1])
            rev_w.append('-')
            i -= 1
        elif direction == 'W':
            rev_v.append('-')
            rev_w.append(w[j-1])
            j -= 1
        elif direction == 'NW':
            rev_v.append(v[i-1])
            rev_w.append(w[j-1])
            i -= 1
            j -= 1
        elif direction == 'U':
            curr = backtrack_u
        elif direction == 'L':
            curr = backtrack_l
        elif direction == 'M': # Switch from lower/upper to middle grid
            if curr is backtrack_l:
                rev_v.append(v[i-1])
                rev_w.append('-')
                i -= 1
            elif curr is backtrack_u:
                rev_v.append('-')
                rev_w.append(w[j-1])
                j -= 1
            curr = backtrack_m
        else:
            raise ValueError('Unexpected if branch')

    align_v = ''.join(rev_v[::-1])
    align_w = ''.join(rev_w[::-1])
    return final_score, align_v, align_w


