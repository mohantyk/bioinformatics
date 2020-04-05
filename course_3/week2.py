import pandas as pd
import numpy as np

def global_alignment(v, w, indel_penalty=5):
    n = len(v)
    m = len(w)

    score_table = pd.read_csv('BLOSUM62.csv')
    longest_path = np.zeros((n+1, m+1), dtype=int)
    backtrack = np.empty((n+1, m+1), dtype=str)
    backtrack[0, :] = 'L'
    backtrack[:, 0] = 'U'

    for i in range(1, n+1):
        for j in range(1, m+1):
            down_path = longest_path[i-1, j] + indel_penalty
            right_path = longest_path[i, j-1] + indel_penalty
            diag_path = longest_path[i-1, j-1] + score_table.loc[v[i-1], w[j-1]]

            longest_path[i, j] = max(down_path, right_path, diag_path)
            if longest_path[i, j] == down_path:
                backtrack[i, j] = 'U'
            elif longest_path[i, j] == right_path:
                backtrack[i, j] = 'L'
            else:
                backtrack[i, j] = 'D'

    rev_v = []
    rev_w = []
    i, j = len(v), len(w)

    while (i, j) != (0, 0):
        direction = backtrack[i, j]
        if direction == 'U':
            rev_v.append(v[i-1])
            rev_w.append('-')
            i -= 1
        elif direction == 'L':
            rev_v.append('-')
            rev_w.append(w[j-1])
            j -= 1
        else:
            rev_v.append(v[i-1])
            rev_w.append(w[j-1])
            i -= 1
            j -= 1

    align_v = ''.join(rev_v[::-1])
    align_w = ''.join(rev_w[::-1])
    return align_v, align_w
