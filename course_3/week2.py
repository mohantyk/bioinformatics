import string

import pandas as pd
import numpy as np

BLOSUM62 = pd.read_csv('BLOSUM62.csv')
PAM250 = pd.read_csv('PAM250.csv')

def global_alignment(v, w, score_table, indel_penalty=5, local_match=False):
    '''
    score_table : pandas table
    '''
    n = len(v)
    m = len(w)

    longest_path = np.zeros((n+1, m+1), dtype=int)
    longest_path[0, :] = np.arange(0, (m+1))*(-indel_penalty)
    longest_path[:, 0] = np.arange(0, (n+1))*(-indel_penalty)

    backtrack = np.empty((n+1, m+1), dtype=str)
    backtrack[0, :] = 'L'
    backtrack[:, 0] = 'U'

    for i in range(1, n+1):
        for j in range(1, m+1):
            down_path = longest_path[i-1, j] - indel_penalty
            right_path = longest_path[i, j-1] - indel_penalty
            diag_path = longest_path[i-1, j-1] + score_table.loc[v[i-1], w[j-1]]

            if local_match:
                longest_path[i, j] = max(0, down_path, right_path, diag_path)
            else:
                longest_path[i, j] = max(down_path, right_path, diag_path)

            if longest_path[i, j] == down_path:
                backtrack[i, j] = 'U'
            elif longest_path[i, j] == right_path:
                backtrack[i, j] = 'L'
            elif longest_path[i, j] == diag_path:
                backtrack[i, j] = 'D'
            elif local_match and longest_path[i, j] == 0:
                backtrack[i, j] = 'S' # Source

    rev_v = []
    rev_w = []
    if local_match:
        final_score = np.amax(longest_path)
        max_index = np.argmax(longest_path)
        (i, j) = np.unravel_index(max_index, longest_path.shape)
    else:
        final_score = longest_path[n, m]
        i, j = n, m

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
        elif direction == 'D':
            rev_v.append(v[i-1])
            rev_w.append(w[j-1])
            i -= 1
            j -= 1
        elif direction == 'S':
            break

    align_v = ''.join(rev_v[::-1])
    align_w = ''.join(rev_w[::-1])
    return final_score, align_v, align_w


def edit_distance(v, w):
    data = -np.ones((26, 26), int)
    np.fill_diagonal(data, 0)
    chars = list(string.ascii_uppercase)
    metrics = pd.DataFrame(data, index=chars, columns=chars)

    n, _, _ = global_alignment(v, w, metrics, 1)
    return -n
