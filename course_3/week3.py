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
                        backtrack_m[i, j] = 'D' # Diagonal
    '''
    print('lower')
    print(backtrack_l)
    print(lower)

    print('middle')
    print(backtrack_m)
    print(middle)

    print('upper')
    print(backtrack_u)
    print(upper)
    '''

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
        elif direction == 'D':
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


def calculate_next_col(prev_col, j, v, w, indel_penalty, score_table):
    '''
    j: index of column to be calculated
    '''
    col = np.zeros_like(prev_col, dtype=int)
    for i in range(len(col)):
            if i == 0:
                col[i] = prev_col[i] - indel_penalty
            else:
                col[i] = max(col[i-1] - indel_penalty, prev_col[i]-indel_penalty,
                            prev_col[i-1] + score_table.loc[v[i-1], w[j-1]])
    return col


def get_col(v, w, middle, indel_penalty=5, score_table=BLOSUM62):
    '''
    middle : index of column whose scores are to be returned
    '''
    n = len(v)
    j = 0
    col = np.arange(0, (n+1))*(-indel_penalty)
    while j < middle:
        j += 1
        prev_col = col
        col = calculate_next_col(prev_col, j, v, w, indel_penalty, score_table)
    return col


def get_middle_node(v, w, indel_penalty=5, score_table=BLOSUM62):
    m = len(w)
    middle = m//2
    from_source = get_col(v, w, middle, indel_penalty, score_table)
    to_sink_rev = get_col(v[::-1], w[::-1], m-middle, indel_penalty, score_table)
    to_sink = to_sink_rev[::-1]
    middle_col = from_source + to_sink
    i_max = np.argmax(middle_col)
    middle_node = (i_max, middle)
    return middle_node


def get_middle_edge(v, w, indel_penalty=5, score_table=BLOSUM62):
    m = len(w)
    middle = m//2
    from_source = get_col(v, w, middle, indel_penalty, score_table)
    from_source_after_middle = calculate_next_col(from_source, middle+1, v, w, indel_penalty, score_table)

    to_sink_before_middle_rev = get_col(v[::-1], w[::-1], m-middle-1, indel_penalty, score_table)
    to_sink_rev = calculate_next_col(to_sink_before_middle_rev, m-middle, v[::-1], w[::-1], indel_penalty, score_table)

    to_sink = to_sink_rev[::-1]
    middle_col = from_source + to_sink
    mid_max = np.argmax(middle_col)
    start = (mid_max, middle)

    # final_score = middle_col[mid_max]
    # print(f'final_score : {final_score}')

    after_middle_col = from_source_after_middle + to_sink_before_middle_rev[::-1]
    possible_end_nodes = (  (mid_max, middle+1),    # Horizontal
                            (mid_max+1, middle+1),  # Diagonal
                            (mid_max+1, middle))    # Vertical
    values = np.array([after_middle_col[mid_max], after_middle_col[mid_max+1], middle_col[mid_max+1]])
    end_node_idx = np.argmax(values)
    end = possible_end_nodes[end_node_idx]
    return (start, end)

def decode_path(v, w, path):
    v_align = []
    w_align = []
    i, j = 0, 0
    for direction in path:
        if direction == 'H':
            w_align.append(w[j])
            v_align.append('-')
            j+= 1
        elif direction == 'V':
            v_align.append(v[i])
            w_align.append('-')
            i+= 1
        elif direction == 'D':
            v_align.append(v[i])
            w_align.append(w[j])
            i+= 1
            j+= 1
        else:
            raise ValueError('Uknown direction')
    v_final = ''.join(v_align)
    w_final = ''.join(w_align)
    return v_final, w_final


def linear_space_align(v, w, indel_penalty=5, score_table=BLOSUM62):
    n, m = len(v), len(w)
    if (n, m) == (0, 0):
        return []
    elif n == 0 and m > 0:
        path = ['H']*m # Horizontal
        return path
    elif m == 0 and n > 0:
        path = ['V']*n # Vertical
        return path
    middle = m//2
    middle_edge = get_middle_edge(v, w, indel_penalty, score_table)
    middle_node = middle_edge[0]
    tail_node = middle_edge[1]

    i = middle_node[0]
    path = linear_space_align(v[:i], w[:middle], indel_penalty, score_table)

    if (middle_node[0] == tail_node[0]) and (middle_node[1]==tail_node[1]-1):
        edge_direction = 'H' # Horizontal
    elif (middle_node[0] == tail_node[0]-1) and (middle_node[1]==tail_node[1]):
        edge_direction = 'V' # Vertical
    elif (middle_node[0] == tail_node[0]-1) and (middle_node[1]==tail_node[1]-1):
        edge_direction = 'D' # Diagonal
    else:
        raise ValueError(f'Wrong direction for middle edge {middle_edge}')
    path.append(edge_direction)

    lower_i = tail_node[0]
    tail_j = tail_node[1]
    tail_path = linear_space_align(v[lower_i:], w[tail_j:], indel_penalty, score_table)
    path.extend(tail_path)

    return path
