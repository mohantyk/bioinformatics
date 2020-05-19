from collections import defaultdict

import pandas as pd
import numpy as np

def suffix_array(text):
    idx_to_suffix = {idx: text[idx:] for idx, _ in enumerate(text)}
    arr = list(sorted(idx_to_suffix, key=idx_to_suffix.get))
    return arr

def burrows_wheeler_transform(text):
    cyclic_shifts = [text[shift:] + text[:shift] for shift in range(len(text))]
    lexic_shifts = sorted(cyclic_shifts)
    bwt = ''.join([pattern[-1] for pattern in lexic_shifts])
    return bwt

def invert_bwt(bwt):
    last_col = bwt
    first_col = ''.join(sorted(bwt))

    first_row = '$'
    last_ltr_count = 0
    for row_idx in range(1, len(bwt)-1):
        last_ltr = first_row[-1]

        # Find col idx of #last_ltr_count occurence of last_ltr in last_col
        count = 0
        for col_idx, ch in enumerate(last_col):
            if ch == last_ltr:
                if count == last_ltr_count:
                    break
                else:
                    count += 1
        nxt_ltr = first_col[col_idx]
        first_row = first_row + nxt_ltr
        #Find count of nxt_ltr
        last_ltr_count = len([ch for ch in first_col[:col_idx] if ch==nxt_ltr])

    first_row = first_row + last_col[0] # Add the last character of first row
    text = first_row[1:] + first_row[0] # Move the $ to the end
    return text


def bw_match(bw_text, pattern):
    '''
    Returns num of places the pattern matches the text encoded by bw_text
    '''
    first_col = ''.join(sorted(bw_text))
    last_col = bw_text

    first_col_indices = defaultdict(list) # Indices of each character in first_col
    for idx, ch in enumerate(first_col):
        first_col_indices[ch].append(idx)

    last_to_first = {} # Mapping from last col index to first col index
    last_col_counts = defaultdict(int)
    for idx, ch in enumerate(last_col):
        count = last_col_counts[ch]
        first_col_idx = first_col_indices[ch][count]
        last_to_first[idx] = first_col_idx
        last_col_counts[ch] += 1

    lpattern = list(pattern)
    top = 0
    bottom = len(first_col) - 1
    while top <= bottom:
        if lpattern:
            ltr = lpattern.pop()
            match_against = last_col[top:bottom+1]
            if ltr in match_against:
                first_idx = top+match_against.index(ltr)
                top = last_to_first[first_idx]

                reversed_idx = match_against[::-1].index(ltr)
                last_idx = bottom - reversed_idx
                bottom = last_to_first[last_idx]
            else:
                return 0

        else:
            return bottom - top + 1


def match_patterns(bw_text, patterns):
    num_matches = [bw_match(bw_text, pattern) for pattern in patterns]
    return num_matches

def viterbi(emitted, alphabet, states, transitions, emissions):
    transition_prob = pd.DataFrame( data=transitions, columns=list(states), index=list(states) )
    emission_prob = pd.DataFrame( data=emissions, columns=list(alphabet), index=list(states))

    graph = np.zeros((len(states), len(emitted)), dtype=float)
    backtrack = np.empty((len(states), len(emitted)), dtype=str)
    for idx, ltr in enumerate(emitted):
        if idx == 0:
            graph[:, idx] = np.log(emission_prob.loc[:, ltr]*1/len(states))
            continue
        for s_idx, state in enumerate(states):
            all_links = graph[:,idx-1] + np.log(transition_prob.loc[:,state].values) + np.log(emission_prob.loc[state, ltr])
            graph[s_idx, idx] = np.max(all_links)
            backtrack[s_idx, idx] = states[np.argmax(all_links)]

    path = []
    s_idx = np.argmax(graph[:, -1])
    state = states[s_idx]
    path.append(state)
    for idx in reversed(range(len(emitted))):
        s_idx = states.index(state)
        state = backtrack[s_idx, idx]
        path.append(state)

    final = ''.join(path[::-1])
    return final


def outcome_likelihood(emitted, alphabet, states, transitions, emissions):
    transition_prob = pd.DataFrame( data=transitions, columns=list(states), index=list(states) )
    emission_prob = pd.DataFrame( data=emissions, columns=list(alphabet), index=list(states))

    graph = np.zeros((len(states), len(emitted)), dtype=float)
    for idx, ltr in enumerate(emitted):
        if idx == 0:
            graph[:, idx] = emission_prob.loc[:, ltr]*(1/len(states))
            continue
        for s_idx, state in enumerate(states):
            graph[s_idx, idx] = np.dot(graph[:,idx-1], transition_prob.loc[:,state].values) * emission_prob.loc[state, ltr]
    final = np.sum(graph[:,-1])
    return final
