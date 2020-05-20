from collections import defaultdict
from itertools import chain
import pandas as pd

def create_profile_hmm(alignment, threshold, alphabet):
    star_alignment = set()
    num_patterns = len(alignment)
    for idx, _ in enumerate(alignment[0]):
        deletions = sum(1 for pattern in alignment
                            if pattern[idx] == '-')
        if deletions <= threshold * num_patterns:
            star_alignment.add(idx)

    nodes = ['S', 'I0'] + list(chain.from_iterable(('M'+str(count+1), 'D'+str(count+1), 'I'+str(count+1))
                                                    for count in range(len(star_alignment)))) + ['E']

    paths = []
    emissions = pd.DataFrame(0, columns=list(alphabet), index=nodes, dtype=float)
    for pattern in alignment:
        path = ['S']
        curr = 0
        for idx, ltr in enumerate(pattern):
            if idx in star_alignment:
                curr += 1
                if ltr != '-': # Letter emitted
                    node = 'M'+str(curr)
                    emissions.loc[node, ltr] += 1
                    path.append(node)
                else:
                    path.append('D'+str(curr)) # Deletion

            else:
                if ltr != '-': # Letter emitted
                    node = 'I'+str(curr)
                    emissions.loc[node, ltr] += 1
                    path.append(node)
        path.append('E')
        paths.append(path)
    # Calculate final emissions probabilities
    emissions = emissions.div(emissions.sum(axis=1), axis=0).fillna(0).round(3)

    incoming = defaultdict(int); incoming['S'] = len(paths)
    outgoing = defaultdict(dict)
    for path in paths:
        for idx, src in enumerate(path[:-1]):
            dst = path[idx+1]
            incoming[dst] += 1
            if dst in outgoing[src]:
                outgoing[src][dst] += 1
            else:
                outgoing[src][dst] = 1

    # Calculate final transition probabilities
    transitions = pd.DataFrame(0, columns=nodes, index=nodes, dtype=float)
    for node in nodes:
        denom = incoming[node]
        if denom == 0:
            continue
        for out, edges in outgoing[node].items():
            transitions.loc[node, out] = edges/denom

    return transitions, emissions
