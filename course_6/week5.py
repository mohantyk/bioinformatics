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
    for pattern in alignment:
        path = ['S']
        curr = 0
        for idx, ltr in enumerate(pattern):
            if idx in star_alignment:
                curr += 1
                if ltr != '-':
                    path.append('M'+str(curr)) # Letter emitted
                else:
                    path.append('D'+str(curr)) # Deletion

            else:
                if ltr != '-':
                    path.append('I'+str(curr)) # Insertion
        path.append('E')
        paths.append(path)

    incoming = defaultdict(int)
    outgoing = defaultdict(dict)
    for path in paths:
        for idx, src in enumerate(path[:-1]):
            dst = path[idx+1]
            incoming[dst] += 1
            if dst in outgoing[src]:
                outgoing[src][dst] += 1
            else:
                outgoing[src][dst] = 1
    print(incoming)
    print(outgoing)

    transitions = pd.DataFrame(0, columns=nodes, index=nodes, dtype=float)
    for node in nodes:
        denom = incoming[node]
        if denom == 0:
            continue
        for out, edges in outgoing[node].items():
            transitions.loc[node, out] = edges/denom

    return transitions, None



if __name__ == '__main__':
    threshold = 0.289
    alphabet = 'ABCDE'
    multiple_alignment = ['EBA', 'E-D', 'EB-', 'EED', 'EBD', 'EBE', 'E-D', 'E-D']
    create_profile_hmm(multiple_alignment, threshold, multiple_alignment)