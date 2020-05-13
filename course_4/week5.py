from math import inf

from week4 import DirectedGraph

import sys
sys.path.append('..')
from course_2.week3 import AMINO_MASS, MASS_2_AMINO

def peptide_identification(spectral_vec, proteome, amino_mass=AMINO_MASS):
    graph = DirectedGraph()
    graph.add_node(0)
    graph.add_weight(0, 0)
    for idx, s_i in enumerate(spectral_vec):
        i = idx+1
        graph.add_node(i)
        graph.add_weight(i, s_i)
    m = i

    best_score = -inf
    best_path = []
    for start_idx, _ in enumerate(proteome):
        curr_node = 0
        score = 0
        idx = start_idx
        path = []

        while True:
            amino = proteome[idx]
            path.append(amino)
            score += graph.weights[curr_node]
            mass = amino_mass[amino]

            if curr_node == m:
                if score > best_score:
                    best_score = score
                    best_path = path[:]
                    break

            curr_node = curr_node + mass
            idx += 1
            if curr_node > m or idx >= len(proteome):
                break

    return ''.join(best_path)




