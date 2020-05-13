from math import inf

from week4 import DirectedGraph

import sys
sys.path.append('..')
from course_2.week3 import AMINO_MASS, MASS_2_AMINO

import logging
logging.basicConfig()
logger = logging.getLogger('week5')
logger.setLevel(logging.ERROR)

def peptide_identification(spectral_vec, proteome, amino_mass=AMINO_MASS):
    graph = DirectedGraph()
    graph.add_node(0)
    graph.add_weight(0, 0)
    for idx, s_i in enumerate(spectral_vec):
        i = idx+1
        graph.add_node(i)
        graph.add_weight(i, s_i)
    m = i
    logger.debug(f'm = {m}')

    best_score = -inf
    best_path = []
    for start_idx, _ in enumerate(proteome):
        curr_node = 0
        score = 0
        idx = start_idx
        path = []
        nodes = [] # Debug

        logger.debug(f'Starting from {start_idx}')
        while True:
            nodes.append(curr_node) # Debug

            amino = proteome[idx]
            path.append(amino)
            score += graph.weights[curr_node]
            mass = amino_mass[amino]

            if curr_node == m:
                mypath = ''.join(path)
                logger.debug(f'{start_idx}, {score}, {mypath}')
                if score > best_score:
                    best_score = score
                    best_path = path[:-1] # One extra amino acid is added while processign last node
                    break

            curr_node = curr_node + mass
            idx += 1
            if curr_node > m or idx >= len(proteome):
                break
        logger.debug(f'{nodes}\n')

    final_peptide = ''.join(best_path)
    return final_peptide, best_score


def psm_search(spectral_vectors, proteome, threshold, amino_mass=AMINO_MASS):
    psm = set()
    for vec in spectral_vectors:
        peptide, score = peptide_identification(vec, proteome, amino_mass)
        if score >= threshold:
            psm.add(peptide)
    return psm


def size_of_spectral_dictionary(spectral_vector, threshold, max_score, amino_mass=AMINO_MASS, prob_factor=1):
    cache = {}
    def size(i, t, vec):
        if i < 0 or t < 0:
            result = 0
        elif i == 0:
            if t==0:
                result = 1
            else:
                result = 0
        else:
            if (i, t) in cache:
                return cache[(i, t)]
            result = prob_factor * sum( size(i-mass, t - vec[-1], vec[:i-mass])
                        for mass in amino_mass.values() )
            cache[(i, t)] = result
        #print(f'Result for {i}, {t} = {result}')
        return result

    m = len(spectral_vector)
    return sum(size(m, t, spectral_vector) for t in range(threshold, max_score+1))
