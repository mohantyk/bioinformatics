#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 29 13:05:42 2020

@author: kaniska
"""

from itertools import product
from math import log2, inf
from week2 import approx_pattern_match, neighbors, hamming


def motif_enumeration(dna, k , d):
    """
    Brute force seach for motifs

    Parameters
    ----------
    dna : list of dna sequences
    k : kmer length
    d : allowed mutations

    Returns
    -------
    motifs : set of motifs

    """
    n = len(dna[0])
    motifs = set()
    for idx in range(n-k):
        kmer = dna[0][idx:idx+k]
        for mutated in neighbors(kmer, d):
            for gene in dna[1:]:
                exists = approx_pattern_match(mutated, gene, d)
                if not exists:
                    break
            else:
                motifs.add(mutated)
                
    return motifs



def entropy(distrib):
    total = 0
    assert( sum(distrib) == 1.0 )
    for prob in distrib:
        try:
            total += prob*log2(prob)
        except ValueError:
            pass
    
    return -total


def distance_from_dna(pattern, dna):
    k = len(pattern)
    total_distance = 0
    for strand in dna:
        min_distance = inf
        for idx in range(len(strand)-k):
            kmer = strand[idx:idx+k]
            distance = hamming(kmer, pattern)
            min_distance = min(distance, min_distance)
        total_distance += min_distance
    return total_distance


def median_string(dna, k):
    min_distance = inf
    result = ''
    for perm in product(['A', 'C', 'T', 'G'], repeat=k):
        pattern = ''.join(perm)
        distance = distance_from_dna(pattern, dna)
        #print(pattern)
        if distance < min_distance:
            result = pattern
            min_distance = distance
    return result
    

def most_probable_kmer(dna, k, probs):
    """
    Finds the most proabable kmer in a dna 

    Parameters
    ----------
    dna : dna string
    k : kmer length
    probs : probability matrix, as a dictionary in hte form
                probs['A'] = [0.3, 0.2, 0.1]

    Returns
    -------
    most probable kmer string

    """
    n = len(dna)
    max_prob = 0
    for idx in range(n-k):
        kmer = dna[idx:idx+k]
        curr = 1.0
        for i, n in enumerate(kmer):
            curr *= probs[n][i]
        if curr > max_prob:
            result = kmer
            max_prob = curr
    return result
