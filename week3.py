#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 29 13:05:42 2020

@author: kaniska
"""

from itertools import product
from math import log2, inf
from collections import defaultdict, Counter

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
        if distance < min_distance:
            result = pattern
            min_distance = distance
    return result



def profile_score(motifs):
    """
    Calculates profile score 

    Parameters
    ----------
    motifs : list of dna strings

    Returns
    -------
    score : number of mismatches with the best motif

    """
    n = len(motifs[0])
    correct = []
    for idx in range(n):
        count = Counter(dna[idx] for dna in motifs)
        correct.append(count.most_common(1)[0][0])
    motif = ''.join(correct)
    
    score = 0
    for dna in motifs:
        score += hamming(motif, dna)
    return score

    

def most_probable_kmer(dna, k, probs):
    """
    Finds the most proabable kmer in a dna 

    Parameters
    ----------
    dna : dna string
    k : kmer length
    probs : probability matrix, as a dictionary in the form
                probs['A'] = [0.3, 0.2, 0.1]

    Returns
    -------
    most probable kmer string

    """
    n = len(dna)
    max_prob = -1
    for idx in range(n-k):
        kmer = dna[idx:idx+k]
        curr = 1.0
        for i, nucleotide in enumerate(kmer):
            curr *= probs[nucleotide][i]
        if curr > max_prob:
            result = kmer
            max_prob = curr
    return result



def create_profile_matrix(profile, pseudocount=True):
    """
    Create probability matrix

    Parameters
    ----------
    profile: list of dna strings
    pseudocount: set to True to add 1 to all counts

    Returns
    -------
    probs : probability matrix, as a dictionary in hte form
                probs['A'] = [0.3, 0.2, 0.1]

    """
    t = len(profile)
    n = len(profile[0])
    probs = defaultdict(list)
    count = Counter()
    
    if pseudocount:
        fudge = 1
    else:
        fudge = 0
    
    for idx in range(n):
        counter = Counter(dna[idx] for dna in profile)
        for nucleotide in ('A', 'C', 'G', 'T'):
            count = counter[nucleotide]
            probs[nucleotide].append( (count+fudge)/t )
    return probs



def greedy_motif_search(dna, k, pseudocount=True):
    """
    Finds the best profile matrix for a list of dna strings

    Parameters
    ----------
    dna : list of dna strings
    k : length of motif
    pseudocount : set to True to use pseudocounts


    Returns
    -------
    profile matrix of the best motifs found

    """
    best_motifs = [gene[:k] for gene in dna]    
    n = len(dna[0])
    
    for idx in range(n-k):
        kmer0 = dna[0][idx:idx+k]
        profile = [kmer0]
        for i in range(1, len(dna)):
            probs = create_profile_matrix( profile, pseudocount )
            probable_kmer = most_probable_kmer(dna[i], k, probs)
            profile.append(probable_kmer)
        best_score = profile_score(best_motifs)
        curr_score = profile_score(profile)
        if curr_score < best_score:
            best_motifs = profile
        
    return best_motifs