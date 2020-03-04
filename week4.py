#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: kaniska
"""
from random import randint, choice, choices
from math import inf

from week3 import create_profile_matrix, most_probable_kmer, profile_score

def create_random_motifs( dna, k ):
    """
    Create random motif
    """
    n = len(dna[0])
    
    motifs = []
    for curr in dna:
        idx = randint(0, n-k)
        kmer = curr[idx:idx+k]
        motifs.append( kmer )
    return motifs



def randomized_motif_search( dna, k ):
    best_motifs = create_random_motifs(dna, k)
    best_score = profile_score( best_motifs )
    
    while True:
        motifs = []
        profile_matrix = create_profile_matrix(best_motifs)
        for curr in dna:
            kmer = most_probable_kmer(curr, k, profile_matrix)
            motifs.append( kmer )
            
        score = profile_score( motifs )
        if score < best_score:
            best_motifs = motifs
            best_score = score
        else:
            return best_motifs
        
  

def kmer_probability( kmer, probs ):
    final = 1.0
    for i, nucleotide in enumerate(kmer):
        final *= probs[nucleotide][i]
    return final



def gibbs_sampler(dna, k, t, N):
    assert(len(dna) == t)
    n = len(dna[0])
    
    motifs = create_random_motifs(dna, k)
    best_motifs = motifs.copy()
    best_score = profile_score( best_motifs )
    
    for _ in range(N):
        i = choice(range(t))
        new_motif = motifs[:i] + motifs[i+1:]
        profile_matrix = create_profile_matrix(new_motif)
        
        curr = dna[i]
        probs = [kmer_probability(curr[idx:idx+k], profile_matrix) 
                 for idx in range(n-k)]   

        [idx] = choices(range(n-k), probs)
        kmer = curr[idx: idx+k]
        motifs[i] = kmer
        
        score = profile_score( motifs )
        if score < best_score:
            best_motifs = motifs.copy()
            best_score = score
            
    return best_motifs



def multiple_runs(num, func, *args):
    best_score = inf
    best_motifs = None
    
    for _ in range(num):
        motifs = func(*args)
        score = profile_score(motifs)
        
        if score < best_score:
            best_score = score
            best_motifs = motifs.copy()
            
    return best_motifs
