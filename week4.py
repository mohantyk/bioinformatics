#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: kaniska
"""
from random import randint

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
        profile_matrix = create_profile_matrix(best_motifs, True)
        for curr in dna:
            kmer = most_probable_kmer(curr, k, profile_matrix)
            motifs.append( kmer )
            
        score = profile_score( motifs )
        if score < best_score:
            best_motifs = motifs
            best_score = score
        else:
            return best_motifs
        
        

        