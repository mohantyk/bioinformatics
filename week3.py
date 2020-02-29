#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 29 13:05:42 2020

@author: kaniska
"""

from week2 import approx_pattern_match, neighbors


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