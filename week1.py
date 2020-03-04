#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: kaniska
"""

def composition(dna, k):
    n = len(dna)
    kmers = []
    for idx in range(n-k+1):
        kmer = dna[idx:idx+k]
        kmers.append(kmer)
        
    return kmers
        