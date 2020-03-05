#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: kaniska
"""

from collections import defaultdict
from itertools import product

def composition(dna, k):
    n = len(dna)
    kmers = []
    for idx in range(n-k+1):
        kmer = dna[idx:idx+k]
        kmers.append(kmer)
        
    return kmers
        

def str_from_graph(kmers):
    for idx, kmer in enumerate(kmers):
        if idx == 0:
            genome = kmer
        else:
            genome += kmer[-1]
            
    return genome



def suffix(kmer):
    return kmer[:-1]

def prefix(kmer):
    return kmer[1:]



def create_overlap_graph(kmers):
    adjacency = defaultdict( list )

    n = len(kmers)
    for i in range(n):
        head = kmers[i]
        for j in range(n):
            if i == j: continue
            tail = kmers[j]
            if prefix(head) == suffix(tail):
                adjacency[head].append(tail)
    
    return adjacency



def universal_k_str(k):
    n = k + 2**k-1
    for combo in product(('0', '1'), repeat=n):
        bstr = ''.join(combo)
        final = set()
        for i in range(n-k+1):
            final.add( bstr[i:i+k])
            if len(final) == 2**k:
                return bstr
            


def de_bruijn_graph(dna, k):
    adjacency = defaultdict(list)
    prev = None
    n = len(dna)
    for i in range(n-(k-1)+1):
        kmer = dna[i:i+k-1]
        if prev is not None:
            adjacency[prev].append(kmer)
        prev = kmer
    return adjacency

