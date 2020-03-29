#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: kaniska
"""

from collections import defaultdict

def find_kmer(dna, kmer):
    """
    Find all indices where kmer occurs in dna

    Parameters
    ----------
    dna : dna string
    kmer : kmer string

    Returns
    -------
    list of indices

    """
    k = len(kmer)
    n = len(dna)
    found = []
    for idx in range(n-k+1):
        window = dna[idx:idx+k]
        if window == kmer:
            found.append(idx)
    return found



def kmer_frequency(dna, k):
    mycounter = defaultdict(int)
    n = len(dna)
    for idx in range(n-k+1):
        kmer = dna[idx:idx+k]
        mycounter[kmer] += 1
    return mycounter



def most_frequent_kmer(dna, k):
    mycounter = kmer_frequency(dna, k)
    max_count = max(mycounter.values())
    max_kmers = []
    for kmer, count in mycounter.items():
        if count == max_count:
            max_kmers.append(kmer)
    return max_kmers



def count_kmer(dna, kmer):
    k = len(kmer)
    n = len(dna)
    count = 0
    for idx in range(n-k+1):
        window = dna[idx:idx+k]
        if window == kmer:
            count += 1
    return count



def filter_counts(counter, threshold):
    kmers = set()
    for kmer, count in counter.items():
        if count >= threshold:
            kmers.add(kmer)
    return kmers



def find_clumps(dna, k, L, t):
    n = len(dna)
    initial_counter = kmer_frequency(dna[:L], k)
    final = filter_counts(initial_counter, t)
    
    for idx in range(n-L-1):
        remove_kmer = dna[idx: idx+k]
        initial_counter[remove_kmer] -= 1
        
        new_end = idx + L + 1
        add_kmer = dna[new_end-k:new_end]
        initial_counter[add_kmer] += 1
        
        if initial_counter[add_kmer] >= t:
            final.add( add_kmer )
    return final
    