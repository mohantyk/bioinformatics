#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 20:36:08 2020

@author: kaniska
"""
from collections import Counter

NUCLEOTIDES = {'A', 'C', 'G', 'T'}
COMPLEMENTS = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

def skew(genome):
    result = [0]   
    for n in genome:
        if n=='G':
            value = result[-1] + 1
        elif n=='C':
            value = result[-1] - 1
        else:
            value = result[-1]
        result.append(value)   
    return result
      


def min_skew(genome):
    skew_values = skew(genome)
    min_val = min(skew_values)
    indices =  []
    for idx, val in enumerate(skew_values):
        if val == min_val:
            indices.append(idx)
    return indices
            

def hamming(g1, g2):
    if len(g1) != len(g2):
        raise ValueError
    
    count = 0
    for n1, n2 in zip(g1, g2):
        count += (n1 != n2)
    return count



def approx_pattern_match(pattern, text, d):
    k = len(pattern)
    n = len(text) 
    indices = []
    for idx in range(n):
        if idx + k > n: break
        window = text[idx:idx+k]
        if hamming(pattern, window) <= d:
            indices.append(idx)
            
    return indices


def neighbors(pattern, d):
    if d == 0:
        return {pattern}
    if len(pattern) == 1:
        return NUCLEOTIDES
    
    first = pattern[0]
    suffix = pattern[1:]
    neighborhood = set()
    suffix_neighbors = neighbors(suffix, d) 
    for neighbor in suffix_neighbors:
        if hamming(suffix, neighbor) == d:
            neighborhood.add(first + neighbor)
        else:
            for nucleotide in NUCLEOTIDES:
                neighborhood.add(nucleotide + neighbor)
    return neighborhood



def reverse_complement(pattern):
    table = str.maketrans(COMPLEMENTS)
    complement_pattern = pattern.translate(table)
    return complement_pattern[::-1]


def approximate_pattern_count(text, k, d, check_complements=False):
    count = Counter()
    n = len(text)

    for idx in range(n-k):
        kmer = text[idx:idx+k]
        kmer_neighbors = neighbors(kmer, d)
        count.update(kmer_neighbors)
        
        if check_complements:
            reversed_kmer = reverse_complement(kmer)
            reversed_neighbors = neighbors(reversed_kmer, d)
            count.update(reversed_neighbors)
    
    max_count = max(count.values())
    final = [kmer for kmer in count if count[kmer]==max_count]
    return final
