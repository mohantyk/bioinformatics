#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: kaniska
"""

from collections import defaultdict, deque
from itertools import product


# Helper functions
def prefix(kmer):
    return kmer[:-1]

def suffix(kmer):
    return kmer[1:]

def adjacency_to_file(adjacency, filename):
    with open(filename, 'w') as f:        
        for k, v in adjacency.items():
            f.write(f'{k} -> {", ".join(v)}\n')
            
            

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



def create_overlap_graph(kmers):
    adjacency = defaultdict( list )

    n = len(kmers)
    for i in range(n):
        head = kmers[i]
        for j in range(n):
            if i == j: continue
            tail = kmers[j]
            if suffix(head) == prefix(tail):
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


def de_bruijn_from_kmers( kmers ):
    adjacency = defaultdict(list)
    for kmer in kmers:
        adjacency[ prefix(kmer) ].append( suffix(kmer) )
    return adjacency



def euler_cycle(adjacency):
    start = 0
    path = deque([start])
    visited = {start}    
    
    while True:     
        next_node = adjacency[start].pop()
        
        while True:  # Add a new cycle
            path.append(next_node)
            visited.add(next_node)
            curr = path[-1]
            try:
                next_node = adjacency[curr].pop()
            except IndexError: # No neighbors left
                break
    
        if len(visited) == len(adjacency): # All done
            break
        
        path.pop()  # Remove the end node ( same as start node )
        start = path[0]
        # Rotate cycle until we find a node with outgoing paths
        while not len( adjacency[start]): 
            node = path.popleft()
            path.append(node)
            start = path[0]
        path.append( start ) # Complete cycle
        
    return path
        
        
