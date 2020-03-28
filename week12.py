#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: kaniska
"""

from collections import defaultdict, deque, Counter
from itertools import product
from copy import deepcopy
from random import choice

import logging
logging.basicConfig()
my_logger = logging.getLogger('MyLogger')
my_logger.setLevel(logging.ERROR)


# Helper functions
def prefix(kmer):
    if isinstance(kmer, tuple): # Read pair
        return tuple(prefix(pattern) for pattern in kmer)
    return kmer[:-1]

def suffix(kmer):
    if isinstance(kmer, tuple): # Read pair
        return tuple(suffix(pattern) for pattern in kmer)
    return kmer[1:]


def composition(dna, k):
    n = len(dna)
    kmers = []
    for idx in range(n-k+1):
        kmer = dna[idx:idx+k]
        kmers.append(kmer)
    return kmers


def read_pairs(genome, k, d):
    reads = []
    n = len(genome)
    for i in range(n-(2*k+d)+1):
        pat1 = genome[i:i+k]
        pat2 = genome[i+k+d:i+2*k+d]
        reads.append((pat1, pat2))
    return reads


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


def euler_cycle(adjacency_list):
    start = choice( list(adjacency_list.keys()) )
    adjacency = deepcopy(adjacency_list) # Copy adjacency matrix
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
            except (IndexError, KeyError): # No neighbors left
                break

        if len(visited) == len(adjacency): # All done
            break

        path.pop()  # Remove the end node ( same as start node )
        start = path[0]
        # Rotate cycle until we find a node with outgoing paths
        while not len( adjacency[start]):
            path.rotate(-1)
            start = path[0]
        path.append(start) # Complete cycle

    return path


def node_degrees(adjacency):
    """
    Calculates in- and out-degrees of each node

    Parameters
    ----------
    adjacency : adjacency list

    Returns
    -------
    degrees : dictionary of {node: (in, out)} pairs

    """
    incoming = Counter()
    out = {}
    for k, v in adjacency.items():
        incoming.update(v)
        out[k] = len(v)

    for k in (out.keys() - incoming.keys()):
        incoming[k] = 0
    for k in incoming.keys() - out.keys():
        out[k] = 0

    degrees = {k: (incoming[k], out[k]) for k in incoming}
    return degrees


def euler_path(adjacency):
    degrees = node_degrees(adjacency)
    unbalanced = set()
    start = []
    stop = []
    for k, (inc, out) in degrees.items():
        if inc != out :
            unbalanced.add(k)
            if inc > out :
                stop.append(k)
            else:
                start.append(k)

    my_logger.debug(f'Unbalanced nodes: {unbalanced}')
    if len(unbalanced) not in (0,2):
        raise ValueError(f'Number of unbalanced nodes is {len(unbalanced)}.')

    balanced = deepcopy(adjacency)
    if unbalanced:
        begin = start[0]
        end = stop[0]
        if end not in balanced:
            balanced[end] = []
        balanced[end].append( begin )

    path = euler_cycle(balanced)
    if unbalanced:
        path.pop()
        while not (path[0]==begin and path[-1]==end):
            path.rotate(-1)

    return path


def genome_from_path(path, d=None):
    path = list(path)
    if not isinstance(path[0], tuple): # Simple kmers
        genome = path[0]
        for node in path[1:]:
            genome += node[-1]
    else :                              # Paired reads, each node in the path is a tuple
        # TODO: Check if the path is a valid solution
        genome = path[0][0]             # Pattern1 of first node
        for node in path[1:]:
            genome += node[0][-1]
        my_logger.debug(f'Genome constructed only from pattern1 : {genome}')
        for node in path[-(d+2):-1]: # Fill in d missing nucleotides in the final pair
            genome += node[1][0]
        my_logger.debug(f'Genome after filling in the missing parts of the final read : {genome}')
        genome += path[-1][1]
    return genome


def find_contigs(graph):
    interior = []
    degrees = node_degrees(graph)
    for node, deg in degrees.items():
        if deg == (1,1):
            interior.append(node)

    paths = []
    graph = deepcopy(graph)
    for start_node in graph:
        while graph[start_node]:
            node = start_node
            path = [node]
            while True:
                nxt = graph[node].pop()
                path.append(nxt)
                node = nxt
                if (nxt not in interior) or (not graph[node]):
                    paths.append(path)
                    break
    my_logger.debug(f'Paths for contigs: {paths}')
    contigs = [genome_from_path(path) for path in paths]
    return contigs
