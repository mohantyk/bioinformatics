import numpy as np
from math import inf
from itertools import combinations

def distances(n, adjacency):
    '''
    inputs:
        n: n is the number of leaf nodes (0...n-1)
        adjacency : weighted adjacency matrix node: [(nghbr0, weight0), (nghbr1, weight1)]
    outputs:
        distance matrix
    '''
    # Create initial distance graph - single edges
    max_node = max(adjacency.keys())
    distances = np.zeros((max_node+1, max_node+1), int)
    for node in adjacency:
        for (nghbr, weight) in adjacency[node]:
            distances[node, nghbr] = weight

    # Iteratively remove node and update weights to its neighbors
    for node in reversed(range(n, max_node+1)):
        nghbrs = np.where(distances[node,:])[0] # Find all its neighbors
        for nghbr in nghbrs:
            for idx in range(node):
                if distances[nghbr, idx] == 0 and idx in nghbrs: # Update only if weight is unknown
                    distances[nghbr, idx] = distances[node, idx] + distances[nghbr, node]
        distances = distances[:node, :node] # Remove the node
        np.fill_diagonal(distances, 0)

    return distances


def limb_length(num_leafs, node, distances):
    '''
    inputs:
        num_leafs: number of leaf nodes
        node: index of node whose limb length is to be calculated
        distances: numpy array with distances between each node pair
    outputs:
        limb length of leaf node
    '''
    limb_length = inf
    for i in range(num_leafs):
        for k in range(num_leafs):
            if i == node or k == node:
                continue
            value = (distances[i][node] + distances[node][k] - distances[i][k])/2
            if value <= limb_length:
                limb_length = value
    return limb_length

def create_bald_matrix(node, distances):
    '''
    inputs:
        node: index of node whose limb length is to be calculated
        distances: numpy array with distances between each node pair
    outputs:
        numpy 2-D array of distances when node limb is set to zero
    '''
    num_leafs = len(distances)
    limb = limb_length(num_leafs, node, distances)
    bald = np.array(distances, dtype=float)
    bald[:, node] -= limb
    bald[node,:] -= limb
    bald[node, node] = 0
    return bald

def find_insertion_end_points(node, distances):
    '''
    distances : should be a bald matrix
    '''
    valid_nodes = list(range(node)) + list(range(node+1, len(distances)))
    for (i, k) in combinations(valid_nodes, 2):
        if distances[i, k] == distances[i, node] + distances[k, node]:
            return (i, k)
    raise ValueError('Can not find insertion point')


