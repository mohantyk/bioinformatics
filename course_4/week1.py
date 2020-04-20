import numpy as np
from math import inf
from itertools import combinations

import logging
logging.basicConfig()
logger = logging.getLogger('week1')
logger.setLevel(logging.ERROR)

def calculate_distances(n, adjacency):
    '''
    inputs:
        n: n is the number of leaf nodes (0...n-1)
        adjacency : weighted adjacency {node: {nghbr0: weight0, nghbr1: weight1}, ... }
    outputs:
        distance matrix
    '''
    # Create initial distance graph - single edges
    max_node = max(adjacency.keys())
    distances = np.zeros((max_node+1, max_node+1), int)
    for node in adjacency:
        for nghbr, weight in adjacency[node].items():
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

def find_path(adjacency, src, dest):
    '''
    weighted adjacency data structure : {node0 : {node1: distance1, node2: distance2 }}
    src, dest : nodes
    '''
    if dest in adjacency[src]:
        return [src, dest]

    for nghbr in adjacency[src]:
        back_to_src = adjacency[nghbr].pop(src) # Remove edge back to src
        new_path = find_path(adjacency, nghbr, dest)
        adjacency[nghbr][src] = back_to_src

        if new_path:
            return [src] + new_path

    return None # If no path found

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
        limb length
        numpy 2-D array of distances when node limb is set to zero
    '''
    num_leafs = len(distances)
    limb = limb_length(num_leafs, node, distances)
    bald = np.array(distances, dtype=float)
    bald[:, node] -= limb
    bald[node,:] -= limb
    bald[node, node] = 0
    return limb, bald

def find_insertion_end_points(node, distances):
    '''
    distances : should be a bald matrix
    '''
    valid_nodes = list(range(node)) + list(range(node+1, len(distances)))
    for (i, k) in combinations(valid_nodes, 2):
        if distances[i, k] == distances[i, node] + distances[k, node]:
            return (i, k)
    raise ValueError('Can not find insertion point')

def trim_distances(node, distances):
    '''
    Removes the column and row for node
    '''
    mask = np.ones(len(distances), dtype=bool)
    mask[node] = False
    trimmed = distances[mask,:][:, mask]
    return trimmed

def generate_new_node(num_leafs):
    '''Generator to yield new nodes'''
    new_node = num_leafs
    while True:
        yield new_node
        new_node += 1

def additive_phylogeny(num_leafs, distances, new_node_generator=None):
    '''
    output:
        weighted adjacency dict
    '''
    assert distances.shape == (num_leafs, num_leafs)
    # Trivial case
    if distances.shape == (2, 2):
        assert distances[0, 1] == distances[1, 0]
        weight = distances[0, 1]
        return {0: {1: weight},
                1: {0: weight} }

    # Create a generator to yield new nodes when needed
    if new_node_generator is None:
        new_node_generator = generate_new_node(num_leafs)

    node = num_leafs - 1
    limb, bald = create_bald_matrix(node, distances)
    # Find insertion point and distance
    (i, k) = find_insertion_end_points(node, bald)
    distance_from_i = bald[i, node]
    # Trim the matrix
    trimmed = trim_distances(node, bald)
    # Create a tree from the trimmed matrix recursively
    base_tree = additive_phylogeny(num_leafs-1, trimmed, new_node_generator)
    logger.debug(f'Obtained tree : {base_tree}')
    # Insert node in the tree
    path = find_path(base_tree, i, k)
    distance_to_node = 0
    for idx, curr_node in enumerate(path):
        if idx == 0: continue # Ignore the src
        prev_node = path[idx-1]
        prev_distance = distance_to_node
        distance_to_node += base_tree[prev_node][curr_node]
        if distance_to_node == distance_from_i:
            insertion_node = curr_node # Can insert at existing node
            break
        elif distance_to_node > distance_from_i:
            # Add a new node
            insertion_node = next(new_node_generator)
            base_tree[insertion_node] = {}
            # Break the current edge
            edge_weight = base_tree[prev_node].pop(curr_node)
            base_tree[curr_node].pop(prev_node)
            base_tree[prev_node][insertion_node] = distance_from_i - prev_distance
            base_tree[insertion_node][prev_node] = base_tree[prev_node][insertion_node]

            base_tree[insertion_node][curr_node] = edge_weight - base_tree[prev_node][insertion_node]
            base_tree[curr_node][insertion_node] = base_tree[insertion_node][curr_node]
            break
        else:
            continue
    logger.debug(f'Inserting at node {insertion_node} in path from {i} -> {k} at distance {distance_from_i} from {i}')
    base_tree[insertion_node][node] = limb
    base_tree[node] = {}
    base_tree[node][insertion_node] = limb
    return base_tree
