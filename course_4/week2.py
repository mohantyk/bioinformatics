from itertools import product, count
from math import inf

import numpy as np
import pandas as pd

import logging
logging.basicConfig()
logger = logging.getLogger('neighbor')
logger.setLevel(logging.ERROR)

class Tree:
    def __init__(self, adjacency=None):
        '''
        adjacency : graph (dict of dicts)
        '''
        self.adjacency = adjacency if adjacency is not None else {}

    def add_edge(self, node0, node1, weight=0):
        if node0 not in self.adjacency:
            self.adjacency[node0] = {}
        if node1 not in self.adjacency:
            self.adjacency[node1] = {}
        self.adjacency[node0][node1] = weight
        self.adjacency[node1][node0] = weight

    def del_edge(self, node0, node1):
        del self.adjacency[node0][node1]
        del self.adjacency[node1][node0]

    def add_node(self, node):
        if node not in self.adjacency:
            self.adjacency[node] = {}
        else:
            raise ValueError(f'Node {node} already exists')
        return node


def cluster_distance(cluster0, cluster1, distances):
    '''
    cluster0/cluster1 : sets
    distances : distance matrix
    '''
    size0 = len(cluster0)
    size1 = len(cluster1)
    total = 0
    for (i, j) in product(cluster0, cluster1):
        total += distances[i, j]
    final_distance = total/(size0*size1)
    return final_distance


def upgma(n, distances):
    clusters = set()
    age = {}
    tree = Tree()
    next_node_idx = 0
    cluster_to_node = {}

    for leaf in range(n):
        tree.add_node(leaf)
        cluster = frozenset({leaf})
        cluster_to_node[cluster] = leaf
        clusters.add(cluster)
        age[leaf] = 0
        next_node_idx += 1

    while len(clusters) > 1:
        # Look for clusters with minimum distance
        min_distance = inf
        for (cl0, cl1) in product(clusters, repeat=2):
            if cl0 == cl1: continue
            dist = cluster_distance(cl0, cl1, distances)
            if dist < min_distance:
                min_distance = dist
                min_clusters = {cl0, cl1}

        # Merge clusters with minimum distance
        cluster0 = min_clusters.pop(); clusters.remove(cluster0)
        cluster1 = min_clusters.pop(); clusters.remove(cluster1)
        new_cluster = cluster0.union(cluster1)
        clusters.add(new_cluster)
        # Add a new node for the new cluster
        new_node = tree.add_node(next_node_idx); next_node_idx += 1
        cluster_to_node[new_cluster] = new_node
        age[new_node] = min_distance/2
        # Connect new node to the child nodes
        child0 = cluster_to_node[cluster0]
        child1 = cluster_to_node[cluster1]
        tree.add_edge(new_node, child0, age[new_node] - age[child0])
        tree.add_edge(new_node, child1, age[new_node] - age[child1])

    return tree.adjacency

def create_d_star(distances):
    total_distance = distances.sum(axis=1)
    d_star = np.zeros_like(distances, dtype=float)
    shape = distances.shape
    n = shape[0]
    for i in range(shape[0]):
        for j in range(shape[1]):
            if i==j: continue
            d_star[i, j] = (n-2)*distances[i, j] - total_distance[i] - total_distance[j]
    return d_star


def neighbor_joining(distances, node_counter=None):
    n = len(distances)
    if node_counter is None:
        node_counter = count(n)

    if not isinstance(distances, pd.DataFrame):
        df = pd.DataFrame(data=distances, index=list(range(n)), columns=list(range(n)))
    else:
        df = distances

    if n == 2:
        leaf0, leaf1 = df.columns
        edge = df.to_numpy()[0,1]
        adjacency = {leaf0: {leaf1: edge}, leaf1: {leaf0: edge}}
        return Tree(adjacency)

    d_star = create_d_star(df.to_numpy())
    np.fill_diagonal(d_star, np.inf)
    (idx0, idx1) = np.unravel_index(np.argmin(d_star), d_star.shape)

    i = df.columns[idx0]
    j = df.columns[idx1]
    total_distance = df.sum(axis=1)
    delta = (total_distance[i] - total_distance[j]) / (n-2)

    limb_i = 0.5*(df.loc[i, j] + delta)
    limb_j = 0.5*(df.loc[i,j] - delta)
    new_node = next(node_counter)

    distance_from_new_node = [0.5*(df.loc[i, node] + df.loc[j, node] - df.loc[i, j])
                              for node in df.columns]
    df[new_node] = distance_from_new_node # Add new column
    df.loc[new_node] = distance_from_new_node + [0] # Add new row
    df = df.drop(index=[i,j], columns=[i,j]) # Drop old rows and columns

    tree = neighbor_joining(df, node_counter)
    tree.add_edge(new_node, i, limb_i)
    tree.add_edge(new_node, j, limb_j)
    return tree

