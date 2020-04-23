from itertools import product
from math import inf

class Tree:
    def __init__(self, adjacency=None):
        '''
        adjacency : graph (dict of dicts)
        '''
        self.adjacency = adjacency if adjacency is not None else {}

    def add_edge(self, node0, node1, weight):
        if node0 not in self.adjacency:
            self.adjacency[node0] = {}
        if node1 not in self.adjacency:
            self.adjacency[node1] = {}
        self.adjacency[node0][node1] = weight
        self.adjacency[node1][node0] = weight

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
        # Connect new node to the cluster nodes
        child0 = cluster_to_node[cluster0]
        child1 = cluster_to_node[cluster1]
        tree.add_edge(new_node, child0, min_distance/2)
        tree.add_edge(new_node, child1, min_distance/2)

    return tree.adjacency



