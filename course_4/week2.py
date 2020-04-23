from itertools import product

class Tree:
    def __init__(self, adjacency=None):
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
    for leaf in range(n):
        tree.add_node(leaf)
        clusters.add(frozenset({leaf}))
        age[leaf] = 0

    while len(clusters) > 1:
        break

