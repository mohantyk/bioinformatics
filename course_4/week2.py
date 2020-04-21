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


def upgma(n, distances):
    pass