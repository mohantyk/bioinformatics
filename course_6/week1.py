from collections import defaultdict
import itertools

import sys
sys.path.append('..')
from course_4.week4 import DirectedGraph

###------------------------------
### Useful classes
class Node:
    def __init__(self, val=None):
        self.val = val
        self.children = {}

    def is_leaf(self):
        return len(self.children) == 0


class Trie:
    def __init__(self):
        self.counter = itertools.count(0)
        self.root = Node(next(self.counter))

    def add(self, str):
        curr = self.root
        for ch in str:
            nxt = curr.children.get(ch, None)
            if nxt is None:
                nxt = Node(next(self.counter))
                curr.children[ch] = nxt
            curr = nxt

    def match_prefix(self, text):
        curr = self.root
        for ch in text:
            if ch not in curr.children:
                return False
            curr = curr.children[ch]
            if curr.is_leaf():
                return True
        return False

    @property
    def adjacency(self):
        graph = self.create_adjacency(self.root)
        return graph.adjacency

    def create_adjacency(self, node, graph=None):
        if graph is None:
            graph = DirectedGraph()
        for ltr, child in node.children.items():
            graph.add_edge(node.val, child.val, ltr)
            self.create_adjacency(child, graph)
        return graph

###--------------------------------
### Functions
def create_trie(patterns):
    trie = Trie()
    for pattern in patterns:
        trie.add(pattern)
    return trie

def match_trie(text, patterns):
    trie = create_trie(patterns)
    matching_indices = []
    for idx, _ in enumerate(text):
        if trie.match_prefix(text[idx:]):
            matching_indices.append(idx)
    return matching_indices