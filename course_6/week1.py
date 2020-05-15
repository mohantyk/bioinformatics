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
        self.position = {}

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


class ModifiedSuffixTrie(Trie):
    def __init__(self, text):
        self.root = Node()
        self.add_text(text)

    def add_text(self, text):
        l = len(text)
        for i, _ in enumerate(text):
            curr = self.root
            for j in range(i,l):
                ch = text[j]
                if ch not in curr.children:
                    curr.children[ch] = Node()
                    curr.position[ch] = j
                curr = curr.children[ch]
            if curr.is_leaf():
                curr.val = i

class SuffixTree:
    def __init__(self, text):
        self.text = text
        mod_suffix_trie = ModifiedSuffixTrie(text)
        self.root = mod_suffix_trie.root
        self.replace_non_branching_paths(self.root)


    def replace_non_branching_paths(self, node):
        '''
        Replace any non-branching paths with a single edge
        '''
        edges = list(node.children)
        for ltr in edges: # Can not iterate over dict while changing it
            child = node.children[ltr]
            num_edges, final = self.find_non_branch_edges(child)
            pos = node.position[ltr]

            del node.children[ltr]
            node.children[(pos, num_edges+1)] = final
            self.replace_non_branching_paths(final)
        node.positions = {} # Reset positions, not needed

    @property
    def edges(self):
        '''
        Returns all edges in suffix tree
        '''
        return self._edges_at_node(self.root)

    def _edges_at_node(self, node):
        '''
        Returns all edges in sub-tree starting from node
        '''
        edges = []
        for (pos, num_edges), child in node.children.items():
            edges.append(self.text[pos:pos+num_edges])
            edges += self._edges_at_node(child)
        return edges



    def find_non_branch_edges(self, node):
        '''
        output:
            num_edges: number of edges on non-branching path starting from node
            final: final node on non-branching path
        '''
        curr = node
        path = []
        while len(curr.children) == 1:
            ltr = list(curr.children)[0]
            path.append(ltr)
            curr = curr.children[ltr]
        final = curr
        num_edges = len(path)
        return num_edges, final




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


if __name__ == '__main__':
    text = 'ATAAATG$'
    suffix_tree = SuffixTree(text)