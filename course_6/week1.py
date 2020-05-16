from collections import defaultdict
from math import inf
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
        self.color = None

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
        '''
        Children keys are now (pos, len) tuples
        '''
        if text[-1] != '$': # Use $ as a text-end marker
            text += '$'
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

    def longest_repeat(self):
        '''
        Find longest repeated substring
        '''
        return self._longest_path_to_fork(self.root)


    def _longest_path_to_fork(self, node):
        '''
        Find path to farthest fork from node
        '''
        longest = -inf
        best_path = ''
        for (pos, num_edges), child in node.children.items():
            edge = self.text[pos:pos+num_edges]
            if not child.is_leaf():
                child_path = self._longest_path_to_fork(child)
                path_len = num_edges + len(child_path)
                if path_len > longest:
                    longest = path_len
                    best_path = edge + child_path
        return best_path


class ColoredSuffixTree(SuffixTree):
    def color_nodes(self, node=None):
        sep = self.text.index('#') # Use '#' as separator between two texts
        assert sep != -1

        if node == None:
            node = self.root

        if node.is_leaf():
            node.color = 'blue' if node.val <= sep else 'red'
            return node.color

        child_colors = set()
        for child in node.children.values():
            if child.color is None:
                self.color_nodes(child)
            child_colors.add(child.color)

        assert None not in child_colors
        if child_colors == {'blue'}:
            node.color = 'blue'
        elif child_colors == {'red'}:
            node.color = 'red'
        else:
            node.color = 'purple'
        return node.color






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


def longest_shared_substring(text1, text2):
    concatenated = text1 + '#' + text2 + '$'
    suffix_tree = SuffixTree(concatenated)

