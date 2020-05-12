from collections import deque
from math import inf

from week2 import Tree

import sys
sys.path.append('..')
from course_2.week3 import AMINO_MASS, MASS_2_AMINO
from course_3.week5 import pairwise

class DirectedGraph(Tree):
    def add_edge(self, node0, node1, edge_weight=0):
        if node0 not in self.adjacency:
            self.adjacency[node0] = {}
        self.adjacency[node0][node1] = edge_weight

    def del_edge(self, node0, node1):
        del self.adjacency[node0][node1]

    def add_weight(self, node, weight):
        try:
            self.weights[node] = weight
        except AttributeError:
            self.weights = {}
            self.weights[node] = weight

    def incoming_edges(self):
        incoming = {}
        for node in self.adjacency:
            for nghbr, edge in self.adjacency[node].items():
                if nghbr not in incoming:
                    incoming[nghbr] = {}
                incoming[nghbr][node] = edge
        return incoming

def graph_from_spectrum(spectrum):
    graph = DirectedGraph()
    if spectrum[0] != 0:
        spectrum = [0] + spectrum
    for idx, mass in enumerate(spectrum):
        for next_mass in spectrum[idx:]:
            diff = next_mass - mass
            if diff in MASS_2_AMINO:
                graph.add_edge(mass, next_mass, MASS_2_AMINO[diff])
    return graph


def calculate_mass(peptide):
    return sum(AMINO_MASS[amino] for amino in peptide)

def ideal_spectrum(peptide):
    n = len(peptide)
    masses = []
    for idx in range(n):
        prefix = peptide[0:idx]
        masses.append(calculate_mass(prefix))
    for idx in range(n):
        suffix = peptide[idx:]
        masses.append(calculate_mass(suffix))
    return sorted(masses)


def find_all_paths(source, sink, graph):
    if source == sink:
        return [[sink]]
    paths = []
    adjacency = graph.adjacency
    for nghbr in adjacency[source]:
        paths_from_nghbr = find_all_paths(nghbr, sink, graph)
        for path in paths_from_nghbr:
            paths.append([source] + path)
    return paths


def decode_ideal_spectrum(spectrum):
    graph = graph_from_spectrum(spectrum)
    paths = find_all_paths(0, max(spectrum), graph)
    for path in paths:
        peptide_elems = []
        for idx, src in enumerate(path[:-1]):
            nghbr = path[idx+1]
            peptide_elems.append(graph.adjacency[src][nghbr])
        peptide = ''.join(peptide_elems)
        if ideal_spectrum(peptide)[1:] == spectrum: # Don't include 0 mass
            return peptide


def peptide_vector(peptide):
    vector = []
    for amino in peptide:
        mass = AMINO_MASS[amino]
        amino_vec = [0]*mass
        amino_vec[-1] = 1
        vector.extend(amino_vec)
    return vector

def vec_to_peptide(vector):
    start = -1
    masses = []
    for idx, val in enumerate(vector):
        if val == 1:
            masses.append(idx-start)
            start = idx
    return ''.join(MASS_2_AMINO[mass] for mass in masses)

def peptide_sequencing(spectral_vec, mass_to_amino=MASS_2_AMINO):
    def score(path, graph):
        total = sum(graph.weights[node] for node in path)
        return total

    graph = DirectedGraph()
    graph.add_node(0)
    graph.add_weight(0, 0)
    for idx, s_i in enumerate(spectral_vec):
        i = idx+1
        graph.add_node(i)
        graph.add_weight(i, s_i)
        for mass in mass_to_amino:
            j = i - mass
            if j >= 0:
                graph.add_edge(j, i, mass_to_amino[mass])
    m = i

    # Dynamic programming
    incoming = graph.incoming_edges()
    best_score = {m: 0}
    best_path = {m: [m]}
    to_process = deque([m])
    while to_process:
        node = to_process.popleft()
        node_weight = graph.weights[node]
        for nghbr in incoming.get(node, {}):
            new_path = best_score[node] + node_weight
            current_best = best_score.get(node, -inf)
            if new_path > current_best:
                best_score[nghbr] = new_path
                to_process.append(nghbr)
                best_path[nghbr] = [nghbr] + best_path[node]

    final_path = best_path[0]
    peptide = []
    for idx, node in enumerate(final_path[:-1]):
        nxt = final_path[idx+1]
        peptide.append(graph.adjacency[node][nxt])
    return ''.join(peptide)