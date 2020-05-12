
from week2 import Tree

import sys
sys.path.append('..')
from course_2.week3 import AMINO_MASS, MASS_2_AMINO

class DirectedGraph(Tree):
    def add_edge(self, node0, node1, weight=0):
        if node0 not in self.adjacency:
            self.adjacency[node0] = {}
        self.adjacency[node0][node1] = weight

    def del_edge(self, node0, node1):
        del self.adjacency[node0][node1]


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


