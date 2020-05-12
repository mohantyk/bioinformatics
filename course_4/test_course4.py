import numpy as np
from numpy.testing import assert_array_equal
from pytest import approx
import pytest
from week1 import *
from week2 import *
from week3 import *
from week4 import *

class TestWeek1:
    def test_distances(self):
        adjacency = {   0: {4: 11},
                        1: {4: 2},
                        2: {5: 6},
                        3: {5: 7},
                        4: {0: 11, 1: 2, 5: 4},
                        5: {4: 4, 3: 7, 2: 6} }
        n = 4
        expected_distances = np.array([ [0, 13, 21, 22],
                                        [13, 0, 12, 13],
                                        [21, 12, 0, 13],
                                        [22, 13, 13, 0]])
        assert_array_equal(calculate_distances(n, adjacency), expected_distances)

    def test_limb_length(self):
        n = 4
        j = 1
        distances = [[0, 13, 21, 22],
                     [13, 0, 12, 13],
                     [21, 12, 0, 13],
                     [22, 13, 13, 0]]
        assert limb_length(n, j, distances) == 2

    def test_bald_matrix(self):
        j = 3
        distances = np.array([  [0, 13, 21, 22],
                                [13, 0, 12, 13],
                                [21, 12, 0, 13],
                                [22, 13, 13, 0]] )
        expected = np.array([[0, 13, 21, 15],
                         [13, 0, 12, 6],
                         [21, 12, 0, 6],
                         [15, 6, 6, 0]])
        limb, bald = create_bald_matrix(j, distances)
        assert limb == 7.0
        assert_array_equal(bald, expected)

    def test_insertion(self):
        j = 3
        bald = np.array([[0, 13, 21, 15],
                         [13, 0, 12, 6],
                         [21, 12, 0, 6],
                         [15, 6, 6, 0]])
        assert find_insertion_end_points(j, bald) == (0, 2)

    def test_trim(self):
        j = 3
        bald = np.array([[0, 13, 21, 15],
                         [13, 0, 12, 6],
                         [21, 12, 0, 6],
                         [15, 6, 6, 0]])

        trimmed = np.array([[0, 13, 21],
                            [13, 0, 12],
                            [21, 12, 0]])

        assert_array_equal(trim_distances(j, bald), trimmed)

    def test_trivial_additive_phylogeny(self):
        n = 2
        distances = np.array([[0, 42], [42, 0]])
        adjacency = {0: {1: 42},
                     1: {0: 42}}
        assert additive_phylogeny(n, distances) == adjacency

    def test_additive_phylogeny(self):
        n = 4
        distances = np.array([  [0, 13, 21, 22],
                                [13, 0, 12, 13],
                                [21, 12, 0, 13],
                                [22, 13, 13, 0]] )

        adjacency = {   0: {4: 11},
                        1: {4: 2},
                        2: {5: 6},
                        3: {5: 7},
                        4: {0: 11, 1: 2, 5: 4},
                        5: {4: 4, 3: 7, 2: 6} }
        assert additive_phylogeny(n, distances) == adjacency

    def test_path(self):
        adjacency = {   0: {4: 11},
                        1: {4: 2},
                        2: {5: 6},
                        3: {5: 7},
                        4: {0: 11, 1: 2, 5: 4},
                        5: {4: 4, 3: 7, 2: 6} }
        assert find_path(adjacency, 0, 3) == [0, 4, 5, 3]


class TestWeek2:

    def test_tree(self):
        adjacency = {0: {1: 10},
                     1: {0: 10}}
        tree = Tree(adjacency)
        tree.add_edge(0, 2, 20)
        final_adjacency = { 0: {1: 10, 2: 20},
                            1: {0: 10},
                            2: {0: 20}}
        assert tree.adjacency == final_adjacency

    def test_cluster_distance(self):
        cluster0 = {0,1}
        cluster1 = {2,3}
        distances = np.array([
            [0,3,4,3],
            [3,0,4,5],
            [4,4,0,2],
            [3,5,2,0]
        ])
        assert cluster_distance(cluster0, cluster1, distances) == 4


    def test_upgma(self):
        n = 4
        distances = np.array([  [0, 20, 17, 11],
                        [20, 0, 20, 13],
                        [17, 20, 0, 10],
                        [11, 13, 10, 0]])
        adjacency = {0: {5: 7.0},
                     1: {6: approx(8.833, abs=1e-3)},
                     2: {4: 5.0},
                     3: {4: 5.0},
                     4: {2: 5.0, 3: 5.0, 5: 2.0},
                     5: {0: 7.0, 4: 2.0, 6: approx(1.833, abs=1e-3)},
                     6: {5: approx(1.833, abs=1e-3), 1: approx(8.833, abs=1e-3)}}
        assert upgma(n, distances) == adjacency

    def test_upgma_example(self):
        n = 4
        distances = np.array([
            [0,3,4,3],
            [3,0,4,5],
            [4,4,0,2],
            [3,5,2,0]
        ])
        adjacency = {0: {5: 1.5},
                     1: {5: 1.5},
                     2: {4: 1.0},
                     3: {4: 1.0},
                     5: {0: 1.5, 1: 1.5, 6: 0.5},
                     4: {2: 1.0, 3: 1.0, 6: 1.0},
                     6: {4: 1.0, 5: 0.5}}
        assert upgma(n, distances) == adjacency

    def test_d_star(self):
        distances = np.array( [ [0, 13, 21, 22],
                        [13, 0, 12, 13],
                        [21, 12, 0, 13],
                        [22, 13, 13, 0]])
        d_star = np.array([[  0, -68, -60, -60],
                            [-68,   0, -60, -60],
                            [-60, -60,   0, -68],
                            [-60, -60, -68,   0]])

        assert_array_equal(create_d_star(distances), d_star)

    def test_neighbor_simple(self):
        distances = np.array([[0, 4], [4, 0]])
        tree = neighbor_joining(distances)
        assert tree.adjacency == {0: {1: 4}, 1: {0: 4}}


    def test_neighbor_joining(self):
        distances = np.array([  [0, 23, 27, 20],
                                [23, 0, 30, 28],
                                [27, 30, 0, 30],
                                [20, 28, 30, 0]])
        adjacency = {0: {4: 8.0},
                     1: {5: 13.50},
                     2: {5: 16.50},
                     3: {4: 12.0},
                     4: {0: 8.0, 3: 12.0, 5: 2.0},
                     5: {1: 13.5, 2: 16.5, 4: 2.0}}
        tree = neighbor_joining(distances)
        assert tree.adjacency == adjacency


class TestWeek3:

    def test_create_binary_tree(self):
        values = [0, 1, 2, 3, 4, 5, 6]
        root = create_binary_tree(values)
        assert root.right.left.val == 5

    def test_small_parsimony_score(self):
        values = ['', '', '', '', '', '', '', 'C', 'C', 'A', 'C', 'G', 'G', 'T', 'C']
        root = create_binary_tree(values)
        score = small_parsimony_score(root, 0)
        assert score == {'A': 5, 'C': 3, 'G': 4, 'T': 4}

    def test_small_parsimony(self):
        values = ['', '', '', 'CAAATCCC', 'ATTGCGAC', 'CTGCGCTG', 'ATGGACGA']
        root = create_binary_tree(values)
        assert small_parsimony(root) == 16
        assert root.val == 'ATAGACAA'

    def test_nearest_neighbor(self):
        graph = Tree()
        graph.add_edge(0, 4)
        graph.add_edge(1, 4)
        graph.add_edge(2, 5)
        graph.add_edge(3, 5)
        graph.add_edge(4, 5)

        nghbr0 = Tree()
        nghbr0.add_edge(1, 4)
        nghbr0.add_edge(0, 5)
        nghbr0.add_edge(3, 4)
        nghbr0.add_edge(2, 5)
        nghbr0.add_edge(4, 5)

        nghbr1 = Tree()
        nghbr1.add_edge(1, 5)
        nghbr1.add_edge(0, 4)
        nghbr1.add_edge(3, 4)
        nghbr1.add_edge(2, 5)
        nghbr1.add_edge(4, 5)

        nghbrs = nearest_tree_neighbor(graph, 5, 4)
        assert nghbrs == [nghbr0, nghbr1]


    def test_unrooted_small_parsimony(self):
        data = '''GCAGGGTA->5
                TTTACGCG->5
                CGACCTGA->6
                GATTCCAC->6
                5->TTTACGCG
                5->GCAGGGTA
                5->7
                TCCGTAGT->7
                7->5
                7->6
                7->TCCGTAGT
                6->GATTCCAC
                6->CGACCTGA
                6->7'''
        score, _ = solve_unrooted_small_parsimony(data.splitlines())
        assert score == 23

    @pytest.mark.skip(reason='Unstable. Fails sometimes, passes sometimes')
    def test_large_parsimony(self):
        data = '''5
                GCAGGGTA->5
                TTTACGCG->5
                CGACCTGA->6
                GATTCCAC->6
                5->TTTACGCG
                5->GCAGGGTA
                5->7
                TCCGTAGT->7
                7->5
                7->6
                7->TCCGTAGT
                6->GATTCCAC
                6->CGACCTGA
                6->7'''
        results = large_parsimony_interchange_heuristic(data.splitlines())
        scores = [step[0] for step in results]
        assert scores == [22, 21]


class TestWeek4:
    def test_spectrum_graph(self):
        spectrum = [57, 71, 154, 185, 301, 332, 415, 429, 486]
        graph = graph_from_spectrum(spectrum)
        expected = {   0: {57: 'G', 71: 'A'},
                            57: {154: 'P', 185: 'Q'},
                            71: {185: 'N'},
                            154: {301: 'F'},
                            185: {332: 'F'},
                            301: {415: 'N', 429: 'Q'},
                            332: {429: 'P'},
                            415: {486: 'A'},
                            429: {486: 'G'}
                        }
        assert graph.adjacency == expected

    def test_ideal_spectrum(self):
        assert ideal_spectrum('GPG') == [0, 57, 57, 154, 154, 211]

    def test_find_all_paths(self):
        graph = DirectedGraph({   0: {57: 'G', 71: 'A'},
                            57: {154: 'P', 185: 'Q'},
                            71: {185: 'N'},
                            154: {301: 'F'},
                            185: {332: 'F'},
                            301: {415: 'N', 429: 'Q'},
                            332: {429: 'P'},
                            415: {486: 'A'},
                            429: {486: 'G'} })
        assert find_all_paths(301, 486, graph) == [[301, 415, 486], [301, 429, 486]]


    def test_decode_ideal_spectrum(self):
        spectrum = [57, 71, 154, 185, 301, 332, 415, 429, 486]
        assert decode_ideal_spectrum(spectrum) == 'GPFNA'

    def test_peptide_sequencing(self):
        spectral_vector = [0,0,0,4,-2,-3,-1,-7,6,5,3,2,1,9,3,-8,0,3,1,2,1,8]
        example_amino = {'X': 4, 'Z': 5}
        assert peptide_sequencing(spectral_vector, example_amino) == 'XZZXX'

