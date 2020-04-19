from numpy.testing import assert_array_equal
from week1 import *

class TestWeek1:
    def test_distances(self):
        adjacency = {0: [(4, 11)],
                    1: [(4, 2)],
                    2: [(5, 6)],
                    3: [(5, 7)],
                    4: [(0, 11), (1, 2), (5, 4)],
                    5: [(4, 4), (3, 7), (2, 6)]}
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
        adjacency = {0: [(1, 42)],
                     1: [(0, 42)]}
        assert additive_phylogeny(n, distances) == adjacency

    def test_additive_phylogeny(self):
        n = 4
        distances = np.array([  [0, 13, 21, 22],
                                [13, 0, 12, 13],
                                [21, 12, 0, 13],
                                [22, 13, 13, 0]] )

        adjacency = {   0: [(4, 11)],
                        1: [(4, 2)],
                        2: [(5, 6)],
                        3: [(5, 7)],
                        4: [(0, 11), (1, 2), (5, 4)],
                        5: [(4, 4), (3, 7), (2, 6)]  }
        assert additive_phylogeny(n, distances) == adjacency