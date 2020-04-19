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
        assert_array_equal(distances(n, adjacency), expected_distances)

    def test_limb_length(self):
        n = 4
        j = 1
        distances = [[0, 13, 21, 22],
                     [13, 0, 12, 13],
                     [21, 12, 0, 13],
                     [22, 13, 13, 0]]
        assert limb_length(n, j, distances) == 2
