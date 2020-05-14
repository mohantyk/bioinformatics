from week1 import *

class TestWeek1:
    def test_farthest_first_traversal(self):
        k = 3
        m = 2
        points = [(0, 0), (5, 5), (0, 5), (1, 1), (2, 2), (3, 3), (1, 2)]

        expected_centers  = {(0,0), (5,5), (0,5)}
        assert farthest_first_traversal(k, m, points) == expected_centers