from pytest import approx

from week1 import *

class TestWeek1:
    def test_farthest_first_traversal(self):
        k = 3
        m = 2
        points = [(0, 0), (5, 5), (0, 5), (1, 1), (2, 2), (3, 3), (1, 2)]

        expected_centers  = {(0,0), (5,5), (0,5)}
        assert farthest_first_traversal(k, m, points) == expected_centers


    def test_squared_error_distortion(self):
        k = 2
        m = 2
        centers = [(2.31, 4.55), (5.96, 9.08)]
        points = [  (3.42, 6.03), (6.23, 8.25), (4.76, 1.64),
                    (4.47, 4.33), (3.95, 7.61), (8.93, 2.97),
                    (9.74, 4.03), (1.73, 1.28), (9.72, 5.01), (7.27, 3.77)]
        assert squared_error_distortion(k, m, centers, points) == approx(18.246, abs=1e-3)


    def test_lloyd(self):
        k = 2
        m = 2
        points = [(1.3, 1.1), (1.3, 0.2), (0.6, 2.8), (3.0, 3.2),
                    (1.2, 0.7), (1.4, 1.6), (1.2, 1.0), (1.2, 1.1),
                    (0.6, 1.5), (1.8, 2.6), (1.2, 1.3), (1.2, 1.0), (0.0, 1.9)]
        book_centers = {(1.800, 2.867), (1.060, 1.140)}
        approx_centers = {(1.044, 1.156), (1.8, 2.867)}
        assert lloyd(k, m, points) == approx_centers
