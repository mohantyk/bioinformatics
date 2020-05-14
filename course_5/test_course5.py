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
        assert squared_error_distortion(k, m, centers, points) == 18.246
