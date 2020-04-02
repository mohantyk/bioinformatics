from math import inf
import numpy as np

def dp_change(money, coins):
    '''
    Dynamic programming implementation for calculating minimum coins
    to make change for money.
    '''
    results = [inf]*(money+1)
    results[0] = 0
    for idx in range(1, money+1):
        for coin in coins:
            if coin <= idx:
                result = 1+results[idx-coin]
                if result < results[idx]:
                    results[idx] = result
    return results[money]


def manhattan_tourist(n, m, down_matrix, right_matrix, diag_matrix=None):
    assert right_matrix.shape == (n+1, m)
    assert down_matrix.shape == (n, m+1)
    assert diag_matrix is None or diag_matrix.shape == (n, m)

    longest_path = np.zeros((n+1, m+1), dtype=int)
    for j in range(1,m+1):
        longest_path[0, j] = longest_path[0, j-1] + right_matrix[0, j-1]
    for i in range(1, n+1):
        longest_path[i, 0] = longest_path[i-1, 0] + down_matrix[i-1, 0]
    for i in range(1, n+1):
        for j in range(1, m+1):
            from_left = longest_path[i, j-1] + right_matrix[i, j-1]
            from_top = longest_path[i-1, j] + down_matrix[i-1, j]
            if diag_matrix is not None:
                from_corner = longest_path[i-1, j-1] + diag_matrix[i-1, j-1]
            else:
                from_corner = 0

            longest_path[i, j] = max(from_left, from_top, from_corner)

    return longest_path[n, m]
