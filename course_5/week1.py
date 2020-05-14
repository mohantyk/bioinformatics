from math import sqrt, inf
import random

def euclidean_dist(p1, p2):
    assert len(p1) == len(p2)
    dist = sqrt(sum((x1-x2)**2 for x1, x2 in zip(p1, p2)))
    return dist

def farthest_first_traversal(k, m, points):
    '''
    Returns centers of k clusters
    '''
    assert len(points[0]) == m
    centers = {points[0]}
    while len(centers) < k:
        farthest_dist = 0
        for pt in points:
            if pt in centers:
                continue
            dist = min(euclidean_dist(pt, center) for center in centers)
            if dist > farthest_dist:
                farthest_dist = dist
                farthest_pt = pt
        centers.add(farthest_pt)
    return centers

def squared_error_distortion(k, m, centers, points):
    pass
