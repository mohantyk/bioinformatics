from math import sqrt, inf
import random


def euclidean_dist(p1, p2):
    assert len(p1) == len(p2)
    dist = sqrt(sum((x1-x2)**2 for x1, x2 in zip(p1, p2)))
    return dist

def center_of_gravity(points):
    n = len(points)
    m = len(points[0])
    center = []
    for idx in range(m):
        center.append(sum(pt[idx] for pt in points)/n)
    return tuple(center)


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
    assert len(points[0]) == m
    n = len(points)
    total_dist = 0
    for pt in points:
        dist = min(euclidean_dist(pt, center) for center in centers)**2
        total_dist += dist
    return total_dist/n


def lloyd(k, m, points):
    centers = set(points[:k])

    while True:
        clusters = {center: set() for center in centers}
        for pt in points:
            chosen = min(centers, key = lambda center: euclidean_dist(center, pt))
            clusters[chosen].add(pt)
        new_centers = {center_of_gravity(list(points)) for points in clusters.values()}
        if new_centers == centers:
            break
        else:
            centers = new_centers

    return list(centers)


