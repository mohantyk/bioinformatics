from collections import deque
from math import inf

import sys
sys.path.append('..')
from course_3.week5 import pairwise

alphabet = ['A', 'C', 'G', 'T']
class Node:
    def __init__(self, val, left=None, right=None):
        self.val = val
        self.left = left
        self.right = right

    def is_leaf(self):
        return (not self.left) and (not self.right)


def create_binary_tree(values):
    root = Node(values[0])
    nodes = deque([root])
    level = 1
    while True:
        indices = slice(2**level-1, 2**(level+1)-1)
        if indices.stop > len(values):
            break

        for lval, rval in pairwise(values[indices]):
            node = nodes.popleft()
            node.left = Node(lval)
            nodes.append(node.left)

            node.right = Node(rval)
            nodes.append(node.right)
        level += 1

    return root


def small_parsimony_score(node, idx):
    score = {letter:inf for letter in alphabet}
    if node.is_leaf():
        ch = node.val[idx]
        score[ch] = 0
        return score

    lscore = small_parsimony_score(node.left, idx)
    rscore = small_parsimony_score(node.right, idx)
    for k in alphabet:
        score[k] = min(lscore[i] + k!=i for i in lscore) + min(rscore[j] + k!=j for j in lscore)
    return score




