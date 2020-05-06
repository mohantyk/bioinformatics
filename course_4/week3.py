from collections import deque
from math import inf

import sys
sys.path.append('..')
from course_3.week5 import pairwise

alphabet = ['A', 'C', 'G', 'T']
class Node:
    ''' Node in a directed binary tree '''
    def __init__(self, val, left=None, right=None):
        self.val = val
        self.left = left
        self.right = right

    def is_leaf(self):
        return (not self.left) and (not self.right)

    def add_char(self, ch):
        '''
        Adds a character to the name (val) of the node
        '''
        if not self.is_leaf():
            self.val = self.val + ch

    def store_backtrack_data(self, data):
        '''
        data : dict of {char: (lchar, rchar)}
        '''
        self.backtrack_data = data

    def backtrack(self, ch):
        '''
        Adds character to name, and then backtracks to child nodes
        '''
        self.add_char(ch)
        if not self.is_leaf():
            lchar, rchar = self.backtrack_data[ch]
            self.left.backtrack(lchar)
            self.right.backtrack(rchar)


def create_binary_tree(values):
    '''
    Creates a binary tree
    input:
        values : array of values
    '''
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


def parsimony_backtrack(score, k):
    '''
    inputs:
        score: dictionary
        k : letter in alphabet
    output:
        curr_min : min score
        min_char : key in score that minimizes the score
    '''
    curr_min = inf
    for i in score:
        curr = score[i] + int(k!=i)
        if curr < curr_min:
            curr_min = curr
            min_char = i
    return curr_min, min_char


def small_parsimony_score(node, idx):
    '''
    Returns a dictionary of letter : score for each letter in alphabet
    '''
    score = {letter:inf for letter in alphabet}
    if node.is_leaf():
        ch = node.val[idx]
        score[ch] = 0
        return score

    lscore = small_parsimony_score(node.left, idx)
    rscore = small_parsimony_score(node.right, idx)
    backtrack = {}
    for k in alphabet:
        lmin, lchar = parsimony_backtrack(lscore, k)
        rmin, rchar = parsimony_backtrack(rscore, k)
        score[k] = min(lscore[i] + int(k!=i) for i in lscore) + min(rscore[j] + int(k!=j) for j in rscore)
        assert score[k] == lmin + rmin # Same thing calculated in two different ways
        backtrack[k] = (lchar, rchar)
        node.store_backtrack_data(backtrack)

    return score

def small_parsimony(root):
    '''
    Returns the final small parsimony score for the tree
    '''
    head = root
    while not head.is_leaf():
        head = head.left
    n = len(head.val)

    total_score = 0
    for idx in range(n): # idx-th tree
        score = small_parsimony_score(root, idx)
        ch = min(score, key=score.get)
        total_score += score[ch]
        root.backtrack(ch)
    return total_score
