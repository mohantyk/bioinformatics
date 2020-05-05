from math import inf

alphabet = ['A', 'C', 'G', 'T']
class Node:
    def __init__(self, val, left=None, right=None):
        self.val = val
        self.left = left
        self.right = right

    def is_leaf(self):
        return (not self.left) and (not self.right)



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




