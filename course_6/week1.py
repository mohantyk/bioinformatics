from collections import defaultdict
import itertools

###------------------------------
### Useful classes
class Node:
    def __init__(self, val=None):
        self.val = val
        self.children = {}


class Trie:
    def __init__(self):
        self.counter = itertools.count(0)
        self.root = Node(next(self.counter))

    def add(self, str):
        curr = self.root
        for ch in str:
            nxt = curr.children.get(ch, None)
            if nxt is None:
                nxt = Node(next(self.counter))
                curr.children[ch] = nxt
            curr = nxt

###--------------------------------
### Functions
def create_trie(patterns):
    trie = Trie()
    for pattern in patterns:
        trie.add(pattern)
    return trie