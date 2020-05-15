from week1 import *

class TestWeek1:
    def test_trie_creation(self):
        patterns = ['ATAGA', 'ATC', 'GAT']
        trie = create_trie(patterns)
        adjacency = {   0: {1: 'A', 7: 'G'},
                        1: {2: 'T'},
                        2: {3: 'A', 6: 'C'},
                        3: {4: 'G'},
                        4: {5: 'A'},
                        7: {8: 'A'},
                        8: {9: 'T'} }
        assert trie.adjacency == adjacency

    def test_trie_matching(self):
        text = 'AATCGGGTTCAATCGGGGT'
        patterns = ['ATCG','GGGT']
        matching_indices = [1, 4, 11, 15]
        assert match_trie(text, patterns) == matching_indices