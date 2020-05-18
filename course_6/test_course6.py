from collections import Counter

from week1 import *
from week2 import *

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

    def test_suffix_tree(self):
        text = 'ATAAATG$'
        suffix_tree = SuffixTree(text)
        edges = ['AAATG$', 'G$', 'T', 'ATG$', 'TG$', 'A', 'A', 'AAATG$', 'G$', 'T', 'G$', '$']
        assert Counter(suffix_tree.edges) == Counter(edges)

    def test_longest_repeat(self):
        text = 'ATATCGTTTTATCGTT'
        suffix_tree = SuffixTree(text)
        assert suffix_tree.longest_repeat() == 'TATCGTT'

    def test_longest_shared_substring(self):
        text_1 = 'TCGGTAGATTGCGCCCACTC'
        text_2 = 'AGGGGCTCGCAGTGTAAGAA'
        assert longest_shared_substring(text_1, text_2) == 'TCG' # AGA is another solution

class TestWeek2:
    def test_suffix_array(self):
        text = 'AACGATAGCGGTAGA$'
        expected = [15, 14, 0, 1, 12, 6, 4, 2, 8, 13, 3, 7, 9, 10, 11, 5]
        assert suffix_array(text) == expected

    def test_burrows_wheeler_transform(self):
        text = 'GCGTGCCTGGTCA$'
        bwt = 'ACTGGCT$TGCGGC'
        assert burrows_wheeler_transform(text) == bwt

    def test_invert_bwt(self):
        bwt = 'TTCCTAACG$A'
        text = 'TACATCACGT$'
        assert invert_bwt(bwt) == text