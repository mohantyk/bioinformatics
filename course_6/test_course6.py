from collections import Counter
import pandas as pd

from week1 import *
from week2 import *
from week4 import *
from week5 import *
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

    def test_bw_matching(self):
        bwt = 'TCCTCTATGAGATCCTATTCTATGAAACCTTCA$GACCAAAATTCTCCGGC'
        patterns = ['CCT', 'CAC', 'GAG', 'CAG', 'ATC']
        num_matches = [2,1,1,0,1]
        assert match_patterns(bwt, patterns) == num_matches


class TestWeek4:
    def test_viterbi(self):
        emitted = 'xyxzzxyxyy'
        alphabet = 'xyz'
        states = 'AB'
        transitions = [ [0.641, 0.359],
                        [0.729, 0.271] ]
        emissions = [   [0.117,0.691,0.192],
                        [0.097,0.42,0.483] ]
        best_path = 'AAABBAAAAA'

        assert viterbi(emitted, alphabet, states, transitions, emissions) == best_path

    def test_outcome_likelihood(self):
        emitted = 'xzyyzzyzyy'
        alphabet = 'xyz'
        states = 'AB'
        transitions = [ [0.303, 0.697],
                        [0.831, 0.169] ]
        emissions = [   [0.533,0.065,0.402],
                        [0.342,0.334,0.324] ]
        likelihood = 1.1005510319694847e-06

        assert outcome_likelihood(emitted, alphabet, states, transitions, emissions) == likelihood


class TestWeek5:
    def test_profile_hmm(self):
        threshold = 0.289
        alphabet = 'ABCDE'
        multiple_alignment = ['EBA', 'E-D', 'EB-', 'EED', 'EBD', 'EBE', 'E-D', 'E-D']

        # Expected HMM
        nodes = ['S', 'I0', 'M1', 'D1', 'I1', 'M2', 'D2', 'I2', 'E']
        transitions = pd.DataFrame(0, columns=nodes, index=nodes, dtype=float)
        transitions.loc['S', 'M1'] = 1.0
        transitions.loc['M1', 'I1'] = 0.625
        transitions.loc['M1', 'M2'] = 0.375
        transitions.loc['I1', 'M2'] = 0.8
        transitions.loc['I1', 'D2'] = 0.2
        transitions.loc['M2', 'E'] = 1.0
        transitions.loc['D2', 'E'] = 1.0

        emissions = pd.DataFrame(0, index=nodes, columns=list(alphabet), dtype=float)
        emissions.loc['M1', 'E'] = 1.0
        emissions.loc['I1', 'B'] = 0.8
        emissions.loc['I1', 'E'] = 0.2
        emissions.loc['M2', 'A'] = 0.143
        emissions.loc['M2', 'D'] = 0.714
        emissions.loc['M2', 'E'] = 0.143

        calc_transitions, calc_emissions = create_profile_hmm(multiple_alignment, threshold, alphabet)
        assert transitions.equals(calc_transitions)
        assert emissions.equals(calc_emissions)