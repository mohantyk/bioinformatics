from week3 import *
from week4 import *

class TestWeek3:
    def test_affine_gap(self):
        v = 'KTGDRLFILHKCQHWEKSHNQWFHLNFRQWRPWWREILEVSPRNYWIISHCPHALWTWGGCFAEPRYWENGATCDSFSPA'
        w = 'KTGFPDVMLRTSYHKSHLQKLIFHMLKLIQCANFRQWRPWWREILSMGSPRNYWIISHCPHALNTWGGCFAEPRYAENGATCDSFSPA'
        score, align_v, align_w = align_with_affine_gap_penalty(v, w, 11, 1)

        assert score == 304
        assert align_v == 'KTGDRLFILHKCQHWEKSHNQW--FHL-------NFRQWRPWWREILEV-SPRNYWIISHCPHALWTWGGCFAEPRYWENGATCDSFSPA'
        assert align_w == 'KTGFPDVMLRTSYH--KSHLQKLIFHMLKLIQCANFRQWRPWWREILSMGSPRNYWIISHCPHALNTWGGCFAEPRYAENGATCDSFSPA'

    def test_middle_node(self):
        v = 'PLEASANTLY'
        w = 'MEASNLY'
        middle, score = get_middle_node(v, w, indel_penalty=5)
        assert middle == (4, 3)
        assert score == 17

    def test_middle_edge(self):
        v = 'PLEASANTLY'
        w = 'MEASNLY'
        middle_edge = get_middle_edge(v, w, indel_penalty=5)
        assert middle_edge[0] == (4, 3)
        assert middle_edge[1] == (5, 4)

    def test_linear_space_align(self):
        v = 'PLEASANTLY'
        w = 'MEANLY'
        path = linear_space_align(v, w)
        assert decode_path(v, w, path) == ('PLEASANTLY', '-MEA--N-LY')

    def test_multiple_lcs_score(self):
        words = ['ATATCCG','TCCGA','ATGTACTG']
        assert multiple_lcs_score(*words) == 3

    def test_multiple_lcs(self):
        words = ['ATATCCG','TCCGA','ATGTACTG']
        _, aligned = multiple_lcs(*words)
        assert score_alignment(*aligned) == 3


class TestWeek4:
    def test_greedy_sort(self):
        steps = greedy_sort_permutation([-3, 4, 1, 5, -2])
        expected = [[-1, -4, 3, 5, -2],
                    [1, -4, 3, 5, -2],
                    [1, +2, -5, -3, +4],
                    [1, 2, 3, 5, 4],
                    [1, 2, 3, -4, -5],
                    [1, 2, 3, 4, -5],
                    [1, 2, 3, 4, 5]]
        assert steps == expected
        
    def test_num_breakpoints(self):
        perm = [+3, +4, +5, -12, -8, -7, -6, +1, +2, +10, +9, -11, +13, +14]
        assert num_breakpoints(perm) == 8