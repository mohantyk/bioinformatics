from week3 import align_with_affine_gap_penalty, get_middle_node

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
        middle = get_middle_node(v, w, indel_penalty=5)
        assert middle == (4, 3)