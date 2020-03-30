from collections import Counter

def spectrum_score(peptide_spectrum, ref_spectrum):
    p_count = Counter(peptide_spectrum)
    r_count = Counter(ref_spectrum)
    score = 0
    for mass in p_count:
        score += min(p_count[mass], r_count[mass])
    return score
