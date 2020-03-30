from collections import Counter

from week3 import spectrum

def spectrum_score(peptide_spectrum, ref_spectrum):
    p_count = Counter(peptide_spectrum)
    r_count = Counter(ref_spectrum)
    score = 0
    for mass in p_count:
        score += min(p_count[mass], r_count[mass])
    return score

def trim(peptides, ref_spectrum, N):
    scores = Counter()
    for peptide in peptides:
        l_spectrum = spectrum(peptide, False)
        score = spectrum_score(l_spectrum, ref_spectrum)
        scores[peptide] = score
    trimmed = [x[0] for x in scores.most_common(N)] # Does not account for ties
    return trimmed
