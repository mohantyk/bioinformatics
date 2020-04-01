from collections import Counter

import logging
logging.basicConfig()
week4_logger = logging.getLogger('week4')
week4_logger.setLevel(logging.ERROR)

from week3 import spectrum, consistent
from week3 import MASS_2_AMINO

def spectrum_score(peptide_spectrum, ref_spectrum):
    p_count = Counter(peptide_spectrum)
    r_count = Counter(ref_spectrum)
    score = 0
    for mass in p_count:
        score += min(p_count[mass], r_count[mass])
    return score


def trim(peptides, ref_spectrum, N):
    if not peptides: # Empty set
        return []
    scores = {}
    for peptide in peptides:
        l_spectrum = spectrum(peptide, False)
        score = spectrum_score(l_spectrum, ref_spectrum)
        scores[peptide] = score
    ranked = sorted(scores, key=scores.get, reverse=True)

    trimmed = ranked[:N]
    lowest_rank = scores[trimmed[-1]]
    for peptide in ranked[N:]:
        if scores[peptide] < lowest_rank:
            break
        else:
            trimmed.append(peptide)
    return trimmed


def leaderboard_cyclopeptide_sequencing(mass_spectrum, N):
    leaderboard = set([''])
    leader_peptide = ''
    all_leaders = []
    parent_mass = mass_spectrum[-1]

    while leaderboard:
        leaderboard = set(candidate+aa for candidate in leaderboard
                                        for aa in MASS_2_AMINO.values())
        leader_spectrum = spectrum(leader_peptide, cyclic=True)
        candidates = leaderboard.copy()
        for candidate in candidates:
            c_spectrum = spectrum(candidate, cyclic=True)
            l_spectrum = spectrum(candidate, cyclic=False)
            peptide_mass = l_spectrum[-1]
            if peptide_mass == parent_mass:
                candidate_score = spectrum_score(c_spectrum, mass_spectrum)
                leader_score = spectrum_score(leader_spectrum, mass_spectrum)
                week4_logger.debug(f'{candidate} : {candidate_score} <> {leader_peptide} : {leader_score}')

                if candidate_score > leader_score:
                    leader_peptide = candidate
                    leader_spectrum = spectrum(leader_peptide, cyclic=True)
                    all_leaders.clear()
                    all_leaders.append(leader_peptide)
                elif candidate_score == leader_score:
                    all_leaders.append(candidate)

            elif peptide_mass > parent_mass:
                leaderboard.remove(candidate)
        leaderboard = trim(leaderboard, mass_spectrum, N)

    return all_leaders


def convolution(spectra):
    diffs = []
    for idx, mass1 in enumerate(spectra):
        for mass2 in spectra[:idx]:
            if mass2 == mass1: continue
            diffs.append(mass1-mass2)
    return diffs


def convolution_cyclopeptide_sequencing(mass_spectrum, M, N):
    ALLOWED_MASS_LOW = 57
    ALLOWED_MASS_HIGH = 200
    # Find alphabet by doing a spectrum convolution and choosing M most frequent masses
    conv_masses = convolution(mass_spectrum)
    allowed_masses = [mass for mass in conv_masses
                            if ALLOWED_MASS_LOW <= mass <= ALLOWED_MASS_HIGH]
    mass_counts = Counter(allowed_masses)
    ranked = sorted(mass_counts, key=mass_counts.get, reverse=True)

    alphabet = ranked[:M]
    lowest_count = mass_counts[alphabet[-1]]
    for mass in ranked[M:]:
        if mass_counts[mass] < lowest_count:
            break
        else:
            alphabet.append(mass)

    # Start sequencing algorithm
    leaderboard = set([''])
    leader_peptide = ''
    all_leaders = []
    parent_mass = mass_spectrum[-1]

    connector = '' # Connector for masses, empty at first, then hyphen
    while leaderboard:
        leaderboard = set(candidate+connector+str(mass) for candidate in leaderboard
                                                    for mass in alphabet)
        connector = '-' # Subsequent masses connected with a hyphen
        leader_spectrum = spectrum(leader_peptide, cyclic=True)
        candidates = leaderboard.copy()
        for candidate in candidates:
            c_spectrum = spectrum(candidate, cyclic=True)
            l_spectrum = spectrum(candidate, cyclic=False)
            peptide_mass = l_spectrum[-1]
            if peptide_mass == parent_mass:
                candidate_score = spectrum_score(c_spectrum, mass_spectrum)
                leader_score = spectrum_score(leader_spectrum, mass_spectrum)
                week4_logger.debug(f'{candidate} : {candidate_score} <> {leader_peptide} : {leader_score}')

                if candidate_score > leader_score:
                    leader_peptide = candidate
                    leader_spectrum = spectrum(leader_peptide, cyclic=True)
                    all_leaders.clear()
                    all_leaders.append(leader_peptide)
                elif candidate_score == leader_score:
                    all_leaders.append(candidate)

            elif peptide_mass > parent_mass:
                leaderboard.remove(candidate)
        leaderboard = trim(leaderboard, mass_spectrum, N)

    return all_leaders

