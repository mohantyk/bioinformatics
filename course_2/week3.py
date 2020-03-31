from collections import Counter

import sys
sys.path.append('..')

from course_1.week2 import reverse_complement

def get_rna_table():
    rna_amino_file = 'RNA_codon_table_1.txt'
    codon_to_amino = {}
    with open(rna_amino_file, 'r') as f:
        for line in f:
            try:
                codon, rna = line.split()
            except ValueError: # Stop codon
                codon = line.strip()
                rna = None
            codon_to_amino[codon] = rna
    return codon_to_amino

def get_amino_mass():
    mass_table_file = 'integer_mass_table.txt'
    mass_table = {}
    with open(mass_table_file, 'r') as f:
        for line in f:
            amino, mass = line.strip().split()
            mass_table[amino] = int(mass)
    return mass_table

CODON_2_AMINO = get_rna_table()
AMINO_MASS = get_amino_mass()
MASS_2_AMINO = {v:k for k,v in AMINO_MASS.items()}
AMINO_NAMES = {'Leu': 'L', 'Arg': 'R', 'Pro': 'P', 'Gln': 'Q', 'His': 'H',
                'Met': 'M', 'Ile': 'I', 'Ser': 'S', 'Thr': 'T', 'Lys': 'K', 'Asn': 'N',
                'Phe': 'F', 'Trp': 'W', 'Cys': 'C', 'Tyr': 'Y',
                'Val': 'V', 'Gly': 'G', 'Ala': 'A', 'Glu': 'E', 'Asp':'D'}
AMINO_SYMBOLS = {v: k for k, v in AMINO_NAMES.items()}

def name_to_int(peptide):
    return '-'.join(str(AMINO_MASS[aa]) for aa in peptide)

def name_to_str(masses):
    assert isinstance(masses, str)
    return ''.join(MASS_2_AMINO[int(mass)] for mass in masses.split('-'))

def rna_to_peptide(rna):
    n = len(rna)
    peptide = []
    for idx in range(0, n, 3):
        codon = rna[idx:idx+3]
        amino_acid = CODON_2_AMINO[codon]
        if amino_acid:
            peptide.append(amino_acid)
    final_peptide = ''.join(peptide)
    return final_peptide

def dna_to_rna(dna):
    return dna.replace('T','U')

def peptide_to_dna(dna, peptide):
    k = 3*len(peptide)
    n = len(dna)
    candidates = []
    for i in range(n-k+1):
        kmer = dna[i:i+k]
        reverse_kmer = reverse_complement(kmer)
        if ((rna_to_peptide(dna_to_rna(kmer)) == peptide) or
            (rna_to_peptide(dna_to_rna(reverse_kmer)) == peptide)):
            candidates.append(kmer)
    return candidates


def spectrum(peptide, cyclic=True):
    n = len(peptide)
    cumsum = [0]*(len(peptide)+1)
    for i, aa in enumerate(peptide):
        cumsum[i+1] = cumsum[i] + AMINO_MASS[aa]
    peptide_mass = cumsum[-1]

    spectrum = [0]
    for start in range(1, n+1):
        for stop in range(start, n+1): # stop is included in kmer
            kmer_mass = cumsum[stop] - cumsum[start-1]
            spectrum.append(kmer_mass)
            if cyclic and (start > 1 and stop < n):
                cyclic_kmer_mass = peptide_mass - kmer_mass
                spectrum.append(cyclic_kmer_mass)
    final_spectrum = sorted(spectrum)
    return final_spectrum


def consistent(peptide_spectrum, mass_spectrum):
    p_count = Counter(peptide_spectrum)
    m_count = Counter(mass_spectrum)
    for mass in peptide_spectrum:
        if p_count[mass] > m_count[mass]:
            return False
    return True

def cyclopeptide_sequencing(mass_spectrum):
    candidates = set([''])
    final_peptides = []
    while candidates:
        candidates = set(candidate+aa for candidate in candidates
                                    for aa in MASS_2_AMINO.values())
        new_candidates = candidates.copy() # Can not change candidates while iterating over it
        for candidate in new_candidates:
            c_spectrum = spectrum(candidate, cyclic=True)
            l_spectrum = spectrum(candidate, cyclic=False)
            if c_spectrum == mass_spectrum:
              final_peptides.append(candidate)
            elif not consistent(l_spectrum, mass_spectrum):
                candidates.remove(candidate)
    return final_peptides
