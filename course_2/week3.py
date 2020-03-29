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

CODON_2_AMINO = get_rna_table()


def rna_to_protein(rna):
    n = len(rna)
    peptide = []
    for idx in range(0, n, 3):
        codon = rna[idx:idx+3]
        amino_acid = CODON_2_AMINO[codon]
        if amino_acid:
            peptide.append(amino_acid)
    protein = ''.join(peptide)
    return protein

def dna_to_rna(dna):
    return dna.replace('T','U')

def peptide_to_dna(dna, peptide):
    k = 3*len(peptide)
    n = len(dna)
    candidates = []
    for i in range(n-k+1):
        kmer = dna[i:i+k]
        reverse_kmer = reverse_complement(kmer)
        if ((rna_to_protein(dna_to_rna(kmer)) == peptide) or
            (rna_to_protein(dna_to_rna(reverse_kmer)) == peptide)):
            candidates.append(kmer)
    return candidates

