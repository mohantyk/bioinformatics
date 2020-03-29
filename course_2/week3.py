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
