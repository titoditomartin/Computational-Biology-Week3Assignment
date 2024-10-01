codon_map = {
    'AUG': ('Methionine', 'M'), 'UUU': ('Phenylalanine', 'F'), 'UUC': ('Phenylalanine', 'F'),
    'UUA': ('Leucine', 'L'), 'UUG': ('Leucine', 'L'), 'UCU': ('Serine', 'S'), 'UCC': ('Serine', 'S'),
    'UCA': ('Serine', 'S'), 'UCG': ('Serine', 'S'), 'UAU': ('Tyrosine', 'Y'), 'UAC': ('Tyrosine', 'Y'),
    'UGU': ('Cysteine', 'C'), 'UGC': ('Cysteine', 'C'), 'UGG': ('Tryptophan', 'W'), 'CUU': ('Leucine', 'L'),
    'CUC': ('Leucine', 'L'), 'CUA': ('Leucine', 'L'), 'CUG': ('Leucine', 'L'), 'CCU': ('Proline', 'P'),
    'CCC': ('Proline', 'P'), 'CCA': ('Proline', 'P'), 'CCG': ('Proline', 'P'), 'CAU': ('Histidine', 'H'),
    'CAC': ('Histidine', 'H'), 'CAA': ('Glutamine', 'Q'), 'CAG': ('Glutamine', 'Q'), 'CGU': ('Arginine', 'R'),
    'CGC': ('Arginine', 'R'), 'CGA': ('Arginine', 'R'), 'CGG': ('Arginine', 'R'), 'AUU': ('Isoleucine', 'I'),
    'AUC': ('Isoleucine', 'I'), 'AUA': ('Isoleucine', 'I'), 'GUU': ('Valine', 'V'), 'GUC': ('Valine', 'V'),
    'GUA': ('Valine', 'V'), 'GUG': ('Valine', 'V'), 'GCU': ('Alanine', 'A'), 'GCC': ('Alanine', 'A'),
    'GCA': ('Alanine', 'A'), 'GCG': ('Alanine', 'A'), 'GAU': ('Aspartic Acid', 'D'), 'GAC': ('Aspartic Acid', 'D'),
    'GAA': ('Glutamic Acid', 'E'), 'GAG': ('Glutamic Acid', 'E'), 'GGU': ('Glycine', 'G'), 'GGC': ('Glycine', 'G'),
    'GGA': ('Glycine', 'G'), 'GGG': ('Glycine', 'G'), 'UAA': ('Stop', None), 'UAG': ('Stop', None),
    'UGA': ('Stop', None), 'ACU': ('Threonine', 'T'), 'ACC': ('Threonine', 'T'), 'ACA': ('Threonine', 'T'),
    'ACG': ('Threonine', 'T'), 'AAU': ('Asparagine', 'N'), 'AAC': ('Asparagine', 'N'), 'AAA': ('Lysine', 'K'),
    'AAG': ('Lysine', 'K'), 'AGU': ('Serine', 'S'), 'AGC': ('Serine', 'S'), 'AGA': ('Arginine', 'R'),
    'AGG': ('Arginine', 'R')
}

codon_possibilities = {
    'Stop': ['UAA', 'UAG', 'UGA'], 'M': ['AUG'], 'F': ['UUU', 'UUC'],
    'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'], 'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
    'Y': ['UAU', 'UAC'], 'C': ['UGU', 'UGC'], 'W': ['UGG'], 'P': ['CCU', 'CCC', 'CCA', 'CCG'],
    'H': ['CAU', 'CAC'], 'Q': ['CAA', 'CAG'], 'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'I': ['AUU', 'AUC', 'AUA'], 'V': ['GUU', 'GUC', 'GUA', 'GUG'], 'A': ['GCU', 'GCC', 'GCA', 'GCG'],
    'D': ['GAU', 'GAC'], 'E': ['GAA', 'GAG'], 'G': ['GGU', 'GGC', 'GGA', 'GGG'],
    'T': ['ACU', 'ACC', 'ACA', 'ACG'], 'N': ['AAU', 'AAC'], 'K': ['AAA', 'AAG']
}


def is_valid_dna_sequence(dna_seq):
    return len(dna_seq) % 3 == 0 and all(base in "ATCG" for base in dna_seq)


def generate_complementary_dna(dna_seq):
    return ''.join({'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}[base] for base in dna_seq)


def convert_to_mrna(dna_seq):
    complement_dna = generate_complementary_dna(dna_seq)
    return complement_dna.replace('T', 'U')


def translate_to_amino_acids(dna_seq):
    mrna_seq = convert_to_mrna(dna_seq)
    amino_acids = []
    for i in range(0, len(mrna_seq), 3):
        codon = codon_map.get(mrna_seq[i:i + 3])
        if codon and codon[0] != 'Stop':
            amino_acids.append(f"{codon[0]} ({codon[1]})")
        else:
            print("Stop Codon detected")
            break
    return " - ".join(amino_acids)


def display_translation_steps(dna_seq):
    comp_dna = generate_complementary_dna(dna_seq)
    print(f"Complement DNA: {comp_dna}")
    print(f"mRNA: {convert_to_mrna(dna_seq)}")
    print(f"Amino Acids: {translate_to_amino_acids(dna_seq)}")


def is_valid_amino_sequence(amino_seq):
    return all(codon in codon_possibilities for codon in amino_seq)


def count_rna_codons(rna_seq):
    codon_counts = {}
    for i in range(0, len(rna_seq), 3):
        codon = rna_seq[i:i + 3]
        codon_counts[codon] = codon_counts.get(codon, 0) + 1

    print(f"\nmRNA Sequence: {rna_seq}")

   
    codons_with_c = {}
    codons_without_c = {}

    for codon, count in codon_counts.items():
        if 'C' in codon:
            codons_with_c[codon] = count
        else:
            codons_without_c[codon] = count

    
    for codon, count in codons_without_c.items():
        print(f"{codon}: {count}")
    for codon, count in codons_with_c.items():
        print(f"{codon}: {count}")


def generate_rna_sequences(amino_seq, current_rna=""):
    if not amino_seq:
        count_rna_codons(current_rna)
    else:
        for codon in codon_possibilities[amino_seq[0]]:
            generate_rna_sequences(amino_seq[1:], current_rna + codon)


while True:
    user_choice = input("\nPlease type 1 or 2 : 1) DNA to Amino Acids, 2) Amino Acids to mRNA: ").strip()
    
    if user_choice == '1':
        dna_input = input("Input DNA sequence: ").strip()
        while not is_valid_dna_sequence(dna_input):
            print("Invalid DNA sequence. Ensure the length is a multiple of 3 and contains valid characters.")
            dna_input = input("Input DNA sequence: ").strip()

        print()
        display_translation_steps(dna_input)

    elif user_choice == '2':
        amino_input = input("Input Amino Acid sequence (Max 3 characters): ").strip().upper()
        while len(amino_input) > 3 or not is_valid_amino_sequence(amino_input):
            print("Invalid amino acid sequence. Ensure it contains valid codons and is at most 3 characters.")
            amino_input = input("Input Amino Acid sequence (Max 3 characters): ").strip().upper()

        generate_rna_sequences(amino_input)
    
    else:
        print("Invalid input. Please choose 1 or 2.")
