#!/usr/bin/env python3

# Simple Gene to Amino Acid Converter
# Translates DNA or RNA sequences directly to amino acid sequences

import re

# Complete DNA codon table
DNA_CODON_TABLE = {
    # Phenylalanine (F)
    'TTT': 'F', 'TTC': 'F',

    # Leucine (L)
    'TTA': 'L', 'TTG': 'L', 'CTT': 'L',
    'CTC': 'L', 'CTA': 'L', 'CTG': 'L',

    # Isoleucine (I)
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I',

    # Methionine/Start (M)
    'ATG': 'M',

    # Valine (V)
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',

    # Serine (S)
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'AGT': 'S', 'AGC': 'S',

    # Proline (P)
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',

    # Threonine (T)
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',

    # Alanine (A)
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',

    # Tyrosine (Y)
    'TAT': 'Y', 'TAC': 'Y',

    # Histidine (H)
    'CAT': 'H', 'CAC': 'H',

    # Glutamine (Q)
    'CAA': 'Q', 'CAG': 'Q',

    # Asparagine (N)
    'AAT': 'N', 'AAC': 'N',

    # Lysine (K)
    'AAA': 'K', 'AAG': 'K',

    # Aspartic Acid (D)
    'GAT': 'D', 'GAC': 'D',

    # Glutamic Acid (E)
    'GAA': 'E', 'GAG': 'E',

    # Cysteine (C)
    'TGT': 'C', 'TGC': 'C',

    # Tryptophan (W)
    'TGG': 'W',

    # Arginine (R)
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGA': 'R', 'AGG': 'R',

    # Glycine (G)
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',

    # Stop codons (*)
    'TAA': '*', 'TAG': '*', 'TGA': '*'
}

# Complete RNA codon table
RNA_CODON_TABLE = {
    # Phenylalanine (F)
    'UUU': 'F', 'UUC': 'F',

    # Leucine (L)
    'UUA': 'L', 'UUG': 'L', 'CUU': 'L',
    'CUC': 'L', 'CUA': 'L', 'CUG': 'L',

    # Isoleucine (I)
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I',

    # Methionine/Start (M)
    'AUG': 'M',

    # Valine (V)
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',

    # Serine (S)
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'AGU': 'S', 'AGC': 'S',

    # Proline (P)
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',

    # Threonine (T)
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',

    # Alanine (A)
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',

    # Tyrosine (Y)
    'UAU': 'Y', 'UAC': 'Y',

    # Histidine (H)
    'CAU': 'H', 'CAC': 'H',

    # Glutamine (Q)
    'CAA': 'Q', 'CAG': 'Q',

    # Asparagine (N)
    'AAU': 'N', 'AAC': 'N',

    # Lysine (K)
    'AAA': 'K', 'AAG': 'K',

    # Aspartic Acid (D)
    'GAU': 'D', 'GAC': 'D',

    # Glutamic Acid (E)
    'GAA': 'E', 'GAG': 'E',

    # Cysteine (C)
    'UGU': 'C', 'UGC': 'C',

    # Tryptophan (W)
    'UGG': 'W',

    # Arginine (R)
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGA': 'R', 'AGG': 'R',

    # Glycine (G)
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',

    # Stop codons (*)
    'UAA': '*', 'UAG': '*', 'UGA': '*'
}


def clean_dna_sequence(sequence):
    """Clean a DNA sequence by removing non-DNA characters."""
    # Keep only A, T, G, C (case insensitive)
    clean = re.sub(r'[^ATGCatgc]', '', sequence)
    return clean.upper()


def clean_rna_sequence(sequence):
    """Clean an RNA sequence by removing non-RNA characters."""
    # Keep only A, U, G, C (case insensitive)
    clean = re.sub(r'[^AUGCaugc]', '', sequence)
    return clean.upper()


def detect_sequence_type(sequence):
    """Automatically detect if a sequence is DNA or RNA."""
    # If it contains T and no U, it's DNA
    if 'T' in sequence.upper() and 'U' not in sequence.upper():
        return 'dna'
    # If it contains U and no T, it's RNA
    elif 'U' in sequence.upper() and 'T' not in sequence.upper():
        return 'rna'
    # If it has neither or both, make an educated guess
    else:
        t_count = sequence.upper().count('T')
        u_count = sequence.upper().count('U')

        if t_count > u_count:
            return 'dna'
        elif u_count > t_count:
            return 'rna'
        else:
            # Default to DNA if we can't determine
            return 'dna'


def translate_sequence(sequence, seq_type=None):
    """
    Directly translate a DNA or RNA sequence to amino acids.

    Args:
        sequence: DNA or RNA sequence
        seq_type: 'dna' or 'rna' (if None, auto-detect)

    Returns:
        Amino acid sequence
    """
    # Auto-detect sequence type if not specified
    if seq_type is None:
        seq_type = detect_sequence_type(sequence)

    # Clean the sequence
    if seq_type == 'dna':
        sequence = clean_dna_sequence(sequence)
        codon_table = DNA_CODON_TABLE
    else:  # RNA
        sequence = clean_rna_sequence(sequence)
        codon_table = RNA_CODON_TABLE

    # Translate codons to amino acids
    protein = ""
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i + 3]

        # Skip incomplete codons at the end
        if len(codon) < 3:
            break

        # Look up the amino acid for this codon
        aa = codon_table.get(codon, 'X')  # 'X' for unknown codons
        protein += aa

    return protein


def get_sequence_input():
    """Get multi-line sequence input from user."""
    print("\nEnter your nucleotide sequence (end with a single '.' on a new line):")
    print("(DNA or RNA sequence with A, T/U, G, C nucleotides)")

    lines = []
    while True:
        line = input()

        if line == '.':
            break
        if line.lower() == 'quit':
            return None

        lines.append(line)

    # Join the lines
    return ''.join(lines)


def format_sequence(sequence, width=60):
    """Format a sequence with line breaks for readability."""
    return '\n'.join([sequence[i:i + width] for i in range(0, len(sequence), width)])


def main():
    print("Simple Gene to Amino Acid Converter")
    print("-" * 50)
    print("This tool translates DNA or RNA sequences into amino acid sequences")
    print("\nInstructions:")
    print("1. Enter your DNA or RNA sequence (can be multiple lines)")
    print("2. Type a single period '.' on a new line to finish")
    print("3. Type 'quit' at any prompt to exit the program")
    print("-" * 50)

    while True:
        # Get sequence
        raw_sequence = get_sequence_input()
        if raw_sequence is None:
            print("Exiting program.")
            break

        if not raw_sequence:
            print("No sequence provided. Please try again.")
            continue

        # Detect sequence type
        seq_type = detect_sequence_type(raw_sequence)

        # Clean sequence
        if seq_type == 'dna':
            cleaned_seq = clean_dna_sequence(raw_sequence)
            print("\nDetected DNA sequence.")
        else:  # RNA
            cleaned_seq = clean_rna_sequence(raw_sequence)
            print("\nDetected RNA sequence.")

        # Display cleaned sequence info
        print(f"\nCleaned sequence ({len(cleaned_seq)} nucleotides):")
        print(format_sequence(cleaned_seq))

        # Translate the sequence
        protein = translate_sequence(cleaned_seq, seq_type)

        # Display the translation
        print(f"\nTranslated amino acid sequence ({len(protein)} amino acids):")
        print(format_sequence(protein))

        # Count stop codons
        stop_codons = protein.count('*')
        if stop_codons > 0:
            print(f"\nNote: Contains {stop_codons} stop codon(s) shown as '*'")

        # Offer to save the results
        print("\nDo you want to save the results? (y/n):")
        save = input("> ").lower().strip()
        if save.startswith('y'):
            try:
                filename = input("Enter a filename (will save as .txt): ")
                if not filename.endswith('.txt'):
                    filename += '.txt'

                with open(filename, 'w') as f:
                    f.write(f">Original {seq_type.upper()} Sequence ({len(cleaned_seq)} nt)\n")
                    f.write(format_sequence(cleaned_seq))
                    f.write(f"\n\n>Translated Amino Acid Sequence ({len(protein)} aa)\n")
                    f.write(format_sequence(protein))

                print(f"Results saved to {filename}")
            except Exception as e:
                print(f"Error saving file: {e}")

        # Ask if user wants to translate another sequence
        print("\nDo you want to translate another sequence? (y/n):")
        again = input("> ").lower().strip()
        if not again.startswith('y'):
            print("Exiting program.")
            break


if __name__ == "__main__":
    main()