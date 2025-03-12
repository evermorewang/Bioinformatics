#!/usr/bin/env python3

# Minimum U Calculator for mRNA Sequences
# This script calculates the theoretical minimum number of U nucleotides
# required in an mRNA sequence encoding a given peptide

import re

# Codon table with min U codons for each amino acid
codon_table = {
    'A': 'GCC',  # Ala: 0 U (GCU, GCC, GCA, GCG)
    'R': 'CGC',  # Arg: 0 U (CGU, CGC, CGA, CGG, AGA, AGG)
    'N': 'AAC',  # Asn: 0 U (AAU, AAC)
    'D': 'GAC',  # Asp: 0 U (GAU, GAC)
    'C': 'UGC',  # Cys: 1 U (UGU, UGC)
    'Q': 'CAG',  # Gln: 0 U (CAA, CAG)
    'E': 'GAG',  # Glu: 0 U (GAA, GAG)
    'G': 'GGC',  # Gly: 0 U (GGU, GGC, GGA, GGG)
    'H': 'CAC',  # His: 0 U (CAU, CAC)
    'I': 'AUC',  # Ile: 1 U (AUU, AUC, AUA)
    'L': 'CUC',  # Leu: 0 U (UUA, UUG, CUU, CUC, CUA, CUG)
    'K': 'AAG',  # Lys: 0 U (AAA, AAG)
    'M': 'AUG',  # Met: 1 U (AUG, start codon only)
    'F': 'UUC',  # Phe: 2 U (UUU, UUC)
    'P': 'CCC',  # Pro: 0 U (CCU, CCC, CCA, CCG)
    'S': 'AGC',  # Ser: 0 U (UCU, UCC, UCA, UCG, AGU, AGC)
    'T': 'ACC',  # Thr: 0 U (ACU, ACC, ACA, ACG)
    'W': 'UGG',  # Trp: 1 U (UGG only)
    'Y': 'UAC',  # Tyr: 1 U (UAU, UAC)
    'V': 'GUC',  # Val: 0 U (GUU, GUC, GUA, GUG)
    '*': 'UAG'  # Stop: 1 U (UAA, UAG, UGA)
}

# Three-letter to one-letter amino acid code mapping
aa_three_to_one = {
    'ALA': 'A',  # Alanine
    'ARG': 'R',  # Arginine
    'ASN': 'N',  # Asparagine
    'ASP': 'D',  # Aspartic acid
    'CYS': 'C',  # Cysteine
    'GLN': 'Q',  # Glutamine
    'GLU': 'E',  # Glutamic acid
    'GLY': 'G',  # Glycine
    'HIS': 'H',  # Histidine
    'ILE': 'I',  # Isoleucine
    'LEU': 'L',  # Leucine
    'LYS': 'K',  # Lysine
    'MET': 'M',  # Methionine
    'PHE': 'F',  # Phenylalanine
    'PRO': 'P',  # Proline
    'SER': 'S',  # Serine
    'THR': 'T',  # Threonine
    'TRP': 'W',  # Tryptophan
    'TYR': 'Y',  # Tyrosine
    'VAL': 'V',  # Valine
    'SEC': 'U',  # Selenocysteine
    'PYL': 'O',  # Pyrrolysine
    'ASX': 'B',  # Aspartic acid or Asparagine
    'GLX': 'Z',  # Glutamic acid or Glutamine
    'XLE': 'J',  # Leucine or Isoleucine
    'XAA': 'X',  # Unknown amino acid
    'TER': '*'  # Stop codon
}


def clean_peptide_sequence(sequence):
    """
    Clean and normalize a peptide sequence by:
    1. Converting three-letter AA codes to one-letter if detected
    2. Removing non-amino acid characters, numbers, etc.
    3. Standardizing to uppercase
    4. Handling FASTA format if detected

    Args:
        sequence: Raw peptide sequence string

    Returns:
        Cleaned peptide sequence in one-letter AA code
    """
    # Check if sequence is in FASTA format and extract just the sequence part
    if sequence.startswith('>'):
        lines = sequence.strip().split('\n')
        sequence = ''.join(lines[1:])

    # Check if this might be three-letter code format
    # Normalize spacing first (replace all whitespace with single spaces)
    normalized = ' '.join(sequence.split())

    # Try to detect if this is a three-letter code sequence
    three_letter_pattern = r'\b(ALA|ARG|ASN|ASP|CYS|GLN|GLU|GLY|HIS|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL|SEC|PYL|ASX|GLX|XAA|TER)\b'

    # Case-insensitive search for three-letter codes
    matches = re.findall(three_letter_pattern, normalized, re.IGNORECASE)

    # If we find significant matches, assume it's a three-letter code sequence
    if len(matches) > 3:  # Arbitrary threshold to detect three-letter format
        # Convert three-letter codes to one-letter
        result = ""
        words = normalized.split()
        for word in words:
            word_upper = word.upper()
            if word_upper in aa_three_to_one:
                result += aa_three_to_one[word_upper]
    else:
        # Assume it's a one-letter code sequence with possible contaminants
        # Keep only valid amino acid characters (case insensitive)
        valid_aa = "ACDEFGHIKLMNPQRSTVWYX*JUO"
        result = ''.join(c for c in sequence if c.upper() in valid_aa)

    return result.upper()


def count_u(sequence):
    """Count the number of 'U' nucleotides in a sequence."""
    return sequence.count('U')


def calc_min_u(aa_sequence):
    """
    Calculate the minimum number of U's needed to encode a protein sequence
    using codons with the lowest U content for each amino acid.

    Args:
        aa_sequence: A string of amino acid letters (like "MEDMPVDPDNEAYF*")

    Returns:
        tuple: (minimum U optimized mRNA sequence, count of U nucleotides)
    """
    min_u_seq = ''
    unknown_aa = set()

    for aa in aa_sequence:
        if aa in codon_table:
            min_u_seq += codon_table[aa]
        else:
            unknown_aa.add(aa)

    if unknown_aa:
        print(f"Warning: Unknown amino acid(s) found and skipped: {', '.join(unknown_aa)}")

    min_u_count = count_u(min_u_seq)
    return min_u_seq, min_u_count


def show_codon_u_content():
    """Display U content for each amino acid's minimum-U codon"""
    print("\nU content in minimum-U codons:")
    print("-" * 40)
    print("Amino Acid | Codon | U Count")
    print("-" * 40)

    for aa, codon in sorted(codon_table.items()):
        u_count = count_u(codon)
        print(f"{aa:^10} | {codon:^5} | {u_count:^7}")


def format_sequence(sequence, width=60):
    """Format a sequence with line breaks for readability."""
    return '\n'.join([sequence[i:i + width] for i in range(0, len(sequence), width)])


def get_sequence_input():
    """Get multi-line sequence input from user."""
    print("\nEnter your amino acid sequence (end with a single '.' on a new line):")
    print("(Can be one-letter codes, three-letter codes, or in FASTA format)")

    lines = []
    while True:
        line = input()

        if line == '.':
            break
        if line.lower() == 'quit':
            return None

        lines.append(line)

    # Join the lines
    return '\n'.join(lines)


def main():
    print("Minimum U Calculator for mRNA Sequences")
    print("=" * 60)
    print("This tool calculates the theoretical minimum number of uracil (U)")
    print("nucleotides required in an mRNA sequence encoding a given peptide.")
    print("\nFeatures:")
    print("• Uses codons with minimum U content for each amino acid")
    print("• Handles various input formats (one-letter, three-letter AA codes)")
    print("• Automatically cleans and normalizes input sequences")
    print("• Provides detailed U content analysis")
    print("\nInstructions:")
    print("1. Enter your amino acid sequence (can be multiple lines)")
    print("2. Type a single period '.' on a new line to finish")
    print("3. Type 'quit' at any prompt to exit the program")
    print("=" * 60)

    while True:
        # Get sequence input
        raw_sequence = get_sequence_input()
        if raw_sequence is None or raw_sequence.lower() == 'quit':
            print("Exiting program.")
            break

        if not raw_sequence:
            print("No sequence provided. Please try again.")
            continue

        # Clean the sequence
        peptide = clean_peptide_sequence(raw_sequence)

        if not peptide:
            print("No valid amino acids found in input. Please try again.")
            continue

        print(f"\nCleaned peptide sequence ({len(peptide)} amino acids):")
        print(format_sequence(peptide))

        # Calculate minimum U content
        min_u_mrna, u_count = calc_min_u(peptide)

        print(f"\nMinimum-U mRNA sequence ({len(min_u_mrna)} nucleotides):")
        print(format_sequence(min_u_mrna))
        print(f"\nTotal U count: {u_count}")
        print(f"Percentage U: {u_count / len(min_u_mrna) * 100:.2f}%")

        # Show amino acid-specific U content
        show_codon_u_content()

        # Offer to save the results
        print("\nDo you want to save the results? (y/n):")
        save = input("> ").lower().strip()
        if save.startswith('y'):
            try:
                filename = input("Enter a filename (will save as .txt): ")
                if not filename.endswith('.txt'):
                    filename += '.txt'

                with open(filename, 'w') as f:
                    f.write("Minimum U Calculator Results\n")
                    f.write("=" * 60 + "\n\n")
                    f.write(f"Original Peptide Sequence ({len(peptide)} aa):\n")
                    f.write(format_sequence(peptide) + "\n\n")
                    f.write(f"Minimum-U mRNA Sequence ({len(min_u_mrna)} nt):\n")
                    f.write(format_sequence(min_u_mrna) + "\n\n")
                    f.write(f"Total U Count: {u_count}\n")
                    f.write(f"Percentage U: {u_count / len(min_u_mrna) * 100:.2f}%\n\n")
                    f.write("U Content by Amino Acid:\n")
                    f.write("-" * 40 + "\n")
                    f.write("Amino Acid | Codon | U Count\n")
                    f.write("-" * 40 + "\n")

                    for aa, codon in sorted(codon_table.items()):
                        u_content = count_u(codon)
                        f.write(f"{aa:^10} | {codon:^5} | {u_content:^7}\n")

                print(f"Results saved to {filename}")
            except Exception as e:
                print(f"Error saving file: {e}")

        # Ask to analyze another sequence
        print("\nDo you want to analyze another sequence? (y/n):")
        again = input("> ").lower().strip()
        if not again.startswith('y'):
            print("Exiting program.")
            break


if __name__ == "__main__":
    main()