#!/usr/bin/env python3

# Minimum U Calculator for mRNA Sequences
# This script calculates the theoretical minimum number of U nucleotides
# required in an mRNA sequence encoding a given peptide
# Upgraded to include DNA sequence and alternative codon options

import re
from collections import defaultdict

# Codon table with min U codons for each amino acid and all alternatives
codon_table = {
    'A': {'min_u': 'GCC', 'all': ['GCU', 'GCC', 'GCA', 'GCG']},  # Ala: 0 U
    'R': {'min_u': 'CGC', 'all': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG']},  # Arg: 0 U
    'N': {'min_u': 'AAC', 'all': ['AAU', 'AAC']},  # Asn: 0 U
    'D': {'min_u': 'GAC', 'all': ['GAU', 'GAC']},  # Asp: 0 U
    'C': {'min_u': 'UGC', 'all': ['UGU', 'UGC']},  # Cys: 1 U
    'Q': {'min_u': 'CAG', 'all': ['CAA', 'CAG']},  # Gln: 0 U
    'E': {'min_u': 'GAG', 'all': ['GAA', 'GAG']},  # Glu: 0 U
    'G': {'min_u': 'GGC', 'all': ['GGU', 'GGC', 'GGA', 'GGG']},  # Gly: 0 U
    'H': {'min_u': 'CAC', 'all': ['CAU', 'CAC']},  # His: 0 U
    'I': {'min_u': 'AUC', 'all': ['AUU', 'AUC', 'AUA']},  # Ile: 1 U
    'L': {'min_u': 'CUC', 'all': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG']},  # Leu: 0 U (for CUC)
    'K': {'min_u': 'AAG', 'all': ['AAA', 'AAG']},  # Lys: 0 U
    'M': {'min_u': 'AUG', 'all': ['AUG']},  # Met: 1 U (start codon only)
    'F': {'min_u': 'UUC', 'all': ['UUU', 'UUC']},  # Phe: 2 U
    'P': {'min_u': 'CCC', 'all': ['CCU', 'CCC', 'CCA', 'CCG']},  # Pro: 0 U
    'S': {'min_u': 'AGC', 'all': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC']},  # Ser: 0 U (for AGC)
    'T': {'min_u': 'ACC', 'all': ['ACU', 'ACC', 'ACA', 'ACG']},  # Thr: 0 U
    'W': {'min_u': 'UGG', 'all': ['UGG']},  # Trp: 1 U (UGG only)
    'Y': {'min_u': 'UAC', 'all': ['UAU', 'UAC']},  # Tyr: 1 U
    'V': {'min_u': 'GUC', 'all': ['GUU', 'GUC', 'GUA', 'GUG']},  # Val: 0 U (for GUC)
    '*': {'min_u': 'UAG', 'all': ['UAA', 'UAG', 'UGA']}  # Stop: 1 U
}

# Group codons by U count
codons_by_u_count = defaultdict(lambda: defaultdict(list))
for aa, data in codon_table.items():
    for codon in data['all']:
        u_count = codon.count('U')
        codons_by_u_count[aa][u_count].append(codon)

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


def rna_to_dna(rna_sequence):
    """Convert an RNA sequence to DNA by replacing U with T."""
    return rna_sequence.replace('U', 'T')


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
            min_u_seq += codon_table[aa]['min_u']
        else:
            unknown_aa.add(aa)

    if unknown_aa:
        print(f"Warning: Unknown amino acid(s) found and skipped: {', '.join(unknown_aa)}")

    min_u_count = count_u(min_u_seq)
    return min_u_seq, min_u_count


def find_alternative_sequences(aa_sequence, max_u_diff=0):
    """
    Find all sequences with U counts close to the minimum.

    Args:
        aa_sequence: A string of amino acid letters
        max_u_diff: Maximum additional U's allowed compared to minimum

    Returns:
        list: Alternative RNA sequences that have U counts within max_u_diff
              of the minimum possible
    """
    min_u_seq, min_u_count = calc_min_u(aa_sequence)

    # Initialize with just the minimum sequence
    if max_u_diff == 0:
        return [min_u_seq]

    # Generate alternative sequences
    from itertools import product

    # For each amino acid, get the codons with U count at or near minimum
    codon_options = []
    for aa in aa_sequence:
        if aa not in codon_table:
            continue

        # Find minimum U count for this amino acid
        min_u_for_aa = min(codons_by_u_count[aa].keys())

        # Collect all codons within max_u_diff of minimum
        valid_codons = []
        for u_count in range(min_u_for_aa, min_u_for_aa + max_u_diff + 1):
            if u_count in codons_by_u_count[aa]:
                valid_codons.extend(codons_by_u_count[aa][u_count])

        codon_options.append(valid_codons)

    # Calculate total combinations (may be very large)
    total_combos = 1
    for options in codon_options:
        total_combos *= len(options)

    # If too many combinations, return just the minimum
    if total_combos > 100:  # Arbitrary limit to prevent excessive computation
        return [min_u_seq]

    # Generate all combinations within the U count constraint
    alternative_seqs = []
    for codon_combo in product(*codon_options):
        rna_seq = ''.join(codon_combo)
        u_count = count_u(rna_seq)

        if u_count <= min_u_count + max_u_diff:
            alternative_seqs.append(rna_seq)

    return alternative_seqs


def show_codon_u_content():
    """Display U content for each amino acid's minimum-U codon"""
    print("\nU content in minimum-U codons:")
    print("-" * 40)
    print("Amino Acid | Codon | U Count")
    print("-" * 40)

    for aa, data in sorted(codon_table.items()):
        codon = data['min_u']
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
    print("• Generates corresponding DNA sequence")
    print("• Finds alternative sequences with minimal U content")
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

        # Generate the DNA sequence
        dna_sequence = rna_to_dna(min_u_mrna)

        print(f"\nMinimum-U mRNA sequence ({len(min_u_mrna)} nucleotides):")
        print(format_sequence(min_u_mrna))
        print(f"\nCorresponding DNA sequence:")
        print(format_sequence(dna_sequence))
        print(f"\nTotal U count: {u_count}")
        print(f"Percentage U: {u_count / len(min_u_mrna) * 100:.2f}%")

        # Show amino acid-specific U content
        show_codon_u_content()

        # Ask if user wants alternative sequences
        print("\nDo you want to see alternative sequences with near-minimal U content? (y/n):")
        alt_seq = input("> ").lower().strip()
        if alt_seq.startswith('y'):
            max_u_diff = 1
            try:
                max_u_diff = int(input("Maximum additional U's allowed (1-3 recommended): "))
            except ValueError:
                print("Using default value of 1")
                max_u_diff = 1

            alt_sequences = find_alternative_sequences(peptide, max_u_diff)

            print(f"\nFound {len(alt_sequences)} alternative sequences:")
            for i, seq in enumerate(alt_sequences, 1):
                u_ct = count_u(seq)
                dna_seq = rna_to_dna(seq)

                print(f"\n---Alternative #{i}---")
                print(f"RNA: {format_sequence(seq)}")
                print(f"DNA: {format_sequence(dna_seq)}")
                print(f"U count: {u_ct} ({u_ct / len(seq) * 100:.2f}%)")

                # Limit display to first 5 sequences if there are many
                if i >= 5 and len(alt_sequences) > 5:
                    remain = len(alt_sequences) - 5
                    print(f"\n...and {remain} more sequences (not displayed)")
                    break

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
                    f.write(f"Corresponding DNA Sequence:\n")
                    f.write(format_sequence(dna_sequence) + "\n\n")
                    f.write(f"Total U Count: {u_count}\n")
                    f.write(f"Percentage U: {u_count / len(min_u_mrna) * 100:.2f}%\n\n")

                    f.write("U Content by Amino Acid:\n")
                    f.write("-" * 40 + "\n")
                    f.write("Amino Acid | Codon | U Count\n")
                    f.write("-" * 40 + "\n")

                    for aa, data in sorted(codon_table.items()):
                        codon = data['min_u']
                        u_content = count_u(codon)
                        f.write(f"{aa:^10} | {codon:^5} | {u_content:^7}\n")

                    # Add alternative sequences if they were generated
                    if alt_seq.startswith('y'):
                        f.write("\n\nAlternative Sequences:\n")
                        f.write("=" * 60 + "\n")

                        for i, seq in enumerate(alt_sequences, 1):
                            u_ct = count_u(seq)
                            dna_seq = rna_to_dna(seq)

                            f.write(f"\n---Alternative #{i}---\n")
                            f.write(f"RNA: {format_sequence(seq)}\n")
                            f.write(f"DNA: {format_sequence(dna_seq)}\n")
                            f.write(f"U count: {u_ct} ({u_ct / len(seq) * 100:.2f}%)\n")

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