#!/usr/bin/env python3

# Codon Optimization Tool for mRNA Sequences
# This script optimizes mRNA sequences based on:
# 1. Highest codon usage frequency
# 2. Minimal U content with high translation probability
# Based on Homo sapiens codon usage data

import re
from collections import defaultdict

# Codon usage frequency data from Homo sapiens [gbpri]
# Format: {codon: [frequency_per_thousand, count]}
codon_usage = {
    # Phenylalanine (F)
    'UUU': [17.6, 714298], 'UUC': [20.3, 824692],
    # Leucine (L)
    'UUA': [7.7, 311881], 'UUG': [12.9, 525688], 'CUU': [13.2, 536515],
    'CUC': [19.6, 796638], 'CUA': [7.2, 290751], 'CUG': [39.6, 1611801],
    # Isoleucine (I)
    'AUU': [16.0, 650473], 'AUC': [20.8, 846466], 'AUA': [7.5, 304565],
    # Methionine (M)
    'AUG': [22.0, 896005],
    # Valine (V)
    'GUU': [11.0, 448607], 'GUC': [14.5, 588138], 'GUA': [7.1, 287712], 'GUG': [28.1, 1143534],
    # Serine (S)
    'UCU': [15.2, 618711], 'UCC': [17.7, 718892], 'UCA': [12.2, 496448], 'UCG': [4.4, 179419],
    'AGU': [12.1, 493429], 'AGC': [19.5, 791383],
    # Proline (P)
    'CCU': [17.5, 713233], 'CCC': [19.8, 804620], 'CCA': [16.9, 688038], 'CCG': [6.9, 281570],
    # Threonine (T)
    'ACU': [13.1, 533609], 'ACC': [18.9, 768147], 'ACA': [15.1, 614523], 'ACG': [6.1, 246105],
    # Alanine (A)
    'GCU': [18.4, 750096], 'GCC': [27.7, 1127679], 'GCA': [15.8, 643471], 'GCG': [7.4, 299495],
    # Tyrosine (Y)
    'UAU': [12.2, 495699], 'UAC': [15.3, 622407],
    # Histidine (H)
    'CAU': [10.9, 441711], 'CAC': [15.1, 613713],
    # Glutamine (Q)
    'CAA': [12.3, 501911], 'CAG': [34.2, 1391973],
    # Asparagine (N)
    'AAU': [17.0, 689701], 'AAC': [19.1, 776603],
    # Lysine (K)
    'AAA': [24.4, 993621], 'AAG': [31.9, 1295568],
    # Aspartic Acid (D)
    'GAU': [21.8, 885429], 'GAC': [25.1, 1020595],
    # Glutamic Acid (E)
    'GAA': [29.0, 1177632], 'GAG': [39.6, 1609975],
    # Cysteine (C)
    'UGU': [10.6, 430311], 'UGC': [12.6, 513028],
    # Tryptophan (W)
    'UGG': [13.2, 535595],
    # Arginine (R)
    'CGU': [4.5, 184609], 'CGC': [10.4, 423516], 'CGA': [6.2, 250760], 'CGG': [11.4, 464485],
    'AGA': [12.2, 494682], 'AGG': [12.0, 486463],
    # Glycine (G)
    'GGU': [10.8, 437126], 'GGC': [22.2, 903565], 'GGA': [16.5, 669873], 'GGG': [16.5, 669768],
    # Stop codons
    'UAA': [1.0, 40285], 'UAG': [0.8, 32109], 'UGA': [1.6, 63237]
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

# One-letter to codon mapping
aa_to_codons = {
    'A': ['GCU', 'GCC', 'GCA', 'GCG'],  # Alanine
    'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],  # Arginine
    'N': ['AAU', 'AAC'],  # Asparagine
    'D': ['GAU', 'GAC'],  # Aspartic acid
    'C': ['UGU', 'UGC'],  # Cysteine
    'Q': ['CAA', 'CAG'],  # Glutamine
    'E': ['GAA', 'GAG'],  # Glutamic acid
    'G': ['GGU', 'GGC', 'GGA', 'GGG'],  # Glycine
    'H': ['CAU', 'CAC'],  # Histidine
    'I': ['AUU', 'AUC', 'AUA'],  # Isoleucine
    'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],  # Leucine
    'K': ['AAA', 'AAG'],  # Lysine
    'M': ['AUG'],  # Methionine
    'F': ['UUU', 'UUC'],  # Phenylalanine
    'P': ['CCU', 'CCC', 'CCA', 'CCG'],  # Proline
    'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],  # Serine
    'T': ['ACU', 'ACC', 'ACA', 'ACG'],  # Threonine
    'W': ['UGG'],  # Tryptophan
    'Y': ['UAU', 'UAC'],  # Tyrosine
    'V': ['GUU', 'GUC', 'GUA', 'GUG'],  # Valine
    '*': ['UAA', 'UAG', 'UGA']  # Stop
}

# Precompute the highest frequency codon for each amino acid
highest_freq_codons = {}
for aa, codons in aa_to_codons.items():
    highest_freq = 0
    best_codon = None
    for codon in codons:
        if codon in codon_usage and codon_usage[codon][0] > highest_freq:
            highest_freq = codon_usage[codon][0]
            best_codon = codon
    if best_codon:
        highest_freq_codons[aa] = best_codon

# Precompute the minimum U content codons for each amino acid
# Strictly prioritize minimum U count, even if it means using less frequent codons
min_u_codons = {}
for aa, codons in aa_to_codons.items():
    min_u_count = float('inf')
    candidates = []

    # First, find minimum U count among all codons for this amino acid
    for codon in codons:
        if codon in codon_usage:
            u_count = codon.count('U')
            if u_count < min_u_count:
                min_u_count = u_count
                candidates = [codon]
            elif u_count == min_u_count:
                candidates.append(codon)

    # Among candidates with minimum U, choose the one with highest frequency
    # but we're ensuring that min_u_count is strictly prioritized
    if candidates:
        best_codon = max(candidates, key=lambda c: codon_usage[c][0] if c in codon_usage else 0)
        min_u_codons[aa] = best_codon


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


def optimize_for_frequency(aa_sequence):
    """
    Optimize an amino acid sequence to use codons with highest frequency in human cells.

    Args:
        aa_sequence: A string of amino acid letters (one-letter code)

    Returns:
        tuple: (optimized RNA sequence, U count, optimized DNA sequence, annotated sequence, overall probability)
    """
    rna_seq = ''
    annotated_seq = []
    unknown_aa = set()
    codon_probabilities = []

    for aa in aa_sequence:
        if aa in highest_freq_codons:
            codon = highest_freq_codons[aa]
            rna_seq += codon

            # Get frequency (convert from per thousand to percentage)
            freq = codon_usage[codon][0] / 10.0  # Convert to percentage
            freq_decimal = freq / 100.0  # Convert to decimal for probability calculation
            codon_probabilities.append(freq_decimal)
            annotated_seq.append(f"{codon}({freq:.2f}%)")
        else:
            unknown_aa.add(aa)

    if unknown_aa:
        print(f"Warning: Unknown amino acid(s) found and skipped: {', '.join(unknown_aa)}")

    # Calculate overall probability (multiply all individual probabilities)
    overall_probability = 1.0
    for prob in codon_probabilities:
        overall_probability *= prob

    u_count = count_u(rna_seq)
    dna_seq = rna_to_dna(rna_seq)
    annotated = " ".join(annotated_seq)

    return rna_seq, u_count, dna_seq, annotated, overall_probability


def optimize_for_min_u(aa_sequence):
    """
    Optimize an amino acid sequence to use codons with minimum U content
    as the absolute priority, then consider frequency as a secondary factor.

    Args:
        aa_sequence: A string of amino acid letters (one-letter code)

    Returns:
        tuple: (optimized RNA sequence, U count, optimized DNA sequence, annotated sequence, overall probability)
    """
    rna_seq = ''
    annotated_seq = []
    unknown_aa = set()
    codon_probabilities = []

    for aa in aa_sequence:
        if aa in min_u_codons:
            codon = min_u_codons[aa]
            rna_seq += codon

            # Get frequency (convert from per thousand to percentage)
            freq = codon_usage[codon][0] / 10.0  # Convert to percentage
            freq_decimal = freq / 100.0  # Convert to decimal for probability calculation
            codon_probabilities.append(freq_decimal)
            annotated_seq.append(f"{codon}({freq:.2f}%)")
        else:
            unknown_aa.add(aa)

    if unknown_aa:
        print(f"Warning: Unknown amino acid(s) found and skipped: {', '.join(unknown_aa)}")

    # Calculate overall probability (multiply all individual probabilities)
    overall_probability = 1.0
    for prob in codon_probabilities:
        overall_probability *= prob

    u_count = count_u(rna_seq)
    dna_seq = rna_to_dna(rna_seq)
    annotated = " ".join(annotated_seq)

    return rna_seq, u_count, dna_seq, annotated, overall_probability


def analyze_codon_usage(sequence, title):
    """Analyze and display the codon usage in a sequence"""
    codon_counts = defaultdict(int)
    total_codons = 0

    # Count the codons
    for i in range(0, len(sequence), 3):
        if i + 3 <= len(sequence):
            codon = sequence[i:i + 3]
            codon_counts[codon] += 1
            total_codons += 1

    if total_codons == 0:
        return

    # Display the codon usage
    print(f"\n{title} Codon Usage Analysis:")
    print("-" * 60)
    print("Codon  | Amino Acid | Count | % of Sequence | Human Freq (/1000)")
    print("-" * 60)

    total_u = 0

    for codon, count in sorted(codon_counts.items()):
        # Determine which amino acid this codon encodes
        aa = "?"
        for a, codons in aa_to_codons.items():
            if codon in codons:
                aa = a
                break

        percentage = count / total_codons * 100
        human_freq = codon_usage.get(codon, [0])[0]
        total_u += codon.count('U') * count

        print(f"{codon:6} | {aa:^10} | {count:5} | {percentage:6.2f}%      | {human_freq:6.1f}")

    print("-" * 60)
    print(f"Total codons: {total_codons}")
    print(f"Total U nucleotides: {total_u}")
    print(f"Percentage U: {total_u / (total_codons * 3) * 100:.2f}%")


def format_sequence(sequence, width=60):
    """Format a sequence with line breaks for readability."""
    return '\n'.join([sequence[i:i + width] for i in range(0, len(sequence), width)])


def format_annotated_sequence(sequence, width=80):
    """Format an annotated sequence with line breaks for readability."""
    words = sequence.split()
    lines = []
    current_line = []
    current_length = 0

    for word in words:
        if current_length + len(word) + 1 > width and current_line:  # +1 for the space
            lines.append(' '.join(current_line))
            current_line = [word]
            current_length = len(word)
        else:
            current_line.append(word)
            current_length += len(word) + (1 if current_line else 0)

    if current_line:
        lines.append(' '.join(current_line))

    return '\n'.join(lines)


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
    print("Codon Optimization Tool for Human Expression")
    print("=" * 60)
    print("This tool optimizes codon usage based on the Homo sapiens genome")
    print("to maximize translation efficiency and/or minimize U content.")
    print("\nFeatures:")
    print("• Optimizes for highest-frequency codons to maximize translation")
    print("• Optimizes for minimum U content while maintaining good translation")
    print("• Handles various input formats (one-letter, three-letter AA codes)")
    print("• Provides detailed codon usage analysis")
    print("• Generates corresponding DNA sequences")
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

        # Optimize for highest frequency
        freq_rna, freq_u_count, freq_dna, freq_annotated, freq_probability = optimize_for_frequency(peptide)

        print(f"\n1. HIGHEST FREQUENCY OPTIMIZATION:")
        print(f"This sequence uses codons most commonly found in human mRNA.")
        print(f"\nRNA sequence ({len(freq_rna)} nucleotides):")
        print(format_sequence(freq_rna))
        print(f"\nDNA sequence:")
        print(format_sequence(freq_dna))
        print(f"\nAnnotated sequence with usage frequencies:")
        print(format_annotated_sequence(freq_annotated))
        print(f"\nU count: {freq_u_count}")
        print(f"Percentage U: {freq_u_count / len(freq_rna) * 100:.2f}%")
        print(f"Overall translation probability: {freq_probability:.10e}")

        # Analyze codon usage in frequency-optimized sequence
        analyze_codon_usage(freq_rna, "Frequency-Optimized")

        # Optimize for minimum U
        min_u_rna, min_u_count, min_u_dna, min_u_annotated, min_u_probability = optimize_for_min_u(peptide)

        print(f"\n2. MINIMUM U OPTIMIZATION:")
        print(f"This sequence uses codons with the absolute minimum U content")
        print(f"and chooses the highest frequency codon only among equal U-count options.")
        print(f"\nRNA sequence ({len(min_u_rna)} nucleotides):")
        print(format_sequence(min_u_rna))
        print(f"\nDNA sequence:")
        print(format_sequence(min_u_dna))
        print(f"\nAnnotated sequence with usage frequencies:")
        print(format_annotated_sequence(min_u_annotated))
        print(f"\nU count: {min_u_count}")
        print(f"Percentage U: {min_u_count / len(min_u_rna) * 100:.2f}%")
        print(f"Overall translation probability: {min_u_probability:.10e}")

        # Analyze codon usage in min-U-optimized sequence
        analyze_codon_usage(min_u_rna, "Min-U-Optimized")

        # Compare the two optimization strategies
        print("\nCOMPARISON OF STRATEGIES:")
        print("-" * 60)
        print(f"Highest Frequency: {freq_u_count} U's ({freq_u_count / len(freq_rna) * 100:.2f}%)")
        print(f"Minimum U: {min_u_count} U's ({min_u_count / len(min_u_rna) * 100:.2f}%)")
        print(
            f"Difference: {freq_u_count - min_u_count} U's ({(freq_u_count - min_u_count) / len(freq_rna) * 100:.2f}%)")
        print(f"\nTranslation probability comparison:")
        print(f"Highest Frequency: {freq_probability:.10e}")
        print(f"Minimum U: {min_u_probability:.10e}")
        if min_u_probability > 0:
            print(f"Ratio: {freq_probability / min_u_probability:.2f}x")
        else:
            print("Ratio: Cannot calculate (minimum U probability is 0)")

        # Count how many codons differ between the two strategies
        different_codons = 0
        for i in range(0, len(freq_rna), 3):
            if i + 3 <= len(freq_rna) and i + 3 <= len(min_u_rna):
                if freq_rna[i:i + 3] != min_u_rna[i:i + 3]:
                    different_codons += 1

        total_codons = len(freq_rna) // 3
        if total_codons > 0:
            percent_different = (different_codons / total_codons) * 100
            print(
                f"\nDifferent codons between strategies: {different_codons}/{total_codons} ({percent_different:.1f}%)")

        # Offer to save the results
        print("\nDo you want to save the results? (y/n):")
        save = input("> ").lower().strip()
        if save.startswith('y'):
            try:
                filename = input("Enter a filename (will save as .txt): ")
                if not filename.endswith('.txt'):
                    filename += '.txt'

                with open(filename, 'w') as f:
                    f.write("Codon Optimization Results\n")
                    f.write("=" * 60 + "\n\n")

                    f.write(f"Original Peptide Sequence ({len(peptide)} aa):\n")
                    f.write(format_sequence(peptide) + "\n\n")

                    f.write("1. HIGHEST FREQUENCY OPTIMIZATION\n")
                    f.write("-" * 60 + "\n")
                    f.write(f"RNA sequence ({len(freq_rna)} nucleotides):\n")
                    f.write(format_sequence(freq_rna) + "\n\n")
                    f.write(f"DNA sequence:\n")
                    f.write(format_sequence(freq_dna) + "\n\n")
                    f.write(f"Annotated sequence with usage frequencies:\n")
                    f.write(format_annotated_sequence(freq_annotated) + "\n\n")
                    f.write(f"U count: {freq_u_count}\n")
                    f.write(f"Percentage U: {freq_u_count / len(freq_rna) * 100:.2f}%\n")
                    f.write(f"Overall translation probability: {freq_probability:.10e}\n\n")

                    f.write("2. MINIMUM U OPTIMIZATION\n")
                    f.write("-" * 60 + "\n")
                    f.write(f"RNA sequence ({len(min_u_rna)} nucleotides):\n")
                    f.write(format_sequence(min_u_rna) + "\n\n")
                    f.write(f"DNA sequence:\n")
                    f.write(format_sequence(min_u_dna) + "\n\n")
                    f.write(f"Annotated sequence with usage frequencies:\n")
                    f.write(format_annotated_sequence(min_u_annotated) + "\n\n")
                    f.write(f"U count: {min_u_count}\n")
                    f.write(f"Percentage U: {min_u_count / len(min_u_rna) * 100:.2f}%\n")
                    f.write(f"Overall translation probability: {min_u_probability:.10e}\n\n")

                    f.write("COMPARISON OF STRATEGIES:\n")
                    f.write("-" * 60 + "\n")
                    f.write(f"Highest Frequency: {freq_u_count} U's ({freq_u_count / len(freq_rna) * 100:.2f}%)\n")
                    f.write(f"Minimum U: {min_u_count} U's ({min_u_count / len(min_u_rna) * 100:.2f}%)\n")
                    f.write(
                        f"Difference: {freq_u_count - min_u_count} U's ({(freq_u_count - min_u_count) / len(freq_rna) * 100:.2f}%)\n\n")
                    f.write(f"Translation probability comparison:\n")
                    f.write(f"Highest Frequency: {freq_probability:.10e}\n")
                    f.write(f"Minimum U: {min_u_probability:.10e}\n")
                    if min_u_probability > 0:
                        f.write(f"Ratio: {freq_probability / min_u_probability:.2f}x\n\n")
                    else:
                        f.write("Ratio: Cannot calculate (minimum U probability is 0)\n\n")

                    # Count how many codons differ between the two strategies
                    different_codons = 0
                    for i in range(0, len(freq_rna), 3):
                        if i + 3 <= len(freq_rna) and i + 3 <= len(min_u_rna):
                            if freq_rna[i:i + 3] != min_u_rna[i:i + 3]:
                                different_codons += 1

                    total_codons = len(freq_rna) // 3
                    if total_codons > 0:
                        percent_different = (different_codons / total_codons) * 100
                        f.write(
                            f"Different codons between strategies: {different_codons}/{total_codons} ({percent_different:.1f}%)\n")

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