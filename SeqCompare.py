# This is to compare Amino Acid, DNA, or RNA Sequence Alignments, with flexible input formats
from Bio import Align
from Bio.Align import substitution_matrices
import re


def clean_amino_acid_sequence(sequence):
    """
    Clean a protein sequence by:
    1. Converting three-letter AA codes to one-letter
    2. Removing numbers and non-amino acid characters
    3. Standardizing whitespace and case

    Returns (cleaned_sequence, changes_made)
    """
    # Original sequence length for comparison
    original_length = len(sequence.replace(" ", "").replace("\n", ""))

    # Define three-letter to one-letter amino acid mapping
    aa_dict = {
        'ALA': 'A',  # Alanine
        'ASX': 'B',  # Asparagine or Aspartic acid
        'CYS': 'C',  # Cysteine
        'ASP': 'D',  # Aspartate
        'GLU': 'E',  # Glutamic acid
        'PHE': 'F',  # Phenylalanine
        'GLY': 'G',  # Glycine
        'HIS': 'H',  # Histidine
        'ILE': 'I',  # Isoleucine
        'LYS': 'K',  # Lysine
        'LEU': 'L',  # Leucine
        'MET': 'M',  # Methionine
        'ASN': 'N',  # Asparagine
        'PRO': 'P',  # Proline
        'GLN': 'Q',  # Glutamine
        'ARG': 'R',  # Arginine
        'SER': 'S',  # Serine
        'THR': 'T',  # Threonine
        'VAL': 'V',  # Valine
        'TRP': 'W',  # Tryptophan
        'TYR': 'Y',  # Tyrosine
        # Non-standard amino acids
        'SEC': 'U',  # Selenocysteine
        'PYL': 'O',  # Pyrrolysine
        'GLX': 'Z',  # Glutamine or Glutamic acid
        'XAA': 'X',  # Unknown amino acid
        'UNK': 'X',  # Unknown
        'XLE': 'J',  # Leucine or Isoleucine
        'TER': '*'  # Termination
    }

    # First check if this looks like a three-letter code sequence
    # Normalize whitespace first
    normalized = ' '.join(sequence.split())

    # Try to detect if this is a three-letter code sequence
    three_letter_pattern = r'\b(ALA|ARG|ASN|ASP|CYS|GLN|GLU|GLY|HIS|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL|SEC|PYL|ASX|GLX|XAA)\b'

    # Case-insensitive search for three-letter codes
    matches = re.findall(three_letter_pattern, normalized, re.IGNORECASE)

    # If we find significant matches, assume it's a three-letter code sequence
    if len(matches) > 5:  # Arbitrary threshold to detect three-letter format
        # Convert three-letter codes to one-letter
        result = ""
        words = normalized.split()
        for word in words:
            word_upper = word.upper()
            if word_upper in aa_dict:
                result += aa_dict[word_upper]
    else:
        # Assume it's a one-letter code sequence with possible contaminants
        # Keep only valid amino acid characters (case insensitive)
        valid_aa = "ACDEFGHIKLMNPQRSTVWYX*JUO"
        result = ''.join(c for c in sequence if c.upper() in valid_aa)

    # Return the cleaned sequence and whether significant changes were made
    cleaned_length = len(result)
    changes_made = original_length != cleaned_length
    return result.upper(), changes_made


def clean_dna_sequence(sequence):
    """
    Clean a DNA sequence by:
    1. Removing numbers, spaces, and non-DNA characters
    2. Standardizing case

    Returns (cleaned_sequence, changes_made)
    """
    original_length = len(sequence.replace(" ", "").replace("\n", ""))

    # Keep only valid DNA characters (A, T, G, C, N for unknown)
    valid_dna = "ATGCNRYKMSWBDHV"
    result = ''.join(c for c in sequence if c.upper() in valid_dna)

    cleaned_length = len(result)
    changes_made = original_length != cleaned_length
    return result.upper(), changes_made


def clean_rna_sequence(sequence):
    """
    Clean an RNA sequence by:
    1. Removing numbers, spaces, and non-RNA characters
    2. Standardizing case

    Returns (cleaned_sequence, changes_made)
    """
    original_length = len(sequence.replace(" ", "").replace("\n", ""))

    # Keep only valid RNA characters (A, U, G, C, N for unknown)
    valid_rna = "AUGCNRYKMSWBDHV"
    result = ''.join(c for c in sequence if c.upper() in valid_rna)

    cleaned_length = len(result)
    changes_made = original_length != cleaned_length
    return result.upper(), changes_made


def align_and_compare(seq1, seq2, seq_type, name1="Seq1", name2="Seq2"):
    """
    Align two sequences, compare them, and detail differences.

    Args:
        seq1 (str): First sequence
        seq2 (str): Second sequence
        seq_type (str): Sequence type: 'protein', 'dna', or 'rna'
        name1 (str): Name of first sequence (default: "Seq1")
        name2 (str): Name of second sequence (default: "Seq2")
    """
    # Initialize aligner with appropriate substitution matrix
    aligner = Align.PairwiseAligner()

    if seq_type == 'protein':
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -1
    elif seq_type in ['dna', 'rna']:
        # For nucleotide sequences, we can use simpler scoring
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -5
        aligner.extend_gap_score = -0.5

    # Perform alignment
    alignments = aligner.align(seq1, seq2)
    best_alignment = next(iter(alignments))

    # Display alignment
    print(f"\nAlignment between {name1} and {name2}:")
    print(best_alignment)
    print(f"Alignment score: {best_alignment.score:.1f}")

    # Compare original sequence lengths
    len_seq1, len_seq2 = len(seq1), len(seq2)
    length_diff = abs(len_seq1 - len_seq2)

    if seq_type == 'protein':
        unit = "residues"
    else:
        unit = "nucleotides"

    print(f"\nOriginal lengths: {name1} = {len_seq1}, {name2} = {len_seq2} {unit}")
    if len_seq2 > len_seq1:
        print(f"{name2} is longer than {name1} by {length_diff} {unit}.")
    elif len_seq1 > len_seq2:
        print(f"{name1} is longer than {name2} by {length_diff} {unit}.")
    else:
        print(f"Sequences are the same length ({len_seq1} {unit}).")

    # Extract aligned sequences and compare
    aligned_seq1, aligned_seq2 = best_alignment.sequences
    print("\nPosition-by-position comparison:")
    mismatches = 0
    gaps = 0
    for pos, (a, b) in enumerate(zip(aligned_seq1, aligned_seq2), 1):
        if a == b:
            symbol = "|"  # Match
        elif a == "-" or b == "-":
            symbol = " "  # Gap
            gaps += 1
        else:
            symbol = "*"  # Mismatch
            mismatches += 1
        print(f"Pos {pos:2d}: {a} {symbol} {b}")

    # Check for additional sequence in seq2
    if len_seq2 > len_seq1:
        extra_start = best_alignment.coordinates[1][-1]  # End index of seq2 in alignment
        if extra_start < len_seq2:
            extra_seq = seq2[extra_start:]
            print(f"\nAdditional sequence in {name2} beyond {name1}:")
            print(f"{unit.capitalize()} {extra_start + 1}-{len_seq2}: {extra_seq}")

    # Summary statistics
    length = len(aligned_seq1)
    percent_identity = (length - (mismatches + gaps)) / length * 100
    print(f"\nSummary:")
    print(f"Total differences:")
    print(f"  Mismatches: {mismatches}")
    print(f"  Gaps: {gaps}")
    print(f"  Length difference: {length_diff}")
    print(f"Alignment length: {length} {unit}")
    print(f"Percent identity: {percent_identity:.2f}%")


def get_sequence_input(seq_number, seq_type):
    """Get multi-line sequence input from user."""
    print(f"\nEnter Sequence {seq_number} (end with a single '.' on a new line):")

    if seq_type == 'protein':
        print("(Can be one-letter codes, three-letter codes, or include numbers/spaces)")
    elif seq_type == 'dna':
        print("(DNA sequence using A, T, G, C nucleotides)")
    elif seq_type == 'rna':
        print("(RNA sequence using A, U, G, C nucleotides)")

    lines = []
    while True:
        line = input()

        if line == '.':
            break
        if line.lower() == 'quit':
            return None

        lines.append(line)

    # Join the lines but don't remove spaces yet
    return ''.join(lines)


def main():
    print("Universal Sequence Alignment and Comparison Tool")
    print("=" * 70)
    print("This tool aligns and compares biological sequences:")
    print("• Protein (amino acid sequences)")
    print("• DNA (nucleotide sequences with A, T, G, C)")
    print("• RNA (nucleotide sequences with A, U, G, C)")
    print("\nFeatures:")
    print("• Handles various input formats (cleans and standardizes sequences)")
    print("• Uses appropriate scoring matrices for each sequence type")
    print("• Provides detailed alignment statistics and comparisons")
    print("\nInstructions:")
    print("1. Select the sequence type you want to compare")
    print("2. Enter each sequence (can be multiple lines)")
    print("3. Type a single period '.' on a new line to finish each sequence")
    print("4. Type 'quit' at any prompt to exit the program")
    print("=" * 70)

    while True:
        # Select sequence type
        print("\nSelect sequence type to compare:")
        print("1. Protein (amino acid sequences)")
        print("2. DNA (nucleotide sequences with A, T, G, C)")
        print("3. RNA (nucleotide sequences with A, U, G, C)")
        print("Type 'quit' to exit")

        seq_choice = input("> ").lower().strip()
        if seq_choice == 'quit':
            print("Exiting program.")
            break

        # Map choice to sequence type
        if seq_choice in ['1', 'protein', 'amino', 'aa', 'p']:
            seq_type = 'protein'
        elif seq_choice in ['2', 'dna', 'd']:
            seq_type = 'dna'
        elif seq_choice in ['3', 'rna', 'r']:
            seq_type = 'rna'
        else:
            print("Invalid choice. Please select 1, 2, 3 or type the sequence type.")
            continue

        print(f"\nYou selected: {seq_type.upper()} sequence comparison")

        # Get sequence names
        print("\nEnter a name for Sequence 1 (or press Enter for 'Seq1'):")
        name1 = input("> ").strip()
        if name1.lower() == 'quit':
            print("Exiting program.")
            break
        if not name1:
            name1 = "Seq1"

        print("\nEnter a name for Sequence 2 (or press Enter for 'Seq2'):")
        name2 = input("> ").strip()
        if name2.lower() == 'quit':
            print("Exiting program.")
            break
        if not name2:
            name2 = "Seq2"

        # Get sequence 1
        raw_seq1 = get_sequence_input(1, seq_type)
        if raw_seq1 is None:
            print("Exiting program.")
            break

        if not raw_seq1:
            print("No sequence provided. Please try again.")
            continue

        # Clean sequence 1 based on type
        if seq_type == 'protein':
            seq1, changes1 = clean_amino_acid_sequence(raw_seq1)
        elif seq_type == 'dna':
            seq1, changes1 = clean_dna_sequence(raw_seq1)
        else:  # RNA
            seq1, changes1 = clean_rna_sequence(raw_seq1)

        if changes1:
            print(f"\nSequence 1 was cleaned and processed.")
            print(f"Final sequence length: {len(seq1)} characters")

        if not seq1:
            print(f"No valid {seq_type} characters found in sequence 1. Please try again.")
            continue

        # Get sequence 2
        raw_seq2 = get_sequence_input(2, seq_type)
        if raw_seq2 is None:
            print("Exiting program.")
            break

        if not raw_seq2:
            print("No sequence provided. Please try again.")
            continue

        # Clean sequence 2 based on type
        if seq_type == 'protein':
            seq2, changes2 = clean_amino_acid_sequence(raw_seq2)
        elif seq_type == 'dna':
            seq2, changes2 = clean_dna_sequence(raw_seq2)
        else:  # RNA
            seq2, changes2 = clean_rna_sequence(raw_seq2)

        if changes2:
            print(f"\nSequence 2 was cleaned and processed.")
            print(f"Final sequence length: {len(seq2)} characters")

        if not seq2:
            print(f"No valid {seq_type} characters found in sequence 2. Please try again.")
            continue

        # Run the alignment and comparison
        try:
            align_and_compare(seq1, seq2, seq_type, name1, name2)
        except Exception as e:
            print(f"Error during alignment: {e}")

        # Offer to save the cleaned sequences
        print("\nDo you want to save the cleaned sequences? (y/n):")
        save = input("> ").lower()
        if save.startswith('y'):
            try:
                filename = input("Enter a filename (will save as .txt): ")
                if not filename.endswith('.txt'):
                    filename += '.txt'
                with open(filename, 'w') as f:
                    f.write(f">{name1} ({seq_type})\n")
                    f.write('\n'.join([seq1[i:i + 60] for i in range(0, len(seq1), 60)]))
                    f.write(f"\n\n>{name2} ({seq_type})\n")
                    f.write('\n'.join([seq2[i:i + 60] for i in range(0, len(seq2), 60)]))
                print(f"Sequences saved to {filename}")
            except Exception as e:
                print(f"Error saving file: {e}")

        # Ask if user wants to do another comparison
        print("\nDo you want to compare another pair of sequences? (y/n):")
        again = input("> ").lower()
        if not again.startswith('y'):
            print("Exiting program.")
            break


if __name__ == "__main__":
    main()