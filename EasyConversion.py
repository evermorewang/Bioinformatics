#!/usr/bin/env python3

"""
Biological Sequence Converter
-----------------------------
A versatile tool for converting between different biological sequence formats:
- DNA to RNA and vice versa
- DNA/RNA to Amino Acid (protein translation)
- Three-letter amino acid codes to one-letter codes and vice versa
"""

import re
import sys
import argparse

# Complete Amino Acid dictionary (three-letter to one-letter codes)
AA_DICT = {
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
    # Ambiguous amino acids
    'ASX': 'B',  # Asparagine or Aspartic acid
    'GLX': 'Z',  # Glutamine or Glutamic acid
    'XLE': 'J',  # Leucine or Isoleucine
    # Non-standard amino acids
    'SEC': 'U',  # Selenocysteine
    'PYL': 'O',  # Pyrrolysine
    # Other
    'XAA': 'X',  # Unknown amino acid
    'UNK': 'X',  # Unknown
    'TER': '*'  # Termination
}

# One-letter to three-letter mapping (reverse of AA_DICT)
ONE_TO_THREE = {value: key for key, value in AA_DICT.items()}

# Complete Codon Table (RNA codons to amino acids)
CODON_TABLE = {
    # Phenylalanine (F)
    'UUU': 'F', 'UUC': 'F',
    # Leucine (L)
    'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    # Isoleucine (I)
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I',
    # Methionine (M) - Start codon
    'AUG': 'M',
    # Valine (V)
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    # Serine (S)
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'AGU': 'S', 'AGC': 'S',
    # Proline (P)
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    # Threonine (T)
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    # Alanine (A)
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    # Tyrosine (Y)
    'UAU': 'Y', 'UAC': 'Y',
    # Stop codons (*)
    'UAA': '*', 'UAG': '*', 'UGA': '*',
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
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
    # Glycine (G)
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

# DNA codon table (just replace U with T)
DNA_CODON_TABLE = {k.replace('U', 'T'): v for k, v in CODON_TABLE.items()}


def clean_sequence(sequence, seq_type):
    """
    Clean a biological sequence based on its type.

    Args:
        sequence (str): The sequence to clean
        seq_type (str): 'dna', 'rna', or 'protein'

    Returns:
        tuple: (cleaned_sequence, original_char_count, cleaned_char_count)
    """
    # Store original character count (excluding whitespace)
    original_no_whitespace = re.sub(r'\s', '', sequence)
    original_count = len(original_no_whitespace)

    # Remove all whitespace and numbers
    cleaned = re.sub(r'[\s\d]', '', sequence)

    if seq_type == 'dna':
        # Keep only valid DNA characters
        valid_seq = ''.join(c for c in cleaned if c.upper() in "ATGCNRYKMSWBDHV").upper()
    elif seq_type == 'rna':
        # Keep only valid RNA characters
        valid_seq = ''.join(c for c in cleaned if c.upper() in "AUGCNRYKMSWBDHV").upper()
    elif seq_type == 'protein':
        # Keep only valid protein characters
        valid_aa = "ACDEFGHIKLMNPQRSTVWY*BJZOUX"
        valid_seq = ''.join(c for c in cleaned if c.upper() in valid_aa).upper()
    else:
        valid_seq = cleaned.upper()

    # Return the cleaned sequence and count information
    return valid_seq, original_count, len(valid_seq)


def convert_aa_three_to_one(text):
    """
    Convert three-letter amino acid codes to one-letter codes

    Returns:
        tuple: (converted_sequence, original_count, converted_count)
    """
    # Normalize whitespace
    normal_text = ' '.join(text.split())

    # Store original count excluding whitespace
    original_no_whitespace = re.sub(r'\s', '', normal_text)
    original_count = len(original_no_whitespace)

    # Try to detect if this is a three-letter code sequence
    three_letter_pattern = r'\b(' + '|'.join(AA_DICT.keys()) + r')\b'

    # Case-insensitive search for three-letter codes
    matches = re.findall(three_letter_pattern, normal_text, re.IGNORECASE)

    # If we don't find matches, assume it might be one-letter already
    if len(matches) < 3:
        # Just clean it as a protein sequence
        cleaned_seq, _, cleaned_count = clean_sequence(normal_text, 'protein')
        return cleaned_seq, original_count, cleaned_count

    # Process each word
    result = ""
    for word in normal_text.split():
        word_upper = word.upper()
        if word_upper in AA_DICT:
            result += AA_DICT[word_upper]

    return result, original_count, len(result)


def convert_aa_one_to_three(sequence):
    """
    Convert one-letter amino acid codes to three-letter codes

    Returns:
        tuple: (converted_sequence, original_count, converted_count)
    """
    # Clean the sequence first
    cleaned, original_count, cleaned_count = clean_sequence(sequence, 'protein')

    # Convert each letter
    result = []
    for aa in cleaned:
        if aa in ONE_TO_THREE:
            result.append(ONE_TO_THREE[aa])
        else:
            result.append('XAA')  # Unknown amino acid

    return ' '.join(result), original_count, len(result)


def dna_to_rna(sequence):
    """
    Convert DNA to RNA (T → U)

    Returns:
        tuple: (converted_sequence, original_count, converted_count)
    """
    cleaned, original_count, cleaned_count = clean_sequence(sequence, 'dna')
    return cleaned.replace('T', 'U'), original_count, cleaned_count


def rna_to_dna(sequence):
    """
    Convert RNA to DNA (U → T)

    Returns:
        tuple: (converted_sequence, original_count, converted_count)
    """
    cleaned, original_count, cleaned_count = clean_sequence(sequence, 'rna')
    return cleaned.replace('U', 'T'), original_count, cleaned_count


def translate_to_protein(sequence, seq_type='rna'):
    """
    Translate RNA or DNA to protein.

    Args:
        sequence (str): RNA or DNA sequence
        seq_type (str): 'rna' or 'dna'

    Returns:
        tuple: (protein_sequence, original_count, protein_count, translation_details)
    """
    # Convert DNA to RNA if necessary
    if seq_type == 'dna':
        rna_sequence, original_count, cleaned_count = dna_to_rna(sequence)
    else:
        rna_sequence, original_count, cleaned_count = clean_sequence(sequence, 'rna')

    # Find start codon
    start_pos = rna_sequence.find('AUG')
    if start_pos == -1:
        return "No start codon (AUG) found.", original_count, 0, {
            "valid_sequence": rna_sequence,
            "nucleotide_count": cleaned_count,
            "error": "No start codon found"
        }

    # Translate codon by codon
    result = ""
    translation_details = {
        "valid_sequence": rna_sequence,
        "nucleotide_count": cleaned_count,
        "start_position": start_pos,
        "codons": []
    }

    for i in range(start_pos, len(rna_sequence) - 2, 3):
        codon = rna_sequence[i:i + 3]

        # Skip if not a complete codon
        if len(codon) < 3:
            continue

        # Get amino acid
        aa = CODON_TABLE.get(codon, 'X')  # X for unknown
        result += aa

        # Store codon info
        translation_details["codons"].append({
            "position": i,
            "codon": codon,
            "amino_acid": aa
        })

        # Stop at stop codon
        if aa == '*':
            break

        # Add feature to count nucleotides/amino acids in a sequence
        elif choice_num == 10:  # Sequence count
            print("\nCount Nucleotides or Amino Acids")
            print("Select sequence type:")
            print("1. DNA")
            print("2. RNA")
            print("3. Protein")
            seq_type_choice = input("> ").strip()

            if seq_type_choice == '1':
                seq_type = 'dna'
                seq = get_sequence_input("\nEnter DNA Sequence to Count")
            elif seq_type_choice == '2':
                seq_type = 'rna'
                seq = get_sequence_input("\nEnter RNA Sequence to Count")
            elif seq_type_choice == '3':
                seq_type = 'protein'
                seq = get_sequence_input("\nEnter Protein Sequence to Count")
            else:
                print("Invalid choice. Please try again.")
                continue

            if seq is None:
                print("Exiting program.")
                break

            if not seq:
                print("No sequence provided. Please try again.")
                continue

            cleaned_seq, original_count, cleaned_count = clean_sequence(seq, seq_type)

            print("\nSequence Analysis:")
            print(f"• Original input length: {len(seq)} characters")
            print(f"• Characters excluding whitespace: {original_count}")
            print(f"• Valid {seq_type.upper()} characters: {cleaned_count}")

            if seq_type in ['dna', 'rna']:
                base_counts = {}
                for base in cleaned_seq:
                    base_counts[base] = base_counts.get(base, 0) + 1

                print("\nNucleotide composition:")
                for base, count in sorted(base_counts.items()):
                    percentage = (count / cleaned_count) * 100
                    print(f"• {base}: {count} ({percentage:.2f}%)")

                gc_content = 0
                if seq_type == 'dna':
                    gc_content = (base_counts.get('G', 0) + base_counts.get('C', 0)) / cleaned_count * 100
                else:  # RNA
                    gc_content = (base_counts.get('G', 0) + base_counts.get('C', 0)) / cleaned_count * 100

                print(f"\nGC content: {gc_content:.2f}%")

            elif seq_type == 'protein':
                aa_counts = {}
                for aa in cleaned_seq:
                    aa_counts[aa] = aa_counts.get(aa, 0) + 1

                print("\nAmino acid composition:")
                for aa, count in sorted(aa_counts.items()):
                    percentage = (count / cleaned_count) * 100
                    three_letter = ONE_TO_THREE.get(aa, 'XXX')
                    print(f"• {aa} ({three_letter}): {count} ({percentage:.2f}%)")

    return result, original_count, len(result), translation_details


def find_open_reading_frames(sequence, seq_type='dna'):
    """
    Find all possible open reading frames (ORFs) in a DNA or RNA sequence.

    Args:
        sequence (str): DNA or RNA sequence
        seq_type (str): 'dna' or 'rna'

    Returns:
        tuple: (list_of_orfs, original_count, cleaned_count, orf_details)
    """
    # Ensure we're working with RNA
    if seq_type == 'dna':
        rna_seq, original_count, cleaned_count = dna_to_rna(sequence)
    else:
        rna_seq, original_count, cleaned_count = clean_sequence(sequence, 'rna')

    orfs = []
    orf_details = {
        "nucleotide_count": cleaned_count,
        "frames": {}
    }

    # Look in all three reading frames
    for frame in range(3):
        frame_orfs = []
        frame_details = []

        # Find all start codons
        i = frame
        while i < len(rna_seq) - 2:
            codon = rna_seq[i:i + 3]

            # If we find a start codon
            if codon == 'AUG':
                # Translate from here
                protein = ""
                codon_details = []
                start_pos = i
                j = i

                while j < len(rna_seq) - 2:
                    codon = rna_seq[j:j + 3]
                    if len(codon) < 3:
                        break

                    aa = CODON_TABLE.get(codon, 'X')
                    protein += aa

                    codon_details.append({
                        "position": j,
                        "codon": codon,
                        "amino_acid": aa
                    })

                    # Stop at stop codon
                    if aa == '*':
                        orfs.append(protein)
                        frame_orfs.append(protein)
                        frame_details.append({
                            "start_position": start_pos,
                            "end_position": j + 2,
                            "length": len(protein),
                            "sequence": protein,
                            "codons": codon_details
                        })
                        break

                    j += 3

                # If we didn't find a stop codon, still add the protein
                if j >= len(rna_seq) - 2 and protein and '*' not in protein:
                    orfs.append(protein)
                    frame_orfs.append(protein)
                    frame_details.append({
                        "start_position": start_pos,
                        "end_position": j - 1,
                        "length": len(protein),
                        "sequence": protein,
                        "codons": codon_details,
                        "note": "No stop codon found"
                    })

            i += 3

        orf_details["frames"][f"frame_{frame}"] = {
            "orf_count": len(frame_orfs),
            "orfs": frame_details
        }

    return orfs, original_count, cleaned_count, orf_details


def reverse_complement(sequence, seq_type='dna'):
    """
    Generate the reverse complement of a DNA or RNA sequence.

    Args:
        sequence (str): DNA or RNA sequence
        seq_type (str): 'dna' or 'rna'

    Returns:
        tuple: (reverse_complement_sequence, original_count, rc_count)
    """
    # Clean the sequence
    if seq_type == 'dna':
        cleaned, original_count, cleaned_count = clean_sequence(sequence, 'dna')
        complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                          'N': 'N', 'R': 'Y', 'Y': 'R', 'K': 'M',
                          'M': 'K', 'S': 'S', 'W': 'W', 'B': 'V',
                          'D': 'H', 'H': 'D', 'V': 'B'}
    else:
        cleaned, original_count, cleaned_count = clean_sequence(sequence, 'rna')
        complement_map = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G',
                          'N': 'N', 'R': 'Y', 'Y': 'R', 'K': 'M',
                          'M': 'K', 'S': 'S', 'W': 'W', 'B': 'V',
                          'D': 'H', 'H': 'D', 'V': 'B'}

    # Get the complement and reverse it
    complement = ''.join(complement_map.get(base, base) for base in cleaned)
    result = complement[::-1]
    return result, original_count, len(result)


def get_sequence_input(prompt_message):
    """
    Get multi-line sequence input from user.

    Args:
        prompt_message (str): Message to display before input

    Returns:
        str: User input or None if user wants to quit
    """
    print(prompt_message)
    print("Enter your sequence (end with a single '.' on a new line, 'quit' to exit):")
    print("(Include whitespace as needed, it will be removed during processing)")

    lines = []
    while True:
        line = input()

        if line == '.':
            break
        if line.lower() == 'quit':
            return None

        lines.append(line)

    return ''.join(lines)


def main():
    print("Biological Sequence Converter")
    print("=" * 70)
    print("This tool converts between different biological sequences:")
    print("• DNA ↔ RNA")
    print("• DNA/RNA → Protein (translation)")
    print("• One-letter ↔ Three-letter amino acid codes")
    print("• Find open reading frames (ORFs)")
    print("• Generate reverse complements")
    print("\nFeatures:")
    print("• Accurate sequence processing with whitespace/invalid character removal")
    print("• Detailed sequence statistics including character counts")
    print("• Robust handling of various input formats")
    print("\nInstructions:")
    print("1. Choose a conversion type")
    print("2. Enter your sequence(s) as prompted")
    print("3. Type a single period '.' on a new line to end input")
    print("4. Type 'quit' at any prompt to exit the program")
    print("=" * 70)

    while True:
        print("\nChoose a conversion type:")
        print("1. DNA → RNA")
        print("2. RNA → DNA")
        print("3. DNA → Protein")
        print("4. RNA → Protein")
        print("5. Three-letter amino acid codes → One-letter")
        print("6. One-letter amino acid codes → Three-letter")
        print("7. Find all open reading frames (ORFs) in DNA/RNA")
        print("8. Generate reverse complement (DNA/RNA)")
        print("9. Display codon table")
        print("Type 'quit' to exit")

        choice = input("> ").lower().strip()
        if choice == 'quit':
            print("Exiting program.")
            break

        try:
            choice_num = int(choice)
        except ValueError:
            print("Invalid choice. Please enter a number between 1-9.")
            continue

        if choice_num == 1:  # DNA → RNA
            dna_seq = get_sequence_input("\nDNA → RNA Conversion")
            if dna_seq is None:
                print("Exiting program.")
                break

            if not dna_seq:
                print("No sequence provided. Please try again.")
                continue

            result, original_count, cleaned_count = dna_to_rna(dna_seq)
            print("\nSequence processing:")
            print(f"• Original input length: {len(dna_seq)} characters")
            print(f"• Characters excluding whitespace: {original_count}")
            print(f"• Valid DNA nucleotides: {cleaned_count}")

            print("\nConverted RNA sequence:")
            for i in range(0, len(result), 60):
                print(result[i:i + 60])
            print(f"Length: {len(result)} nucleotides")

        elif choice_num == 2:  # RNA → DNA
            rna_seq = get_sequence_input("\nRNA → DNA Conversion")
            if rna_seq is None:
                print("Exiting program.")
                break

            if not rna_seq:
                print("No sequence provided. Please try again.")
                continue

            result, original_count, cleaned_count = rna_to_dna(rna_seq)
            print("\nSequence processing:")
            print(f"• Original input length: {len(rna_seq)} characters")
            print(f"• Characters excluding whitespace: {original_count}")
            print(f"• Valid RNA nucleotides: {cleaned_count}")

            print("\nConverted DNA sequence:")
            for i in range(0, len(result), 60):
                print(result[i:i + 60])
            print(f"Length: {len(result)} nucleotides")

        elif choice_num == 3:  # DNA → Protein
            dna_seq = get_sequence_input("\nDNA → Protein Translation")
            if dna_seq is None:
                print("Exiting program.")
                break

            if not dna_seq:
                print("No sequence provided. Please try again.")
                continue

            result, original_count, protein_count, details = translate_to_protein(dna_seq, 'dna')
            print("\nSequence processing:")
            print(f"• Original input length: {len(dna_seq)} characters")
            print(f"• Characters excluding whitespace: {original_count}")
            print(f"• Valid DNA nucleotides: {len(details['valid_sequence'])}")

            if isinstance(result, str) and result.startswith("No start codon"):
                print(f"\n{result}")
                continue

            print("\nTranslated protein sequence (one-letter codes):")
            for i in range(0, len(result), 60):
                print(result[i:i + 60])
            print(f"Length: {len(result)} amino acids")

            if details.get("start_position") is not None:
                print(f"Translation started at position: {details['start_position'] + 1}")

        elif choice_num == 4:  # RNA → Protein
            rna_seq = get_sequence_input("\nRNA → Protein Translation")
            if rna_seq is None:
                print("Exiting program.")
                break

            if not rna_seq:
                print("No sequence provided. Please try again.")
                continue

            result, original_count, protein_count, details = translate_to_protein(rna_seq, 'rna')
            print("\nSequence processing:")
            print(f"• Original input length: {len(rna_seq)} characters")
            print(f"• Characters excluding whitespace: {original_count}")
            print(f"• Valid RNA nucleotides: {len(details['valid_sequence'])}")

            if isinstance(result, str) and result.startswith("No start codon"):
                print(f"\n{result}")
                continue

            print("\nTranslated protein sequence (one-letter codes):")
            for i in range(0, len(result), 60):
                print(result[i:i + 60])
            print(f"Length: {len(result)} amino acids")

            if details.get("start_position") is not None:
                print(f"Translation started at position: {details['start_position'] + 1}")

        elif choice_num == 5:  # Three-letter → One-letter AA
            aa_seq = get_sequence_input("\nThree-letter → One-letter Amino Acid Conversion")
            if aa_seq is None:
                print("Exiting program.")
                break

            if not aa_seq:
                print("No sequence provided. Please try again.")
                continue

            result, original_count, converted_count = convert_aa_three_to_one(aa_seq)
            print("\nSequence processing:")
            print(f"• Original input length: {len(aa_seq)} characters")
            print(f"• Characters excluding whitespace: {original_count}")
            print(f"• Converted amino acids: {converted_count}")

            print("\nConverted amino acid sequence (one-letter codes):")
            for i in range(0, len(result), 60):
                print(result[i:i + 60])
            print(f"Length: {len(result)} amino acids")

        elif choice_num == 6:  # One-letter → Three-letter AA
            aa_seq = get_sequence_input("\nOne-letter → Three-letter Amino Acid Conversion")
            if aa_seq is None:
                print("Exiting program.")
                break

            if not aa_seq:
                print("No sequence provided. Please try again.")
                continue

            result, original_count, converted_count = convert_aa_one_to_three(aa_seq)
            print("\nSequence processing:")
            print(f"• Original input length: {len(aa_seq)} characters")
            print(f"• Characters excluding whitespace: {original_count}")
            print(f"• Valid amino acids: {converted_count}")

            print("\nConverted amino acid sequence (three-letter codes):")
            # Print with line breaks for readability
            words = result.split()
            line_length = 0
            line = []
            for word in words:
                if line_length + len(word) + 1 > 60:  # +1 for space
                    print(' '.join(line))
                    line = [word]
                    line_length = len(word)
                else:
                    line.append(word)
                    line_length += len(word) + 1
            if line:
                print(' '.join(line))

            print(f"Length: {len(words)} amino acids")

        elif choice_num == 7:  # Find ORFs
            print("\nFind Open Reading Frames (ORFs)")
            print("Select sequence type:")
            print("1. DNA")
            print("2. RNA")
            seq_type_choice = input("> ").strip()

            if seq_type_choice == '1':
                seq_type = 'dna'
                seq = get_sequence_input("\nDNA Sequence for ORF Analysis")
            elif seq_type_choice == '2':
                seq_type = 'rna'
                seq = get_sequence_input("\nRNA Sequence for ORF Analysis")
            else:
                print("Invalid choice. Please try again.")
                continue

            if seq is None:
                print("Exiting program.")
                break

            if not seq:
                print("No sequence provided. Please try again.")
                continue

            orfs, original_count, cleaned_count, details = find_open_reading_frames(seq, seq_type)

            print("\nSequence processing:")
            print(f"• Original input length: {len(seq)} characters")
            print(f"• Characters excluding whitespace: {original_count}")
            print(f"• Valid nucleotides: {cleaned_count}")

            if not orfs:
                print("\nNo open reading frames found.")
            else:
                print(f"\nFound {len(orfs)} open reading frames:")
                for i, orf in enumerate(orfs, 1):
                    print(f"\nORF #{i}:")
                    # Print with line breaks for readability
                    for j in range(0, len(orf), 60):
                        print(orf[j:j + 60])
                    print(f"Length: {len(orf)} amino acids")

        elif choice_num == 8:  # Reverse complement
            print("\nGenerate Reverse Complement")
            print("Select sequence type:")
            print("1. DNA")
            print("2. RNA")
            seq_type_choice = input("> ").strip()

            if seq_type_choice == '1':
                seq_type = 'dna'
                seq = get_sequence_input("\nDNA Sequence for Reverse Complement")
            elif seq_type_choice == '2':
                seq_type = 'rna'
                seq = get_sequence_input("\nRNA Sequence for Reverse Complement")
            else:
                print("Invalid choice. Please try again.")
                continue

            if seq is None:
                print("Exiting program.")
                break

            if not seq:
                print("No sequence provided. Please try again.")
                continue

            result, original_count, rc_count = reverse_complement(seq, seq_type)

            print("\nSequence processing:")
            print(f"• Original input length: {len(seq)} characters")
            print(f"• Characters excluding whitespace: {original_count}")
            print(f"• Valid nucleotides: {rc_count}")

            print("\nReverse complement sequence:")
            for i in range(0, len(result), 60):
                print(result[i:i + 60])
            print(f"Length: {len(result)} nucleotides")

        elif choice_num == 9:  # Display codon table
            print("\nCodon Table (RNA):")
            print("=" * 50)
            print("Codon  Amino Acid  Three-Letter Code")
            print("-" * 50)

            # Sort codons by amino acid for better readability
            sorted_codons = sorted(CODON_TABLE.items(), key=lambda x: x[1])

            current_aa = None
            for codon, aa in sorted_codons:
                if aa != current_aa:
                    current_aa = aa
                    print()  # Add line break between amino acids

                three_letter = ONE_TO_THREE.get(aa, 'XXX')
                if aa == '*':
                    three_letter = 'STOP'

                print(f"{codon}    {aa}           {three_letter}")

            print("\nCodon Table (DNA):")
            print("=" * 50)
            print("Codon  Amino Acid  Three-Letter Code")
            print("-" * 50)

            # Sort codons by amino acid for better readability
            sorted_dna_codons = sorted(DNA_CODON_TABLE.items(), key=lambda x: x[1])

            current_aa = None
            for codon, aa in sorted_dna_codons:
                if aa != current_aa:
                    current_aa = aa
                    print()  # Add line break between amino acids

                three_letter = ONE_TO_THREE.get(aa, 'XXX')
                if aa == '*':
                    three_letter = 'STOP'

                print(f"{codon}    {aa}           {three_letter}")

        else:
            print("Invalid choice. Please enter a number between 1-9.")
            continue

        # Ask if user wants to do another conversion
        print("\nDo you want to perform another conversion? (y/n):")
        again = input("> ").lower()
        if not again.startswith('y'):
            print("Exiting program.")
            break


if __name__ == "__main__":
    main()