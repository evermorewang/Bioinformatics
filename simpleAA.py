#!/usr/bin/env python3

import re
import sys
import argparse

AA_DICT = {
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
    'ASN': 'N',  # Asparagine
    'ASP': 'D',  # Aspartic acid
    # Non-standard amino acids
    'SEC': 'U',  # Selenocysteine
    'PYL': 'O',  # Pyrrolysine
    'GLX': 'Z',  # Glutamine or Glutamic acid
    'XAA': 'X',  # Unknown amino acid
    'UNK': 'X',  # Unknown
    'XLE': 'J',  # Leucine or Isoleucine
    'TER': '*'  # Termination
}


def convert_sequence(text):
    """Convert three-letter amino acid codes to one-letter codes"""
    # Normalize whitespace
    normal_text = ' '.join(text.split())

    # Process each word
    result = ""
    for word in normal_text.split():
        word_upper = word.upper()
        if word_upper in AA_DICT:
            result += AA_DICT[word_upper]

    return result


def main():
    print("Three-Letter to One-Letter Amino Acid Converter")
    print("-" * 50)
    print("Instructions:")
    print("1. Enter your amino acid sequence (can be multiple lines)")
    print("2. Type a single period '.' on a new line to process the sequence")
    print("3. Type 'quit' to exit the program")
    print("-" * 50)

    while True:
        print("\nEnter amino acid sequence (end with a single '.' on a new line):")

        # Collect multi-line input
        lines = []
        while True:
            line = input()

            # Check for end marker or quit command
            if line == '.':
                break
            if line.lower() == 'quit':
                print("Exiting program.")
                return

            # Add the line to our collection
            lines.append(line)

        # If no input was provided, prompt again
        if not lines:
            print("No input provided. Try again or type 'quit' to exit.")
            continue

        # Join the lines and process
        full_query = " ".join(lines)
        result = convert_sequence(full_query)

        # Display result
        if result:
            print("\nConverted sequence:")
            print(result)
            print(f"Length: {len(result)} amino acids")
        else:
            print("No valid amino acids found in the input.")


if __name__ == "__main__":
    main()