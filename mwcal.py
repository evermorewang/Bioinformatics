#!/usr/bin/env python3
"""
Amino Acid Molecular Weight Calculator

This script calculates the molecular weight of a protein from its amino acid sequence.
It accepts various input formats and cleans the input to a single-letter amino acid sequence.
"""

import re
import sys

# Amino acid molecular weights in Daltons (g/mol)
# These values include the weight of the peptide backbone
AMINO_ACID_WEIGHTS = {
    'A': 71.08,  # Alanine
    'R': 156.19, # Arginine
    'N': 114.10, # Asparagine
    'D': 115.09, # Aspartic acid
    'C': 103.14, # Cysteine
    'Q': 128.13, # Glutamine
    'E': 129.12, # Glutamic acid
    'G': 57.05,  # Glycine
    'H': 137.14, # Histidine
    'I': 113.16, # Isoleucine
    'L': 113.16, # Leucine
    'K': 128.17, # Lysine
    'M': 131.19, # Methionine
    'F': 147.18, # Phenylalanine
    'P': 97.12,  # Proline
    'S': 87.08,  # Serine
    'T': 101.11, # Threonine
    'W': 186.21, # Tryptophan
    'Y': 163.18, # Tyrosine
    'V': 99.13,  # Valine
}

# Three-letter to one-letter amino acid code mapping
THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
    'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

def clean_sequence(sequence):
    """
    Clean the input sequence to obtain a single-letter amino acid sequence.
    
    Args:
        sequence (str): The input amino acid sequence in any format
        
    Returns:
        str: Cleaned single-letter amino acid sequence
    """
    # Convert to uppercase
    sequence = sequence.upper()
    
    # Remove common non-sequence characters (spaces, numbers, etc.)
    sequence = re.sub(r'[^A-Z]', '', sequence)
    
    # Check if the sequence is in three-letter format
    three_letter_pattern = r'(ALA|ARG|ASN|ASP|CYS|GLN|GLU|GLY|HIS|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL)'
    if re.match(three_letter_pattern * 2, sequence):  # Check if at least 2 consecutive three-letter codes
        # Convert three-letter format to one-letter
        one_letter_seq = ''
        for i in range(0, len(sequence), 3):
            if i + 3 <= len(sequence):
                three_letter = sequence[i:i+3]
                if three_letter in THREE_TO_ONE:
                    one_letter_seq += THREE_TO_ONE[three_letter]
                else:
                    # Skip if not a valid three-letter code
                    pass
        return one_letter_seq
    
    # If already in one-letter format, check for valid amino acids
    cleaned_seq = ''
    for aa in sequence:
        if aa in AMINO_ACID_WEIGHTS:
            cleaned_seq += aa
    
    return cleaned_seq

def calculate_molecular_weight(sequence):
    """
    Calculate the molecular weight of a protein from its amino acid sequence.
    
    Args:
        sequence (str): Single-letter amino acid sequence
        
    Returns:
        float: Molecular weight in Daltons (g/mol)
    """
    weight = 0.0
    
    for aa in sequence:
        if aa in AMINO_ACID_WEIGHTS:
            weight += AMINO_ACID_WEIGHTS[aa]
    
    # Add water (H2O) molecule weight for the last amino acid
    # (peptide bond formation releases water)
    water_weight = 18.02  # g/mol
    weight = weight + water_weight
    
    return weight

def main():
    print("Amino Acid Molecular Weight Calculator")
    print("--------------------------------------")
    print("Enter an amino acid sequence in any format.")
    print("The script will clean it and calculate the molecular weight.")
    print("Type 'exit' to quit the program.")
    print()
    
    while True:
        user_input = input("Enter amino acid sequence: ")
        
        if user_input.lower() == 'exit':
            print("Exiting the program.")
            sys.exit(0)
        
        if not user_input.strip():
            print("Please enter a valid sequence.")
            continue
        
        cleaned_sequence = clean_sequence(user_input)
        
        if not cleaned_sequence:
            print("No valid amino acids found in the input.")
            continue
        
        molecular_weight = calculate_molecular_weight(cleaned_sequence)
        
        print("\nResults:")
        print(f"Cleaned sequence: {cleaned_sequence}")
        print(f"Sequence length: {len(cleaned_sequence)} amino acids")
        print(f"Molecular weight: {molecular_weight:.2f} Da (g/mol)")
        print()

if __name__ == "__main__":
    main()
