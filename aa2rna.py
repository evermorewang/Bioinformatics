from Bio.Seq import Seq
from Bio.Data import CodonTable
import re


def clean_amino_acid_input(input_sequence):
    """
    Clean and standardize amino acid input to single letter code.
    Handles both 3-letter and 1-letter amino acid codes.
    """
    # Dictionary for 3-letter to 1-letter conversion
    aa_dict = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
        'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
        'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
        'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
        'STOP': '*'
    }

    # Check if input is likely 3-letter code (contains letters other than single AA codes)
    if any(c not in "ACDEFGHIKLMNPQRSTVWY*" for c in input_sequence.upper()):
        # Split by spaces, dashes or commas if present
        parts = re.split(r'[\s,-]+', input_sequence)
        # Convert each part to 1-letter code
        cleaned_parts = []
        for part in parts:
            part = part.upper()
            if part in aa_dict:
                cleaned_parts.append(aa_dict[part])
            else:
                # Try to match partial names (e.g., "ala" to "ALA")
                matched = False
                for key in aa_dict:
                    if key.startswith(part) or part.startswith(key):
                        cleaned_parts.append(aa_dict[key])
                        matched = True
                        break
                if not matched:
                    # Keep as is if no match found
                    cleaned_parts.append(part[0] if len(part) > 0 else '')

        return ''.join(cleaned_parts)
    else:
        # Already in 1-letter code, just remove spaces and non-AA characters
        return ''.join(c for c in input_sequence.upper() if c in "ACDEFGHIKLMNPQRSTVWY*")


def amino_acid_to_mrna(aa_sequence):
    """
    Convert amino acid sequence to mRNA sequence.
    Uses the standard codon table and returns multiple possible mRNA sequences.
    """
    # Get the standard codon table
    standard_table = CodonTable.standard_dna_table

    # Build a dictionary of amino acids to codons
    aa_to_codon = {}
    for codon, aa in standard_table.forward_table.items():
        if aa not in aa_to_codon:
            aa_to_codon[aa] = []
        aa_to_codon[aa].append(codon)

    # Add stop codons
    aa_to_codon['*'] = standard_table.stop_codons

    # Generate the mRNA sequence (replace T with U)
    mrna_options = []
    for aa in aa_sequence:
        if aa in aa_to_codon:
            # Show all possible codons
            codons = aa_to_codon[aa]
            mrna_options.append(f"{aa}({', '.join(c.replace('T', 'U') for c in codons)})")
        else:
            mrna_options.append(f"{aa}(?)")

    return ' '.join(mrna_options)


def display_codon_table():
    """
    Display the standard codon table.
    """
    standard_table = CodonTable.standard_dna_table

    print("\nComplete Codon Table:")
    print("==================")

    # Format each codon and its amino acid
    codon_list = []
    for codon in sorted(standard_table.forward_table.keys()):
        aa = standard_table.forward_table[codon]
        codon_list.append(f"{codon.replace('T', 'U')}: {aa}")

    # Add stop codons
    for codon in standard_table.stop_codons:
        codon_list.append(f"{codon.replace('T', 'U')}: STOP(*)")

    # Print in columns
    columns = 4
    rows = (len(codon_list) + columns - 1) // columns

    for i in range(rows):
        row = []
        for j in range(columns):
            idx = i + j * rows
            if idx < len(codon_list):
                row.append(codon_list[idx].ljust(15))
            else:
                row.append("".ljust(15))
        print("".join(row))


def main():
    print("Amino Acid to mRNA Converter")
    print("===========================")
    print("Input can be 1-letter (ACGT) or 3-letter (Ala-Cys-Gly-Thr) amino acid codes.")
    print("Enter 'table' to display the codon table.")
    print("Enter 'quit' to exit.")

    while True:
        query = input("\nEnter amino acid sequence: ")

        if query.lower() == 'quit':
            break
        elif query.lower() == 'table':
            display_codon_table()
            continue

        # Clean the input
        cleaned_sequence = clean_amino_acid_input(query)

        if cleaned_sequence:
            print(f"\nCleaned amino acid sequence: {cleaned_sequence}")
            print(f"Possible mRNA codons:")
            print(amino_acid_to_mrna(cleaned_sequence))
        else:
            print("Invalid input. Please try again.")


if __name__ == "__main__":
    main()