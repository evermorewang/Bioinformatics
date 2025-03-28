from Bio.SeqUtils.ProtParam import ProteinAnalysis
from collections import Counter


def calculate_peptide_charge(sequence, pH=7.0):
    """
    Calculate the net charge of a peptide sequence at a given pH.

    Args:
        sequence (str): Amino acid sequence in one-letter code
        pH (float): pH value at which to calculate the charge (default: 7.0)

    Returns:
        float: Net charge of the peptide at the specified pH
    """
    # Clean the sequence (remove any non-standard AA characters)
    clean_seq = ''.join(aa for aa in sequence if aa.isalpha())

    # Create a ProteinAnalysis object
    protein = ProteinAnalysis(clean_seq)

    # Calculate the charge at the specified pH
    charge = protein.charge_at_pH(pH)

    return charge


def analyze_bonds(sequence):
    """
    Analyze and count different types of bonds in the peptide.

    Args:
        sequence (str): Amino acid sequence in one-letter code

    Returns:
        dict: Count of different bond types
    """
    # Amino acids that can form disulfide bonds
    disulfide_bonds = Counter(aa for aa in sequence if aa == 'C')
    potential_disulfide = disulfide_bonds['C'] // 2  # Each bond requires 2 cysteines

    # Amino acids that can form hydrogen bonds (donors)
    h_bond_donors = Counter(aa for aa in sequence if aa in 'STYNQKRH')

    # Amino acids that can form hydrogen bonds (acceptors)
    h_bond_acceptors = Counter(aa for aa in sequence if aa in 'DEQNYSTH')

    # Amino acids that can form ionic bonds
    positive_ions = Counter(aa for aa in sequence if aa in 'KRH')
    negative_ions = Counter(aa for aa in sequence if aa in 'DE')
    potential_ionic = min(sum(positive_ions.values()), sum(negative_ions.values()))

    # Amino acids that can form hydrophobic interactions
    hydrophobic_residues = Counter(aa for aa in sequence if aa in 'AILMFVPWG')

    return {
        "potential_disulfide_bonds": potential_disulfide,
        "h_bond_donors": dict(h_bond_donors),
        "h_bond_acceptors": dict(h_bond_acceptors),
        "positive_ions": dict(positive_ions),
        "negative_ions": dict(negative_ions),
        "potential_ionic_bonds": potential_ionic,
        "hydrophobic_residues": dict(hydrophobic_residues)
    }


def categorize_amino_acids(sequence):
    """
    Categorize amino acids as hydrophilic or hydrophobic.

    Args:
        sequence (str): Amino acid sequence in one-letter code

    Returns:
        dict: Counts and positions of hydrophilic and hydrophobic amino acids
    """
    # Define hydrophobic and hydrophilic amino acids
    hydrophobic = "AVILMFYWCG"
    hydrophilic = "RHKDESTNQP"

    # Count hydrophobic and hydrophilic amino acids
    hydrophobic_count = sum(1 for aa in sequence if aa in hydrophobic)
    hydrophilic_count = sum(1 for aa in sequence if aa in hydrophilic)

    # Get positions of hydrophobic and hydrophilic amino acids
    hydrophobic_positions = [i + 1 for i, aa in enumerate(sequence) if aa in hydrophobic]
    hydrophilic_positions = [i + 1 for i, aa in enumerate(sequence) if aa in hydrophilic]

    # Calculate percentage
    total_aa = len(sequence)
    hydrophobic_percentage = (hydrophobic_count / total_aa) * 100 if total_aa > 0 else 0
    hydrophilic_percentage = (hydrophilic_count / total_aa) * 100 if total_aa > 0 else 0

    return {
        "hydrophobic": {
            "count": hydrophobic_count,
            "percentage": hydrophobic_percentage,
            "amino_acids": "".join(sorted(set(aa for aa in sequence if aa in hydrophobic))),
            "positions": hydrophobic_positions
        },
        "hydrophilic": {
            "count": hydrophilic_count,
            "percentage": hydrophilic_percentage,
            "amino_acids": "".join(sorted(set(aa for aa in sequence if aa in hydrophilic))),
            "positions": hydrophilic_positions
        }
    }


def highlight_sequence(sequence):
    """
    Create a colored representation of the sequence highlighting
    hydrophobic and hydrophilic amino acids.

    Args:
        sequence (str): Amino acid sequence in one-letter code

    Returns:
        str: Formatted string with hydrophobic and hydrophilic highlighting indicators
    """
    result = ""
    hydrophobic = "AVILMFYWCG"
    hydrophilic = "RHKDESTNQP"

    for aa in sequence:
        if aa in hydrophobic:
            result += f"[H]{aa}[/H]"  # Hydrophobic indicator
        elif aa in hydrophilic:
            result += f"[P]{aa}[/P]"  # Hydrophilic (polar) indicator
        else:
            result += aa

    return result


def main():
    # Get peptide sequence from user
    print("Peptide Analysis Tool")
    print("--------------------")
    sequence = input("Enter the peptide sequence (one-letter amino acid code): ").strip().upper()

    if not sequence:
        print("Error: No sequence provided.")
        return

    # Ask for pH (with default option)
    try:
        pH_input = input("Enter pH value [7.0]: ").strip()
        pH = 7.0 if pH_input == "" else float(pH_input)
    except ValueError:
        print("Invalid pH value. Using default pH 7.0")
        pH = 7.0

    # Calculate charge at specified pH
    charge = calculate_peptide_charge(sequence, pH)
    print(f"\nResults:")
    print(f"Sequence: {sequence}")
    print(f"Length: {len(sequence)} amino acids")
    print(f"Net charge at pH {pH:.1f}: {charge:.2f}")

    # Analyze bonds
    print("\nBond Analysis:")
    print("-------------")
    bonds = analyze_bonds(sequence)

    print(f"Potential disulfide bonds: {bonds['potential_disulfide_bonds']}")
    print(f"Potential ionic bonds: {bonds['potential_ionic_bonds']}")

    print("\nHydrogen bond donors:")
    for aa, count in bonds['h_bond_donors'].items():
        print(f"  {aa}: {count}")

    print("\nHydrogen bond acceptors:")
    for aa, count in bonds['h_bond_acceptors'].items():
        print(f"  {aa}: {count}")

    print("\nPositively charged residues:")
    for aa, count in bonds['positive_ions'].items():
        print(f"  {aa}: {count}")

    print("\nNegatively charged residues:")
    for aa, count in bonds['negative_ions'].items():
        print(f"  {aa}: {count}")

    # Categorize amino acids
    print("\nAmino Acid Categorization:")
    print("------------------------")
    categories = categorize_amino_acids(sequence)

    print(f"Hydrophobic amino acids: {categories['hydrophobic']['count']} "
          f"({categories['hydrophobic']['percentage']:.1f}%)")
    print(f"  Types: {categories['hydrophobic']['amino_acids']}")

    print(f"Hydrophilic amino acids: {categories['hydrophilic']['count']} "
          f"({categories['hydrophilic']['percentage']:.1f}%)")
    print(f"  Types: {categories['hydrophilic']['amino_acids']}")

    # Print highlighted sequence
    print("\nSequence with hydrophobicity highlighting:")
    print(f"[H] = hydrophobic, [P] = hydrophilic (polar)")
    print(highlight_sequence(sequence))

    # Ask if user wants to see charge at different pH values
    show_range = input("\nWould you like to see the charge at different pH values? (y/n): ").strip().lower()

    if show_range == 'y' or show_range == 'yes':
        pH_values = [2.0, 4.0, 6.0, 7.0, 8.0, 10.0, 12.0]
        print("\npH\tCharge")
        print("-----------------")
        for pH in pH_values:
            charge = calculate_peptide_charge(sequence, pH)
            print(f"{pH:.1f}\t{charge:.2f}")


if __name__ == "__main__":
    main()