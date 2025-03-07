from Bio.SeqUtils import seq1
from Bio.SeqUtils.ProtParam import ProteinAnalysis


def extract_cdr_regions(heavy_chain_sequence):
    """
    Extract CDR regions from an antibody heavy chain sequence.

    Note: These are approximate positions based on Kabat numbering.
    Actual CDR boundaries can vary slightly depending on the method used.

    Args:
    heavy_chain_sequence (str): Full heavy chain Fab sequence

    Returns:
    dict: Dictionary containing CDR1, CDR2, and CDR3 sequences
    """
    # Approximate CDR regions for heavy chain (Kabat numbering)
    # These positions are standard but can have slight variations
    cdr_regions = {
        'CDR1': heavy_chain_sequence[26:35],  # Typically around positions 26-35
        'CDR2': heavy_chain_sequence[52:56],  # Typically around positions 52-56
        'CDR3': heavy_chain_sequence[95:102]  # Typically around positions 95-102
    }

    return cdr_regions


# Example usage
def main():
    # Example heavy chain Fab sequence (hypothetical)
    heavy_chain_sequence = (
        "QVQLVESGGGLVKPGGSLRLSCAASGFTFSDYYMSWIRQAPGKGLEWVSYITYSGSTIYYADSVKGRFTISRDNAKSSLYLQMNSLRAEDTAVYYCARDRGTTMVPFDYWGQGTLVTVSS"
    )

    # Extract CDR regions
    cdrs = extract_cdr_regions(heavy_chain_sequence)

    # Print out CDR sequences
    print("CDR Sequences:")
    for cdr, sequence in cdrs.items():
        print(f"{cdr}: {sequence}")

        # Optional: Calculate basic properties of each CDR
        cdr_analysis = ProteinAnalysis(sequence)
        print(f"  Length: {len(sequence)}")
        print(f"  Molecular Weight: {cdr_analysis.molecular_weight():.2f}")
        print(f"  Isoelectric Point: {cdr_analysis.isoelectric_point():.2f}")
        print()


if __name__ == "__main__":
    main()