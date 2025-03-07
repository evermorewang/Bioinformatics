from Bio import Align
from Bio.Align import substitution_matrices


def align_and_compare(seq1, seq2, name1="Seq1", name2="Seq2"):
    """
    Align two sequences with BLOSUM62, compare them, and detail differences.

    Args:
        seq1 (str): First sequence (e.g., amino acid string).
        seq2 (str): Second sequence.
        name1 (str): Name of first sequence (default: "Seq1").
        name2 (str): Name of second sequence (default: "Seq2").
    """
    # Initialize aligner with BLOSUM62 matrix
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -1

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
    print(f"\nOriginal lengths: {name1} = {len_seq1}, {name2} = {len_seq2}")
    if len_seq2 > len_seq1:
        print(f"{name2} is longer than {name1} by {length_diff} residues.")
    elif len_seq1 > len_seq2:
        print(f"{name1} is longer than {name2} by {length_diff} residues.")
    else:
        print("Sequences are the same length.")

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
            print(f"Residues {extra_start + 1}-{len_seq2}: {extra_seq}")

    # Summary statistics
    length = len(aligned_seq1)
    percent_identity = (length - (mismatches + gaps)) / length * 100
    print(f"\nSummary:")
    print(f"Total differences:")
    print(f"  Mismatches: {mismatches}")
    print(f"  Gaps: {gaps}")
    print(f"  Length difference: {length_diff}")
    print(f"Alignment length: {length}")
    print(f"Percent identity: {percent_identity:.2f}%")


# Example usage
if __name__ == "__main__":
    # Sample sequences (seq2 is longer)
    seq1 = "QVQLVESGGGLVKPGGSLRLSCAASGFTFSDYYMSWIRQAPGKGLEWVSYITYSGSTIYYADSVKGRFTISRDNAKSSLYLQMNSLRAEDTAVYYCARDRGTTMVPFDYWGQGTLVTVSSDIQMTQSPSSLSASVGDRVTITCQASQDITNYLNWYQQKPGKAPKLLIYAASNLETGVPSRFSGSGSGTDFTFTISGLQPEDIATYYCQQYDNLPLTFGGGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC"
    seq2 = "QVQLVESGGGLVKPGGSLRLSCAASGFTFSDYYMSWIRQAPGKGLEWVSYITYSGSTIYYADSVKGRFTISRDNAKSSLYLQMNSLRAEDTAVYYCARDRGTTMVPFDYWGQGTLVTVSSDIQMTQSPSSLSASVGDRVTITCQASQDITNYLNWYQQKPGKAPKLLIYAASNLETGVPSRFSGSGSGTDFTFTISGLQPEDIATYYCQQYDNLPLTFGGGTKVEIK"

    # Run the alignment and comparison
    align_and_compare(seq1, seq2, "Seq1", "Seq2")