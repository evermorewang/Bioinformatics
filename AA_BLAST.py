from Bio import pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.Blast import NCBIWWW, NCBIXML
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


def align_sequences(seq1, seq2, seq1_name="Sequence 1", seq2_name="Sequence 2"):
    """Align two protein sequences and display their differences"""

    print(f"Aligning sequences...")
    print(f"{seq1_name}: {len(seq1)} amino acids")
    print(f"{seq2_name}: {len(seq2)} amino acids")

    # Use BLOSUM62 matrix for protein alignment
    matrix = substitution_matrices.load("BLOSUM62")

    # Perform global alignment with pairwise2
    alignments = pairwise2.align.globalds(seq1, seq2, matrix, -10, -0.5)

    # Get the best alignment
    best_alignment = alignments[0]
    aligned_seq1, aligned_seq2, score, start, end = best_alignment

    print(f"\nAlignment score: {score:.1f}")

    # Calculate sequence identity and similarity
    identical = 0
    similar = 0
    aligned = 0

    # Identify differences and similarities
    differences = []
    for i in range(len(aligned_seq1)):
        if aligned_seq1[i] != '-' and aligned_seq2[i] != '-':
            aligned += 1
            if aligned_seq1[i] == aligned_seq2[i]:
                identical += 1
                differences.append((i, aligned_seq1[i], aligned_seq2[i], "identical"))
            elif (aligned_seq1[i], aligned_seq2[i]) in matrix and matrix[(aligned_seq1[i], aligned_seq2[i])] > 0:
                similar += 1
                differences.append((i, aligned_seq1[i], aligned_seq2[i], "similar"))
            else:
                differences.append((i, aligned_seq1[i], aligned_seq2[i], "different"))
        else:
            differences.append((i, aligned_seq1[i], aligned_seq2[i], "gap"))

    identity = identical / aligned * 100 if aligned > 0 else 0
    similarity = (identical + similar) / aligned * 100 if aligned > 0 else 0

    print(f"Sequence identity: {identity:.1f}%")
    print(f"Sequence similarity: {similarity:.1f}%")

    # Print alignment in blocks
    block_size = 60
    for i in range(0, len(aligned_seq1), block_size):
        block_end = min(i + block_size, len(aligned_seq1))

        # Create match line
        match_line = ""
        for j in range(i, block_end):
            if aligned_seq1[j] == aligned_seq2[j]:
                match_line += "|"
            elif (aligned_seq1[j], aligned_seq2[j]) in matrix and matrix[(aligned_seq1[j], aligned_seq2[j])] > 0:
                match_line += "+"
            else:
                match_line += " "

        print(f"\n{seq1_name} {i + 1:4d} {aligned_seq1[i:block_end]}")
        print(f"{'':10s} {match_line}")
        print(f"{seq2_name} {i + 1:4d} {aligned_seq2[i:block_end]}")

    # Visualize the alignment with colorful differences
    visualize_alignment(aligned_seq1, aligned_seq2, differences, seq1_name, seq2_name)

    return differences, identity, similarity


def visualize_alignment(seq1, seq2, differences, seq1_name, seq2_name):
    """Create a graphical representation of the sequence alignment"""

    # Define colors
    identical_color = 'lightgreen'
    similar_color = 'lightblue'
    different_color = 'salmon'
    gap_color = 'lightgray'

    fig, ax = plt.subplots(figsize=(12, 4))

    # Plot the color-coded alignment
    positions = np.arange(len(seq1))
    heights = np.ones(len(seq1))

    # Draw bars for sequence 1 (top)
    colors1 = []
    for i, (_, aa1, aa2, diff_type) in enumerate(differences):
        if diff_type == "identical":
            colors1.append(identical_color)
        elif diff_type == "similar":
            colors1.append(similar_color)
        elif diff_type == "different":
            colors1.append(different_color)
        else:  # gap
            colors1.append(gap_color)

    # Draw bars for sequence 2 (bottom)
    colors2 = []
    for i, (_, aa1, aa2, diff_type) in enumerate(differences):
        if diff_type == "identical":
            colors2.append(identical_color)
        elif diff_type == "similar":
            colors2.append(similar_color)
        elif diff_type == "different":
            colors2.append(different_color)
        else:  # gap
            colors2.append(gap_color)

    # Draw the bars
    ax.bar(positions, heights, width=1.0, color=colors1, edgecolor='gray', linewidth=0.5, alpha=0.8)
    ax.bar(positions, -heights, width=1.0, color=colors2, edgecolor='gray', linewidth=0.5, alpha=0.8)

    # Add sequence labels
    ax.text(-5, 0.5, seq1_name, fontsize=12, verticalalignment='center')
    ax.text(-5, -0.5, seq2_name, fontsize=12, verticalalignment='center')

    # Add legend
    legend_elements = [
        mpatches.Patch(facecolor=identical_color, edgecolor='gray', label='Identical'),
        mpatches.Patch(facecolor=similar_color, edgecolor='gray', label='Similar'),
        mpatches.Patch(facecolor=different_color, edgecolor='gray', label='Different'),
        mpatches.Patch(facecolor=gap_color, edgecolor='gray', label='Gap')
    ]
    ax.legend(handles=legend_elements, loc='upper right')

    # Set the axis limits and labels
    ax.set_xlim(-10, len(seq1) + 10)
    ax.set_ylim(-1.5, 1.5)
    ax.set_xlabel('Alignment Position')
    ax.set_title('Sequence Alignment Comparison')

    # Remove y-axis ticks
    ax.set_yticks([])

    # Show the plot
    plt.tight_layout()
    plt.show()


def run_blast(sequence, database="nr", program="blastp", expect=10, hitlist_size=10):
    """Run a BLAST search against NCBI database"""
    print(f"Running BLAST search (this may take a few minutes)...")
    result_handle = NCBIWWW.qblast(program, database, sequence, expect=expect, hitlist_size=hitlist_size)

    # Parse the results
    blast_record = NCBIXML.read(result_handle)

    # Print BLAST results
    print(f"\nBLAST Results:")
    print(f"Query length: {blast_record.query_length} amino acids")

    for i, alignment in enumerate(blast_record.alignments):
        for hsp in alignment.hsps:
            print(f"\nHit #{i + 1}: {alignment.title[:70]}...")
            print(f"Alignment length: {hsp.align_length}")
            print(f"E-value: {hsp.expect:.2e}")
            print(f"Score: {hsp.score} bits ({hsp.score_bits:.1f} bits)")
            print(f"Identities: {hsp.identities}/{hsp.align_length} ({hsp.identities / hsp.align_length * 100:.1f}%)")
            print(f"Positives: {hsp.positives}/{hsp.align_length} ({hsp.positives / hsp.align_length * 100:.1f}%)")
            print(f"Gaps: {hsp.gaps}/{hsp.align_length} ({hsp.gaps / hsp.align_length * 100:.1f}%)")

            # Print a snippet of the alignment
            print(f"Query: {hsp.query[:60]}...")
            print(f"Match: {hsp.match[:60]}...")
            print(f"Sbjct: {hsp.sbjct[:60]}...")

    return blast_record


# Example usage
if __name__ == "__main__":
    # Example protein sequences (insulin from human and mouse)
    Seq1 = "QVQLVESGGGLVKPGGSLRLSCAASGFTFSDYYMSWIRQAPGKGLEWVSYITYSGSTIYYADSVKGRFTISRDNAKSSLYLQMNSLRAEDTAVYYCARDRGTTMVPFDYWGQGTLVTVSSDIQMTQSPSSLSASVGDRVTITCQASQDITNYLNWYQQKPGKAPKLLIYAASNLETGVPSRFSGSGSGTDFTFTISGLQPEDIATYYCQQYDNLPLTFGGGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC"
    Seq2 = "QVQLVESGGGLVKPGGSLRLSCAASGFTFSDYYMSWIRQAPGKGLEWVSYITYSGSTIYYADSVKGRFTISRDNAKSSLYLQMNSLRAEDTAVYYCARDRGTTMVPFDYWGQGTLVTVSSDIQMTQSPSSLSASVGDRVTITCQASQDITNYLNWYQQKPGKAPKLLIYAASNLETGVPSRFSGSGSGTDFTFTISGLQPEDIATYYCQQYDNLPLTFGGGTKVEIK"

    # You can replace these with your own sequences
    # Either directly:
    # seq1 = "YOURLONGPROTEINSSEQUENCEHERE..."
    # seq2 = "YOURLONGPROTEINSSEQUENCEHERE..."

    # Or load from FASTA files:
    # from Bio import SeqIO
    # seq1 = str(SeqIO.read("sequence1.fasta", "fasta").seq)
    # seq2 = str(SeqIO.read("sequence2.fasta", "fasta").seq)

    # Perform the alignment
    align_sequences(Seq1, Seq2, "Seq1", "Seq2")

    # Uncomment to run a BLAST search against NCBI database (requires internet connection)
    # blast_record = run_blast(human_insulin)