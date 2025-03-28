# GC and AU Rate Calculator for DNA/RNA Sequences
# Based on SeqCompare.py template

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


def dna_to_rna(dna_sequence):
    """
    Convert DNA sequence to RNA by replacing T with U
    """
    return dna_sequence.replace("T", "U")


def calculate_gc_content(sequence):
    """
    Calculate GC content percentage in a given sequence
    """
    gc_count = sequence.count("G") + sequence.count("C")
    total_count = len(sequence)

    # Avoid division by zero
    if total_count == 0:
        return 0.0

    gc_percentage = (gc_count / total_count) * 100
    return gc_percentage


def calculate_au_content(rna_sequence):
    """
    Calculate AU content percentage in an RNA sequence
    """
    au_count = rna_sequence.count("A") + rna_sequence.count("U")
    total_count = len(rna_sequence)

    # Avoid division by zero
    if total_count == 0:
        return 0.0

    au_percentage = (au_count / total_count) * 100
    return au_percentage


def get_sequence_input(seq_type):
    """Get multi-line sequence input from user."""
    print(f"\nEnter {seq_type} Sequence (end with a single '.' on a new line):")

    if seq_type == 'DNA':
        print("(DNA sequence using A, T, G, C nucleotides)")
    elif seq_type == 'RNA':
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


def analyze_sequence(sequence, seq_type):
    """
    Analyze the sequence to calculate GC and AU content
    """
    if seq_type == 'DNA':
        # Convert DNA to RNA for AU calculation
        rna_sequence = dna_to_rna(sequence)
        gc_content = calculate_gc_content(sequence)
        au_content = calculate_au_content(rna_sequence)
        return gc_content, au_content, rna_sequence
    else:  # RNA
        gc_content = calculate_gc_content(sequence)
        au_content = calculate_au_content(sequence)
        return gc_content, au_content, sequence


def display_results(sequence, seq_type, gc_content, au_content, rna_sequence=None):
    """
    Display the analysis results
    """
    print("\n" + "=" * 50)
    print(f"{seq_type} Sequence Analysis Results")
    print("=" * 50)

    # Display the cleaned sequence
    print(f"\nCleaned {seq_type} Sequence:")
    for i in range(0, len(sequence), 60):
        print(sequence[i:i + 60])

    # If DNA was converted to RNA, show the RNA sequence
    if seq_type == 'DNA':
        print(f"\nConverted RNA Sequence:")
        for i in range(0, len(rna_sequence), 60):
            print(rna_sequence[i:i + 60])

    # Display the nucleotide composition
    print("\nNucleotide Composition:")

    if seq_type == 'DNA':
        a_count = sequence.count("A")
        t_count = sequence.count("T")
        g_count = sequence.count("G")
        c_count = sequence.count("C")
        other_count = len(sequence) - (a_count + t_count + g_count + c_count)

        print(f"A: {a_count} ({a_count / len(sequence) * 100:.2f}%)")
        print(f"T: {t_count} ({t_count / len(sequence) * 100:.2f}%)")
        print(f"G: {g_count} ({g_count / len(sequence) * 100:.2f}%)")
        print(f"C: {c_count} ({c_count / len(sequence) * 100:.2f}%)")
        if other_count > 0:
            print(f"Other: {other_count} ({other_count / len(sequence) * 100:.2f}%)")
    else:  # RNA
        a_count = sequence.count("A")
        u_count = sequence.count("U")
        g_count = sequence.count("G")
        c_count = sequence.count("C")
        other_count = len(sequence) - (a_count + u_count + g_count + c_count)

        print(f"A: {a_count} ({a_count / len(sequence) * 100:.2f}%)")
        print(f"U: {u_count} ({u_count / len(sequence) * 100:.2f}%)")
        print(f"G: {g_count} ({g_count / len(sequence) * 100:.2f}%)")
        print(f"C: {c_count} ({c_count / len(sequence) * 100:.2f}%)")
        if other_count > 0:
            print(f"Other: {other_count} ({other_count / len(sequence) * 100:.2f}%)")

    # Display GC and AU content
    print("\nContent Analysis:")
    print(f"GC Content: {gc_content:.2f}%")
    print(f"AU Content: {au_content:.2f}%")

    print("\nSequence Length: {0} nucleotides".format(len(sequence)))
    print("=" * 50)


def main():
    print("GC and AU Rate Calculator for DNA/RNA Sequences")
    print("=" * 70)
    print("This tool analyzes nucleotide sequences to calculate:")
    print("• GC content (percentage of G and C nucleotides)")
    print("• AU content (percentage of A and U nucleotides)")
    print("\nFeatures:")
    print("• Handles both DNA and RNA sequences")
    print("• Automatically converts DNA to RNA for AU calculation")
    print("• Cleans and standardizes input sequences")
    print("• Provides detailed nucleotide composition analysis")
    print("\nInstructions:")
    print("1. Select the sequence type (DNA or RNA)")
    print("2. Enter your sequence (can be multiple lines)")
    print("3. Type a single period '.' on a new line to finish input")
    print("4. Type 'quit' at any prompt to exit the program")
    print("=" * 70)

    while True:
        # Select sequence type
        print("\nSelect sequence type to analyze:")
        print("1. DNA (nucleotide sequences with A, T, G, C)")
        print("2. RNA (nucleotide sequences with A, U, G, C)")
        print("Type 'quit' to exit")

        seq_choice = input("> ").lower().strip()
        if seq_choice == 'quit':
            print("Exiting program.")
            break

        # Map choice to sequence type
        if seq_choice in ['1', 'dna', 'd']:
            seq_type = 'DNA'
        elif seq_choice in ['2', 'rna', 'r']:
            seq_type = 'RNA'
        else:
            print("Invalid choice. Please select 1, 2 or type the sequence type.")
            continue

        print(f"\nYou selected: {seq_type} sequence analysis")

        # Get sequence input
        raw_seq = get_sequence_input(seq_type)
        if raw_seq is None:
            print("Exiting program.")
            break

        if not raw_seq:
            print("No sequence provided. Please try again.")
            continue

        # Clean sequence based on type
        if seq_type == 'DNA':
            clean_seq, changes = clean_dna_sequence(raw_seq)
        else:  # RNA
            clean_seq, changes = clean_rna_sequence(raw_seq)

        if changes:
            print(f"\nSequence was cleaned and processed.")
            print(f"Final sequence length: {len(clean_seq)} nucleotides")

        if not clean_seq:
            print(f"No valid {seq_type} characters found in sequence. Please try again.")
            continue

        # Analyze the sequence
        try:
            gc_content, au_content, rna_seq = analyze_sequence(clean_seq, seq_type)
            display_results(clean_seq, seq_type, gc_content, au_content, rna_seq)
        except Exception as e:
            print(f"Error during analysis: {e}")

        # Offer to save the results
        print("\nDo you want to save the analysis results? (y/n):")
        save = input("> ").lower()
        if save.startswith('y'):
            try:
                filename = input("Enter a filename (will save as .txt): ")
                if not filename.endswith('.txt'):
                    filename += '.txt'
                with open(filename, 'w') as f:
                    f.write(f"{seq_type} Sequence Analysis Results\n")
                    f.write("=" * 50 + "\n")
                    f.write(f"\nSequence:\n")
                    for i in range(0, len(clean_seq), 60):
                        f.write(clean_seq[i:i + 60] + "\n")

                    if seq_type == 'DNA':
                        f.write(f"\nConverted RNA Sequence:\n")
                        for i in range(0, len(rna_seq), 60):
                            f.write(rna_seq[i:i + 60] + "\n")

                    f.write("\nContent Analysis:\n")
                    f.write(f"GC Content: {gc_content:.2f}%\n")
                    f.write(f"AU Content: {au_content:.2f}%\n")
                    f.write(f"\nSequence Length: {len(clean_seq)} nucleotides\n")
                print(f"Results saved to {filename}")
            except Exception as e:
                print(f"Error saving file: {e}")

        # Ask if user wants to do another analysis
        print("\nDo you want to analyze another sequence? (y/n):")
        again = input("> ").lower()
        if not again.startswith('y'):
            print("Exiting program.")
            break


if __name__ == "__main__":
    main()