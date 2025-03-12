#This is to convert 3 letter AA to 1 letter AA
#!/usr/bin/env python3

# Sequence Cleaner - Removes non-standard characters from biological sequences

import re


def clean_dna(sequence):
    """Remove all non-DNA characters and return a clean DNA sequence."""
    # Keep only A, T, G, C (case insensitive)
    clean = re.sub(r'[^ATGCatgc]', '', sequence)
    return clean.upper()


def clean_rna(sequence):
    """Remove all non-RNA characters and return a clean RNA sequence."""
    # Keep only A, U, G, C (case insensitive)
    clean = re.sub(r'[^AUGCaugc]', '', sequence)
    return clean.upper()


def clean_protein(sequence):
    """Remove all non-amino acid characters and return a clean protein sequence."""
    # All standard amino acids in one-letter code
    valid_aa = 'ACDEFGHIKLMNPQRSTVWY'
    # Keep only valid amino acid letters (case insensitive)
    clean = re.sub(r'[^' + valid_aa + valid_aa.lower() + ']', '', sequence)
    return clean.upper()


def format_sequence(sequence, width=60):
    """Format a sequence with line breaks at specified width."""
    return '\n'.join(sequence[i:i + width] for i in range(0, len(sequence), width))


def main():
    print("Sequence Cleaner and Formatter")
    print("-" * 50)
    print("This script removes non-standard characters from biological sequences")
    print("and formats them properly.")
    print("\nInstructions:")
    print("1. Select sequence type (DNA, RNA, or Protein)")
    print("2. Enter your sequence (can be multiple lines)")
    print("3. Type a single period '.' on a new line to process the sequence")
    print("4. Type 'quit' to exit the program")
    print("-" * 50)

    while True:
        # Get sequence type
        print("\nSelect sequence type:")
        print("1. DNA (A, T, G, C)")
        print("2. RNA (A, U, G, C)")
        print("3. Protein (amino acids)")
        print("Type 'quit' to exit")

        choice = input("> ").lower().strip()

        if choice == 'quit':
            print("Exiting program.")
            break

        # Validate choice
        if choice not in ['1', '2', '3', 'dna', 'rna', 'protein']:
            print("Invalid choice. Please select 1, 2, 3 or type the sequence type.")
            continue

        # Map choice to sequence type
        seq_type = {'1': 'dna', '2': 'rna', '3': 'protein',
                    'dna': 'dna', 'rna': 'rna', 'protein': 'protein'}[choice]

        # Get the sequence
        print(f"\nEnter your {seq_type.upper()} sequence (end with a single '.' on a new line):")

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
        full_sequence = "".join(lines)

        # Clean the sequence based on type
        if seq_type == 'dna':
            cleaned = clean_dna(full_sequence)
        elif seq_type == 'rna':
            cleaned = clean_rna(full_sequence)
        else:  # protein
            cleaned = clean_protein(full_sequence)

        # Format and display the result
        formatted = format_sequence(cleaned)

        print("\nCleaned and formatted sequence:")
        print(formatted)
        print(f"Length: {len(cleaned)} characters")

        # Ask if user wants to save to file
        save = input("\nDo you want to save this sequence to a file? (y/n): ").lower()
        if save.startswith('y'):
            filename = input("Enter filename to save as: ")
            if not filename.endswith('.txt'):
                filename += '.txt'
            try:
                with open(filename, 'w') as f:
                    f.write(formatted)
                print(f"Sequence saved to {filename}")
            except Exception as e:
                print(f"Error saving file: {e}")


if __name__ == "__main__":
    main()