"""
ProtTrans for protein analysis on M2 MacBook Air
------------------------------------------------
This script demonstrates how to use ProtTrans models for protein analysis--get embeddings.
"""

# Required packages - Install in PyCharm terminal:
# pip install torch transformers sentencepiece
# pip install matplotlib seaborn scikit-learn  # For visualization
# pip install biopython  # For FASTA handling

import torch
from transformers import T5EncoderModel, T5Tokenizer
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA


def setup_device():
    """Set up computation device (MPS for Apple Silicon, CPU otherwise)"""
    try:
        if torch.backends.mps.is_available():
            device = torch.device("mps")
            print("Using MPS (Metal) device for Apple Silicon")
            return device
        else:
            print("MPS not available, using CPU")
            return torch.device("cpu")
    except Exception as e:
        print(f"Error setting up MPS: {e}")
        print("Falling back to CPU")
        return torch.device("cpu")


def load_prot_t5_model(model_name="Rostlab/prot_t5_xl_half_uniref50-enc"):
    """Load a ProtT5 model from Hugging Face"""
    print(f"Loading {model_name}...")
    tokenizer = T5Tokenizer.from_pretrained(model_name, do_lower_case=False)

    # Use half precision for efficiency on M2 chip
    model = T5EncoderModel.from_pretrained(model_name, torch_dtype=torch.float16)
    model = model.to(device)
    model.eval()
    print(f"{model_name} loaded successfully")
    return model, tokenizer


def get_protein_embeddings(model, tokenizer, sequences, batch_size=2):
    """Get embeddings from a ProtT5 model"""
    # Format sequences with spaces between amino acids
    sequences_processed = [" ".join(list(seq)) for seq in sequences]

    all_embeddings = []

    # Process in batches
    for i in range(0, len(sequences_processed), batch_size):
        batch_seqs = sequences_processed[i:i + batch_size]
        print(f"Processing batch {i // batch_size + 1}/{(len(sequences_processed) + batch_size - 1) // batch_size}")

        # Tokenize sequences and get model output
        ids = tokenizer.batch_encode_plus(batch_seqs, add_special_tokens=True, padding="longest")
        input_ids = torch.tensor(ids['input_ids']).to(device)
        attention_mask = torch.tensor(ids['attention_mask']).to(device)

        with torch.no_grad():
            embedding_repr = model(input_ids=input_ids, attention_mask=attention_mask)

        # Get sequence embeddings (mean pooling)
        for b, seq_len in enumerate(attention_mask.sum(dim=1)):
            seq_emd = embedding_repr.last_hidden_state[b, :seq_len]
            all_embeddings.append(seq_emd.mean(dim=0).cpu().numpy())

    return np.array(all_embeddings)


def load_sequences_from_file(file_path):
    """Load protein sequences from a file"""
    sequences = []
    sequence_ids = []

    # Try to detect file format
    if file_path.endswith(('.fasta', '.fa', '.faa')):
        # Load sequences from FASTA file
        from Bio import SeqIO
        for record in SeqIO.parse(file_path, "fasta"):
            sequences.append(str(record.seq))
            sequence_ids.append(record.id)
    else:
        # Assume it's a text file with one sequence per line
        with open(file_path, 'r') as f:
            for i, line in enumerate(f):
                line = line.strip()
                if line:  # Skip empty lines
                    sequences.append(line)
                    sequence_ids.append(f"Sequence_{i + 1}")

    print(f"Loaded {len(sequences)} sequences from {file_path}")
    return sequences, sequence_ids


def visualize_embeddings(embeddings, labels=None):
    """Visualize protein embeddings with PCA"""
    # Apply PCA for dimensionality reduction
    pca = PCA(n_components=2)
    embedding_2d = pca.fit_transform(embeddings)

    # Create plot
    plt.figure(figsize=(10, 8))
    sns.set_style("whitegrid")

    if labels is None or len(labels) != len(embeddings):
        plt.scatter(embedding_2d[:, 0], embedding_2d[:, 1], s=100, alpha=0.7)
    else:
        # If we have fewer than 10 unique labels, color by label
        unique_labels = np.unique(labels)
        if len(unique_labels) < 10:
            for label in unique_labels:
                mask = np.array(labels) == label
                plt.scatter(embedding_2d[mask, 0], embedding_2d[mask, 1], s=100, alpha=0.7, label=label)
            plt.legend()
        else:
            plt.scatter(embedding_2d[:, 0], embedding_2d[:, 1], s=100, alpha=0.7)
            # Add annotations for points if not too many
            if len(embeddings) <= 30:
                for i, label in enumerate(labels):
                    plt.annotate(label, (embedding_2d[i, 0], embedding_2d[i, 1]), fontsize=8)

    plt.title("Protein Embeddings PCA", fontsize=16)
    plt.xlabel("PC1", fontsize=14)
    plt.ylabel("PC2", fontsize=14)
    plt.tight_layout()
    plt.savefig("protein_embeddings.png")
    print("Visualization saved to protein_embeddings.png")
    plt.show()


def calculate_sequence_similarity(embeddings, sequence_ids=None):
    """Calculate similarity between protein sequences"""
    from scipy.spatial.distance import pdist, squareform

    # Calculate cosine distances between sequences
    distances = pdist(embeddings, metric="cosine")
    # Convert to similarities (1 - distance)
    similarities = 1 - squareform(distances)

    # Visualize similarity matrix
    plt.figure(figsize=(10, 8))
    labels = sequence_ids if sequence_ids and len(sequence_ids) <= 20 else None
    sns.heatmap(similarities, cmap="viridis", vmin=0, vmax=1,
                xticklabels=labels, yticklabels=labels, annot=True if len(embeddings) <= 10 else False)
    plt.title("Protein Sequence Similarity", fontsize=16)
    plt.tight_layout()
    plt.savefig("sequence_similarity.png")
    print("Similarity matrix saved to sequence_similarity.png")
    plt.show()

    return similarities


def example_sequence_analysis():
    """Example: Calculate similarity between protein sequences"""
    # Example protein sequences (Hemoglobin variants)
    sequences = [
        "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
        # Hemoglobin alpha
        "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH",
        # Hemoglobin beta
        "MGHFTEEDKATITSLWGKVNVEDAGGETLGRLLVVYPWTQRFFDSFGNLSSASAIMGNPKVKAHGKKVLTSLGDAIKHLDDLKGTFAQLSELHCDKLHVDPENFKLLGNVLVTVLAIHFGKEFTPEVQASWQKMVTGVASALSSRYH"
        # Hemoglobin delta
    ]
    sequence_ids = ["HBA", "HBB", "HBD"]

    print("Computing embeddings for example sequences...")
    embeddings = get_protein_embeddings(model, tokenizer, sequences)

    # Visualize embeddings
    visualize_embeddings(embeddings, labels=sequence_ids)

    # Calculate similarity
    similarities = calculate_sequence_similarity(embeddings, sequence_ids)

    return embeddings, similarities


def analyze_fasta_file(file_path):
    """Analyze sequences from a FASTA file"""
    # Load sequences from file
    sequences, sequence_ids = load_sequences_from_file(file_path)

    if not sequences:
        print("No sequences found in the file.")
        return

    # Get embeddings
    embeddings = get_protein_embeddings(model, tokenizer, sequences)

    # Save embeddings
    np.save(f"{os.path.splitext(file_path)[0]}_embeddings.npy", embeddings)
    print(f"Embeddings saved to {os.path.splitext(file_path)[0]}_embeddings.npy")

    # Visualize embeddings
    visualize_embeddings(embeddings, labels=sequence_ids)

    # If we have multiple sequences, calculate similarity
    if len(sequences) > 1:
        similarities = calculate_sequence_similarity(embeddings, sequence_ids)

    return embeddings


def main():
    """Main function"""
    global device, model, tokenizer

    # Set up device
    device = setup_device()

    # Load ProtT5 model
    model, tokenizer = load_prot_t5_model()

    # Provide examples or analyze a file
    import sys
    if len(sys.argv) > 1:
        file_path = sys.argv[1]
        print(f"\n--- Analyzing File: {file_path} ---")
        analyze_fasta_file(file_path)
    else:
        print("\n--- Running Example Analysis ---")
        example_sequence_analysis()

    print("\nAnalysis completed!")


if __name__ == "__main__":
    main()