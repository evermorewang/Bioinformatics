"""
Protein Sequence Analysis using ProtBERT
----------------------------------------
A reliable transformer-based protein embedding model
"""

import torch
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from transformers import BertTokenizer, BertModel
from Bio import SeqIO
import argparse
import time
import sys


# Set up device for computation
def setup_device():
    if torch.cuda.is_available():
        return torch.device("cuda")
    elif hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")  # Apple Silicon
    else:
        return torch.device("cpu")


device = setup_device()
print(f"Using device: {device}")


def load_model(model_name="Rostlab/prot_bert"):
    """Load a protein language model (smaller than ProtT5-XL)"""
    print(f"Loading {model_name}...")
    start_time = time.time()

    # Load tokenizer and model
    tokenizer = BertTokenizer.from_pretrained(model_name, do_lower_case=False)
    model = BertModel.from_pretrained(model_name)
    model = model.to(device)
    model.eval()

    elapsed = time.time() - start_time
    print(f"Model loaded in {elapsed:.1f} seconds")
    return model, tokenizer


def embed_proteins(model, tokenizer, sequences, max_length=1024):
    """Generate embeddings for protein sequences"""
    embeddings = []
    ids = []

    for i, seq in enumerate(sequences):
        # Truncate sequence if needed
        if len(seq) > max_length - 2:  # Account for [CLS] and [SEP]
            print(f"Sequence {i + 1} truncated from {len(seq)} to {max_length - 2} residues")
            seq = seq[:max_length - 2]

        # Add special tokens and convert to model inputs
        encoded_input = tokenizer.encode_plus(
            ' '.join(seq),
            add_special_tokens=True,
            max_length=max_length,
            padding='max_length',
            truncation=True,
            return_tensors='pt'
        )

        # Move to device
        encoded_input = {k: v.to(device) for k, v in encoded_input.items()}

        # Generate embeddings
        print(f"Processing sequence {i + 1}/{len(sequences)}")
        with torch.no_grad():
            output = model(**encoded_input)

        # Get the [CLS] token embedding as the sequence representation
        embedding = output.last_hidden_state[:, 0, :].cpu().numpy()
        embeddings.append(embedding[0])
        ids.append(f"Seq_{i + 1}")

    return np.array(embeddings), ids


def load_fasta(file_path):
    """Load sequences from a FASTA file"""
    sequences = []
    ids = []

    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append(str(record.seq))
        ids.append(record.id)

    print(f"Loaded {len(sequences)} sequences from {file_path}")
    return sequences, ids


def visualize_embeddings(embeddings, ids=None):
    """Create a PCA visualization of protein embeddings"""
    # Apply PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(embeddings)

    # Create plot
    plt.figure(figsize=(10, 8))
    plt.scatter(pca_result[:, 0], pca_result[:, 1], s=100, alpha=0.7)

    # Add labels if not too many points
    if ids and len(ids) <= 20:
        for i, txt in enumerate(ids):
            plt.annotate(txt, (pca_result[i, 0], pca_result[i, 1]), fontsize=9)

    plt.title('Protein Embeddings - PCA Visualization', fontsize=15)
    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%} variance)')
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%} variance)')
    plt.grid(alpha=0.3)

    # Save and show plot
    output_file = "protein_embeddings_pca.png"
    plt.savefig(output_file, dpi=300)
    print(f"PCA visualization saved to {output_file}")
    plt.show()

    return pca_result


def calculate_similarities(embeddings, ids=None):
    """Calculate and visualize pairwise similarities between proteins"""
    # Normalize embeddings for cosine similarity
    normalized = embeddings / np.linalg.norm(embeddings, axis=1, keepdims=True)

    # Calculate cosine similarity matrix
    similarity_matrix = np.dot(normalized, normalized.T)

    # Visualize as heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(
        similarity_matrix,
        annot=len(embeddings) <= 15,
        fmt=".2f",
        cmap="viridis",
        xticklabels=ids if ids and len(ids) <= 20 else False,
        yticklabels=ids if ids and len(ids) <= 20 else False
    )
    plt.title("Protein Sequence Similarity", fontsize=15)

    # Save and show plot
    output_file = "protein_similarity_matrix.png"
    plt.savefig(output_file, dpi=300)
    print(f"Similarity matrix saved to {output_file}")
    plt.show()

    return similarity_matrix


def analyze_file(file_path):
    """Analyze protein sequences from a file"""
    # Load sequences
    if file_path.lower().endswith(('.fasta', '.fa', '.faa')):
        sequences, seq_ids = load_fasta(file_path)
    else:
        with open(file_path, 'r') as f:
            sequences = [line.strip() for line in f if line.strip()]
        seq_ids = [f"Seq_{i + 1}" for i in range(len(sequences))]

    if not sequences:
        print("No sequences found in the file.")
        return

    # Generate embeddings
    model, tokenizer = load_model()
    embeddings, _ = embed_proteins(model, tokenizer, sequences)

    # Save embeddings
    output_name = os.path.splitext(os.path.basename(file_path))[0]
    np.save(f"{output_name}_embeddings.npy", embeddings)
    print(f"Embeddings saved to {output_name}_embeddings.npy")

    # Visualize results
    visualize_embeddings(embeddings, seq_ids)

    # Calculate similarities if multiple sequences
    if len(sequences) > 1:
        calculate_similarities(embeddings, seq_ids)


def main():
    parser = argparse.ArgumentParser(description="Protein sequence analysis using transformer models")
    parser.add_argument("file", nargs="?", help="Path to FASTA or text file with protein sequences")
    args = parser.parse_args()

    if args.file:
        analyze_file(args.file)
    else:
        # Example sequences if no file provided
        print("No file provided. Running example analysis...")
        example_sequences = [
            "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGK",  # HBA (partial)
            "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKV",  # HBB (partial)
            "MGHFTEEDKATITSLWGKVNVEDAGGETLGRLLVVYPWTQRFFDSFGNLSSASAIMGNPKV"  # HBD (partial)
        ]
        example_ids = ["HBA", "HBB", "HBD"]

        model, tokenizer = load_model()
        embeddings, _ = embed_proteins(model, tokenizer, example_sequences)
        visualize_embeddings(embeddings, example_ids)
        calculate_similarities(embeddings, example_ids)


if __name__ == "__main__":
    main()