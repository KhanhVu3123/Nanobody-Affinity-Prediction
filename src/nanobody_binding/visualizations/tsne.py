#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
t-SNE visualization of embedding vectors.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE


# -----------------------
# Config
# -----------------------
INPUT_FILE = "data/embeddings/embeddings.csv"
OUTPUT_FILE = "results/tSNEPlot_of_antigen.pdf"


# -----------------------
# Data loading
# -----------------------
def load_embeddings(csv_path):
    """Load embeddings from CSV into NumPy array."""
    df = pd.read_csv(csv_path)

    embeddings = df["Embeddings"].apply(
        lambda x: np.fromstring(x.strip("[]"), sep=",")
    )

    return np.stack(embeddings.values)


# -----------------------
# t-SNE
# -----------------------
def run_tsne(embeddings, perplexity=30, n_iter=300):
    """Run t-SNE dimensionality reduction."""
    tsne = TSNE(
        n_components=2,
        perplexity=perplexity,
        n_iter=n_iter,
        random_state=42
    )
    return tsne.fit_transform(embeddings)


# -----------------------
# Plotting
# -----------------------
def plot_tsne(results, output_path):
    """Create and save t-SNE scatter plot."""
    plt.figure(figsize=(10, 6))
    plt.scatter(results[:, 0], results[:, 1], s=10)

    plt.xlabel("t-SNE Component 1")
    plt.ylabel("t-SNE Component 2")
    plt.title("t-SNE of Embeddings")

    plt.savefig(output_path, format="pdf", bbox_inches="tight")
    plt.show()


# -----------------------
# Main
# -----------------------
def main():
    embeddings = load_embeddings(INPUT_FILE)
    tsne_results = run_tsne(embeddings)
    plot_tsne(tsne_results, OUTPUT_FILE)


if __name__ == "__main__":
    main()