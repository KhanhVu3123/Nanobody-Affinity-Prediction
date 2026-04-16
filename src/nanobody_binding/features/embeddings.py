#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Generate ESM-2 embeddings from sequence files.
"""

import torch
import numpy as np
import gc


# -----------------------
# Config
# -----------------------
MODEL_NAME = "esm2_t36_3B_UR50D"
MAX_LEN = 1024

INPUT_FILES = [
    "data/sequence/New_Wrong_concat_seq.txt",
    "data/sequence/New_concat_seq.txt",
]

OUTPUT_FILES = [
    "data/embeddings/New_Wrong_concat_embeddings.npy",
    "data/embeddings/New_concat_embeddings.npy",
]


# -----------------------
# Load model
# -----------------------
print("Loading ESM model...")
model, alphabet = torch.hub.load("facebookresearch/esm", MODEL_NAME)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = model.to(device)
model.eval()

print("Using device:", device)


# -----------------------
# FASTA parser
# -----------------------
def load_sequences(file_path, max_len=1024):
    sequences = []
    seq = ""

    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()

            if line.startswith(">"):
                if seq and len(seq) <= max_len:
                    sequences.append(seq)
                seq = ""
            else:
                seq += line

        if seq and len(seq) <= max_len:
            sequences.append(seq)

    return sequences


# -----------------------
# Embedding function
# -----------------------
def compute_embeddings(sequences):
    embeddings = []

    for i, seq in enumerate(sequences):
        tokens = alphabet.encode(seq)
        tokens = torch.tensor([tokens], dtype=torch.long).to(device)

        with torch.no_grad():
            results = model(tokens, repr_layers=[36])
            reps = results["representations"][36].squeeze(0)

        embeddings.append(reps.cpu().numpy())

        if i % 10 == 0 or i == len(sequences) - 1:
            print(f"{i+1}/{len(sequences)} sequences processed")

        # cleanup
        del tokens, results, reps
        torch.cuda.empty_cache()
        gc.collect()

    return embeddings


# -----------------------
# Main loop
# -----------------------
for input_path, output_path in zip(INPUT_FILES, OUTPUT_FILES):

    print(f"\nProcessing: {input_path}")

    sequences = load_sequences(input_path, MAX_LEN)
    print(f"Loaded {len(sequences)} sequences")

    embeddings = compute_embeddings(sequences)

    np.save(output_path, embeddings)
    print(f"Saved embeddings to: {output_path}")


print("\nDone.")