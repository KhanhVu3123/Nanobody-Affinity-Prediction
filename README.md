# Nanobody–Antigen Affinity Prediction

This project predicts nanobody–antigen binding using protein sequence embeddings and machine learning. It combines structural bioinformatics, protein language models, and classical + deep learning methods to classify binding interactions.

---

## Overview

The pipeline consists of:

1. Parsing PDB structures (nanobody + antigen chains)
2. Extracting amino acid sequences and CDR regions
3. Generating protein embeddings using **ESM-2**
4. Building datasets of true vs false binding pairs
5. Training machine learning models (ANN + classical ML)
6. Evaluating performance using classification metrics
7. Visualising embedding space using t-SNE

---

## Key Components

### 1. Sequence Extraction
Extracts nanobody and antigen sequences from PDB files and gzipped structures.

- Cα-based sequence reconstruction
- Chain separation (nanobody vs antigen)
- Support for multi-chain PDBs

---

### 2. CDR Detection
CDR regions are identified using framework alignment-based matching:

- Aligns nanobody sequence to canonical framework regions
- Locates CDR1, CDR2, and CDR3 boundaries
- Extracts variable regions for downstream modeling

---

### 3. Embedding Generation (ESM-2)

Protein sequences are encoded using Meta AI’s ESM-2 transformer:

- Model: `esm2_t36_3B_UR50D`
- Outputs residue-level embeddings
- Used as input features for ML models

---

### 4. Machine Learning Models

#### Deep Learning
- Artificial Neural Network (PyTorch)
- Hyperparameter tuning using Optuna
- Dropout regularisation

#### Classical Models
- Logistic Regression
- Support Vector Machine
- Random Forest
- K-Nearest Neighbours
- Naive Bayes
- Decision Tree

---

### 5. Evaluation Metrics

Models are evaluated using:

- Accuracy
- Precision
- Recall
- F1-score
- ROC-AUC

Results are averaged over multiple train-test splits.

---

### 6. Visualisation

- t-SNE projection of embedding space
- Separation of binding vs non-binding pairs

---

## Requirements

Install dependencies:

```bash
pip install torch numpy pandas scikit-learn matplotlib biopython optuna
