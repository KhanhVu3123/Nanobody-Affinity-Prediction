#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Nanobody–Antigen Binding Classifier using Embeddings + PyTorch + Optuna
"""

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim

from torch.utils.data import DataLoader, TensorDataset
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc

import matplotlib.pyplot as plt
import optuna


# -----------------------
# Reproducibility
# -----------------------
def set_seed(seed=42):
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


set_seed()


# -----------------------
# Data Loading
# -----------------------
def load_embeddings(path):
    df = pd.read_csv(path)
    embeddings = []

    for line in df["Embeddings"]:
        values = line.strip()[1:-1].split(", ")
        embeddings.append([float(v) for v in values])

    return np.array(embeddings)


true_path = "data/embeddings/True_Embed.csv"
false_path = "data/embeddings/False_Embed.csv"

X_true = load_embeddings(true_path)
X_false = load_embeddings(false_path)

y_true = np.ones(len(X_true))
y_false = np.zeros(len(X_false))

X = np.vstack((X_true, X_false))
y = np.concatenate((y_true, y_false))


# -----------------------
# Device
# -----------------------
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print("Using device:", device)


# -----------------------
# Train-test split
# -----------------------
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42, stratify=y
)

X_train = torch.tensor(X_train, dtype=torch.float32)
y_train = torch.tensor(y_train, dtype=torch.float32)

X_test = torch.tensor(X_test, dtype=torch.float32)
y_test = torch.tensor(y_test, dtype=torch.float32)

train_loader = DataLoader(
    TensorDataset(X_train, y_train), batch_size=64, shuffle=True
)

val_loader = DataLoader(
    TensorDataset(X_test, y_test), batch_size=64, shuffle=False
)


# -----------------------
# Model
# -----------------------
class Net(nn.Module):
    def __init__(self, dropout_rate=0.5):
        super().__init__()

        self.net = nn.Sequential(
            nn.Linear(1280, 256),
            nn.ReLU(),
            nn.Dropout(dropout_rate),

            nn.Linear(256, 128),
            nn.ReLU(),
            nn.Dropout(dropout_rate),

            nn.Linear(128, 1),
            nn.Sigmoid()
        )

    def forward(self, x):
        return self.net(x).squeeze()


criterion = nn.BCELoss()


# -----------------------
# Train function
# -----------------------
def train_one_epoch(model, loader, optimizer):
    model.train()
    total_loss = 0

    for x, y in loader:
        x, y = x.to(device), y.to(device)

        optimizer.zero_grad()
        preds = model(x)
        loss = criterion(preds, y)

        loss.backward()
        optimizer.step()

        total_loss += loss.item()

    return total_loss / len(loader)


# -----------------------
# Evaluation
# -----------------------
def evaluate(model, loader):
    model.eval()

    y_true, y_scores = [], []

    with torch.no_grad():
        for x, y in loader:
            x = x.to(device)

            preds = model(x).cpu().numpy()

            y_scores.extend(preds)
            y_true.extend(y.numpy())

    y_true = np.array(y_true)
    y_scores = np.array(y_scores)

    preds_binary = (y_scores > 0.5).astype(int)
    accuracy = (preds_binary == y_true).mean()

    fpr, tpr, _ = roc_curve(y_true, y_scores)
    roc_auc = auc(fpr, tpr)

    return accuracy, fpr, tpr, roc_auc


# -----------------------
# Optuna objective
# -----------------------
def objective(trial):
    lr = trial.suggest_float("lr", 1e-5, 1e-2, log=True)
    dropout = trial.suggest_float("dropout_rate", 0.1, 0.5)

    model = Net(dropout).to(device)
    optimizer = optim.Adam(model.parameters(), lr=lr)

    for _ in range(10):
        train_one_epoch(model, train_loader, optimizer)

    acc, _, _, _ = evaluate(model, val_loader)
    return acc


study = optuna.create_study(direction="maximize")
study.optimize(objective, n_trials=50)

best = study.best_params
print("Best params:", best)


# -----------------------
# Final training
# -----------------------
model = Net(best["dropout_rate"]).to(device)
optimizer = optim.Adam(model.parameters(), lr=best["lr"])

for epoch in range(10):
    loss = train_one_epoch(model, train_loader, optimizer)
    print(f"Epoch {epoch+1}: loss={loss:.4f}")


# -----------------------
# Final evaluation
# -----------------------
accuracy, fpr, tpr, roc_auc = evaluate(model, val_loader)

print(f"Accuracy: {accuracy * 100:.2f}%")
print(f"AUC: {roc_auc:.4f}")


# -----------------------
# ROC curve
# -----------------------
plt.figure()
plt.plot(fpr, tpr, lw=2, label=f"ROC (AUC = {roc_auc:.3f})")
plt.plot([0, 1], [0, 1], linestyle="--")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Nanobody–Antigen Binding ROC Curve")
plt.legend()
plt.show()


# -----------------------
# Save model
# -----------------------
torch.save(model.state_dict(), "NbAgBindingModel.pth")