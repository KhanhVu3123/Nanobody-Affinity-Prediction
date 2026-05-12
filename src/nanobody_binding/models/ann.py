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
optuna.logging.set_verbosity(optuna.logging.WARNING)


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

X_train_t = torch.tensor(X_train, dtype=torch.float32)
y_train_t = torch.tensor(y_train, dtype=torch.float32)

X_test_t = torch.tensor(X_test, dtype=torch.float32)
y_test_t = torch.tensor(y_test, dtype=torch.float32)


# -----------------------
# Model
# -----------------------
def build_model(n_layers, hidden_dim, dropout_rate):
    """Build a variable-depth MLP. Architecture is now tunable via Optuna."""
    layers = []
    in_dim = 1280

    for _ in range(n_layers):
        layers += [nn.Linear(in_dim, hidden_dim), nn.ReLU(), nn.Dropout(dropout_rate)]
        in_dim = hidden_dim

    layers.append(nn.Linear(in_dim, 1))
    layers.append(nn.Sigmoid())

    model = nn.Sequential(*layers)
    _kaiming_init(model)
    return model


def _kaiming_init(model):
    """He initialisation for all Linear layers — better convergence with ReLU."""
    for m in model.modules():
        if isinstance(m, nn.Linear):
            nn.init.kaiming_uniform_(m.weight, nonlinearity="relu")
            nn.init.zeros_(m.bias)


criterion = nn.BCELoss()


# -----------------------
# Early stopping
# -----------------------
class EarlyStopping:
    def __init__(self, patience=10):
        self.patience = patience
        self.best_loss = float("inf")
        self.counter = 0
        self.best_state = None

    def step(self, val_loss, model):
        if val_loss < self.best_loss:
            self.best_loss = val_loss
            self.counter = 0
            self.best_state = {k: v.clone() for k, v in model.state_dict().items()}
        else:
            self.counter += 1

        return self.counter >= self.patience

    def restore(self, model):
        if self.best_state is not None:
            model.load_state_dict(self.best_state)


# -----------------------
# Train / eval functions
# -----------------------
def train_one_epoch(model, loader, optimizer):
    model.train()
    total_loss = 0

    for x, y in loader:
        x, y = x.to(device), y.to(device)

        optimizer.zero_grad()
        preds = model(x).squeeze()
        loss = criterion(preds, y)

        loss.backward()
        optimizer.step()

        total_loss += loss.item()

    return total_loss / len(loader)


def eval_loss(model, loader):
    model.eval()
    total_loss = 0

    with torch.no_grad():
        for x, y in loader:
            x, y = x.to(device), y.to(device)
            preds = model(x).squeeze()
            total_loss += criterion(preds, y).item()

    return total_loss / len(loader)


def evaluate(model, loader):
    model.eval()

    y_true_list, y_scores = [], []

    with torch.no_grad():
        for x, y in loader:
            x = x.to(device)
            preds = model(x).squeeze().cpu().numpy()
            y_scores.extend(np.atleast_1d(preds))
            y_true_list.extend(y.numpy())

    y_true_arr = np.array(y_true_list)
    y_scores_arr = np.array(y_scores)

    preds_binary = (y_scores_arr > 0.5).astype(int)
    accuracy = (preds_binary == y_true_arr).mean()

    fpr, tpr, _ = roc_curve(y_true_arr, y_scores_arr)
    roc_auc = auc(fpr, tpr)

    return accuracy, fpr, tpr, roc_auc


# -----------------------
# Optuna objective
# -----------------------
def objective(trial):
    lr = trial.suggest_float("lr", 1e-5, 1e-2, log=True)
    dropout = trial.suggest_float("dropout_rate", 0.1, 0.5)
    n_layers = trial.suggest_int("n_layers", 1, 4)
    hidden_dim = trial.suggest_categorical("hidden_dim", [64, 128, 256, 512])
    batch_size = trial.suggest_categorical("batch_size", [32, 64, 128])

    train_loader = DataLoader(
        TensorDataset(X_train_t, y_train_t), batch_size=batch_size, shuffle=True
    )
    val_loader = DataLoader(
        TensorDataset(X_test_t, y_test_t), batch_size=batch_size, shuffle=False
    )

    model = build_model(n_layers, hidden_dim, dropout).to(device)
    optimizer = optim.Adam(model.parameters(), lr=lr)
    # Halve LR when val loss plateaus for 5 consecutive epochs
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, patience=5, factor=0.5)
    stopper = EarlyStopping(patience=10)

    for _ in range(50):
        train_one_epoch(model, train_loader, optimizer)
        val_loss = eval_loss(model, val_loader)
        scheduler.step(val_loss)

        if stopper.step(val_loss, model):
            break

    stopper.restore(model)
    acc, _, _, _ = evaluate(model, val_loader)
    return acc


study = optuna.create_study(direction="maximize")
study.optimize(objective, n_trials=50)

best = study.best_params
print("Best params:", best)


# -----------------------
# Final training
# -----------------------
EPOCHS = 100

train_loader = DataLoader(
    TensorDataset(X_train_t, y_train_t), batch_size=best["batch_size"], shuffle=True
)
val_loader = DataLoader(
    TensorDataset(X_test_t, y_test_t), batch_size=best["batch_size"], shuffle=False
)

model = build_model(best["n_layers"], best["hidden_dim"], best["dropout_rate"]).to(device)
optimizer = optim.Adam(model.parameters(), lr=best["lr"])
scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, patience=5, factor=0.5)
stopper = EarlyStopping(patience=15)

train_losses, val_losses = [], []

for epoch in range(EPOCHS):
    t_loss = train_one_epoch(model, train_loader, optimizer)
    v_loss = eval_loss(model, val_loader)

    scheduler.step(v_loss)
    train_losses.append(t_loss)
    val_losses.append(v_loss)

    print(f"Epoch {epoch + 1}: train_loss={t_loss:.4f}  val_loss={v_loss:.4f}")

    if stopper.step(v_loss, model):
        print(f"Early stopping at epoch {epoch + 1}")
        break

stopper.restore(model)


# -----------------------
# Final evaluation
# -----------------------
accuracy, fpr, tpr, roc_auc = evaluate(model, val_loader)

print(f"\nAccuracy: {accuracy * 100:.2f}%")
print(f"AUC:      {roc_auc:.4f}")


# -----------------------
# Plots
# -----------------------
fig, axes = plt.subplots(1, 2, figsize=(12, 4))

axes[0].plot(train_losses, label="Train")
axes[0].plot(val_losses, label="Val")
axes[0].set_xlabel("Epoch")
axes[0].set_ylabel("Loss")
axes[0].set_title("Training Curve")
axes[0].legend()

axes[1].plot(fpr, tpr, lw=2, label=f"ROC (AUC = {roc_auc:.3f})")
axes[1].plot([0, 1], [0, 1], linestyle="--")
axes[1].set_xlabel("False Positive Rate")
axes[1].set_ylabel("True Positive Rate")
axes[1].set_title("Nanobody–Antigen Binding ROC Curve")
axes[1].legend()

plt.tight_layout()
plt.show()


# -----------------------
# Save model
# -----------------------
torch.save(model.state_dict(), "NbAgBindingModel.pth")
