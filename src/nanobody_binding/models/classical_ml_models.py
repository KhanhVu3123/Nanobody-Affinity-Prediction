#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compare multiple ML algorithms on nanobody–antigen embedding data.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score,
    f1_score, roc_auc_score
)

from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.tree import DecisionTreeClassifier


# -----------------------
# Config
# -----------------------
TRUE_PATH = "data/embeddings/True_Embed.csv"
FALSE_PATH = "data/embeddings/False_Embed.csv"
N_SPLITS = 10
RESULTS_PATH = "results/classical_ml_results.csv"


# -----------------------
# Data loading
# -----------------------
def load_embeddings(csv_path):
    """Convert string embeddings into NumPy array."""
    df = pd.read_csv(csv_path)

    embeddings = df["Embeddings"].apply(
        lambda x: np.fromstring(x.strip("[]"), sep=",")
    )

    return np.stack(embeddings.values)


def prepare_data(true_path, false_path):
    X_true = load_embeddings(true_path)
    X_false = load_embeddings(false_path)

    y_true = np.ones(len(X_true))
    y_false = np.zeros(len(X_false))

    X = np.vstack((X_true, X_false))
    y = np.concatenate((y_true, y_false))

    return X, y


# -----------------------
# Models
# -----------------------
def get_classifiers():
    # class_weight='balanced' compensates for unequal true/false pair counts
    return {
        "Logistic Regression": LogisticRegression(max_iter=1000, class_weight="balanced"),
        "SVM": SVC(probability=True, class_weight="balanced"),
        "Random Forest": RandomForestClassifier(n_estimators=200, class_weight="balanced"),
        "KNN": KNeighborsClassifier(n_neighbors=7),
        "Naive Bayes": GaussianNB(),
        "Decision Tree": DecisionTreeClassifier(class_weight="balanced"),
    }


# -----------------------
# Evaluation
# -----------------------
def evaluate_models(X, y, n_splits=10):
    classifiers = get_classifiers()

    metrics = {
        "accuracy": {name: [] for name in classifiers},
        "precision": {name: [] for name in classifiers},
        "recall": {name: [] for name in classifiers},
        "f1": {name: [] for name in classifiers},
        "auc": {name: [] for name in classifiers},
    }

    # StratifiedKFold guarantees class ratio is preserved in every fold,
    # unlike repeated train_test_split which can produce imbalanced splits
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)

    for fold, (train_idx, test_idx) in enumerate(skf.split(X, y)):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        # Scale within each fold to avoid data leakage
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train)
        X_test = scaler.transform(X_test)

        for name, model in classifiers.items():
            model.fit(X_train, y_train)
            y_pred = model.predict(X_test)

            metrics["accuracy"][name].append(accuracy_score(y_test, y_pred))
            metrics["precision"][name].append(precision_score(y_test, y_pred, zero_division=0))
            metrics["recall"][name].append(recall_score(y_test, y_pred, zero_division=0))
            metrics["f1"][name].append(f1_score(y_test, y_pred, zero_division=0))

            y_scores = model.predict_proba(X_test)[:, 1]
            metrics["auc"][name].append(roc_auc_score(y_test, y_scores))

        print(f"Fold {fold + 1}/{n_splits} done")

    return metrics


def summarize_metrics(metrics):
    mean = {
        m: {k: np.mean(v) for k, v in metrics[m].items()}
        for m in metrics
    }

    std = {
        m: {k: np.std(v) for k, v in metrics[m].items()}
        for m in metrics
    }

    return mean, std


def save_results(mean_metrics, std_metrics, output_path):
    rows = []

    for model_name in mean_metrics["accuracy"]:
        row = {"model": model_name}
        for metric in mean_metrics:
            row[f"{metric}_mean"] = round(mean_metrics[metric][model_name], 4)
            row[f"{metric}_std"] = round(std_metrics[metric][model_name], 4)
        rows.append(row)

    pd.DataFrame(rows).to_csv(output_path, index=False)
    print(f"Results saved to {output_path}")


# -----------------------
# Plotting
# -----------------------
def plot_metrics(mean_metrics, std_metrics):
    names = list(mean_metrics["accuracy"].keys())
    _, axes = plt.subplots(1, 2, figsize=(14, 5))

    for ax, metric in zip(axes, ["accuracy", "auc"]):
        means = [mean_metrics[metric][n] for n in names]
        stds = [std_metrics[metric][n] for n in names]

        ax.bar(names, means, yerr=stds, capsize=5)
        ax.set_ylabel(f"Mean {metric.upper()}")
        ax.set_title(f"Model Comparison ({metric.upper()})")
        ax.set_xticklabels(names, rotation=45, ha="right")
        ax.set_ylim(0, 1)

    plt.tight_layout()
    plt.show()


# -----------------------
# Main
# -----------------------
def main():
    X, y = prepare_data(TRUE_PATH, FALSE_PATH)

    metrics = evaluate_models(X, y, N_SPLITS)
    mean_metrics, std_metrics = summarize_metrics(metrics)

    for metric in mean_metrics:
        print(f"\n{metric.upper()}")
        for model in mean_metrics[metric]:
            mean_val = mean_metrics[metric][model]
            std_val = std_metrics[metric][model]
            print(f"  {model}: {mean_val:.4f} ± {std_val:.4f}")

    save_results(mean_metrics, std_metrics, RESULTS_PATH)
    plot_metrics(mean_metrics, std_metrics)


if __name__ == "__main__":
    main()
