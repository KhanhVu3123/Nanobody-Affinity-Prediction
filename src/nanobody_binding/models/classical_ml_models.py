#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compare multiple ML algorithms on nanobody–antigen embedding data.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
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

    scaler = StandardScaler()
    X = scaler.fit_transform(X)

    return X, y


# -----------------------
# Models
# -----------------------
def get_classifiers():
    return {
        "Logistic Regression": LogisticRegression(max_iter=1000),
        "SVM": SVC(probability=True),
        "Random Forest": RandomForestClassifier(),
        "KNN": KNeighborsClassifier(),
        "Naive Bayes": GaussianNB(),
        "Decision Tree": DecisionTreeClassifier(),
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

    for i in range(n_splits):
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, random_state=42 + i
        )

        for name, model in classifiers.items():
            model.fit(X_train, y_train)
            y_pred = model.predict(X_test)

            metrics["accuracy"][name].append(accuracy_score(y_test, y_pred))
            metrics["precision"][name].append(precision_score(y_test, y_pred))
            metrics["recall"][name].append(recall_score(y_test, y_pred))
            metrics["f1"][name].append(f1_score(y_test, y_pred))

            y_scores = model.predict_proba(X_test)[:, 1]
            metrics["auc"][name].append(roc_auc_score(y_test, y_scores))

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


# -----------------------
# Plotting
# -----------------------
def plot_accuracy(mean_metrics, std_metrics):
    names = list(mean_metrics["accuracy"].keys())
    means = list(mean_metrics["accuracy"].values())
    stds = list(std_metrics["accuracy"].values())

    plt.figure(figsize=(10, 5))
    plt.bar(names, means, yerr=stds, capsize=5)

    plt.ylabel("Mean Accuracy")
    plt.title("Model Comparison (Accuracy)")
    plt.xticks(rotation=45)

    plt.tight_layout()
    plt.show()


# -----------------------
# Main
# -----------------------
def main():
    X, y = prepare_data(TRUE_PATH, FALSE_PATH)

    metrics = evaluate_models(X, y, N_SPLITS)
    mean_metrics, std_metrics = summarize_metrics(metrics)

    # Print results
    for metric in mean_metrics:
        print(f"\n{metric.upper()}")
        for model in mean_metrics[metric]:
            mean_val = mean_metrics[metric][model]
            std_val = std_metrics[metric][model]
            print(f"{model}: {mean_val:.4f} ± {std_val:.4f}")

    plot_accuracy(mean_metrics, std_metrics)


if __name__ == "__main__":
    main()