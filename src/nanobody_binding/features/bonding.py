#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Utilities for analyzing amino acid bonding in nanobody PDB structures.
"""

import os
from collections import Counter

import numpy as np
import pandas as pd

import ProcessPdbFile
import FindingCDRRegion


# -----------------------
# Helpers
# -----------------------
def add_percentage(df):
    """Convert counts to percentage strings."""
    total = df.to_numpy().sum()
    df = (df / total * 100).round(3)
    return df.astype(str) + "%"


def build_matrix_from_counts(pair_list):
    """Convert pair list into a DataFrame matrix."""
    count_dict = Counter(pair_list)

    rows = sorted({p[0] for p in pair_list})
    cols = sorted({p[1] for p in pair_list})

    df = pd.DataFrame(0, index=rows, columns=cols)

    for pair, count in count_dict.items():
        df.at[pair[0], pair[1]] = count

    return df


# -----------------------
# Core functions
# -----------------------
def find_bonding_pairs(pdb_file, chain, threshold):
    """
    Intermolecular bonding (nanobody vs antigen).
    """
    big_dict = ProcessPdbFile.get_sequence(pdb_file)

    if chain not in big_dict:
        raise ValueError(f"Chain {chain} not found in {pdb_file}")

    # Extract sequence
    aminolst = big_dict[chain]
    sequence = "".join(next(iter(d)) for d in aminolst)

    # Extract CDR regions
    cdr_idx = FindingCDRRegion.find_CDR(sequence)
    aminolst = (
        aminolst[cdr_idx[0]:cdr_idx[1]+1] +
        aminolst[cdr_idx[2]:cdr_idx[3]+1] +
        aminolst[cdr_idx[4]:cdr_idx[5]+1]
    )

    bonding = []

    for other_chain, aminolst2 in big_dict.items():
        if other_chain == chain:
            continue

        for aa1_dict in aminolst:
            aa1 = next(iter(aa1_dict))
            coord1 = ProcessPdbFile.to_float(aa1_dict[aa1])

            for aa2_dict in aminolst2:
                aa2 = next(iter(aa2_dict))
                coord2 = ProcessPdbFile.to_float(aa2_dict[aa2])

                dist_sq = sum((coord1[i] - coord2[i]) ** 2 for i in range(3))

                if dist_sq < threshold ** 2:
                    bonding.append(aa1 + aa2)

    return bonding


def find_self_pairs(pdb_file, threshold):
    """
    Intramolecular bonding (within same chain).
    """
    big_dict = ProcessPdbFile.get_sequence(pdb_file)
    pairs = []

    for chain_data in big_dict.values():
        for aa1_dict in chain_data:
            aa1 = next(iter(aa1_dict))
            coord1 = aa1_dict[aa1]

            for aa2_dict in chain_data:
                aa2 = next(iter(aa2_dict))
                coord2 = aa2_dict[aa2]

                dist = ProcessPdbFile.coordinate_distance(coord1, coord2)

                if 0 < dist <= threshold:
                    pairs.append(aa1 + aa2)

    return pairs


# -----------------------
# Cleaning
# -----------------------
def clean_self_pairs(pair_list):
    """Remove symmetric duplicates (AB == BA)."""
    seen = set()
    cleaned = []

    for p in pair_list:
        if len(p) != 2:
            continue

        if p[::-1] not in seen:
            seen.add(p)
            cleaned.append(p)

    return cleaned


# -----------------------
# Table construction
# -----------------------
def construct_self_bonding_table(base_dir, nanobody_file, threshold=8):
    names = [line[:4] for line in open(nanobody_file)]

    all_pairs = []

    for i, name in enumerate(names):
        print(f"{i+1}/{len(names)}")

        pdb_path = os.path.join(base_dir, f"{name}.pdb")
        all_pairs.extend(find_self_pairs(pdb_path, threshold))

    all_pairs = clean_self_pairs(all_pairs)

    df = build_matrix_from_counts(all_pairs)
    df = add_percentage(df)

    df.to_csv("data/sequence/self_bonding.csv")
    return df


def construct_inter_bonding_table(base_dir, nanobody_file, threshold=8):
    entries = [line.strip() for line in open(nanobody_file)]

    all_pairs = []

    for i, line in enumerate(entries):
        pdb_id = line[:4]
        chain = line[22]

        print(f"{i+1}/{len(entries)}")

        pdb_path = os.path.join(base_dir, f"{pdb_id}.pdb")
        all_pairs.extend(find_bonding_pairs(pdb_path, chain, threshold))

    df = build_matrix_from_counts(all_pairs)
    df = add_percentage(df)

    df.to_csv("data/sequence/intermolecular_bonding.csv")
    return df


# -----------------------
# Main
# -----------------------
def main():
    df = pd.read_csv("data/sequence/intermolecular_bonding.csv", index_col=0)
    df = add_percentage(df)
    df.to_csv("data/sequence/intermolecular_bonding.csv")


if __name__ == "__main__":
    main()