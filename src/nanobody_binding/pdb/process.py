#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PDB processing utilities:
- Extract sequences from PDB / gzipped PDB
- Extract nanobody / antigen chains
- Compute distances between residues
"""

import gzip
import os

import FindingCDRRegion


# -----------------------
# Constants
# -----------------------
AMINO_ACID_MAP = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'UNK': ''
}


# -----------------------
# Helpers
# -----------------------
def clean(s: str) -> str:
    return s.strip()


def parse_atom_line(line):
    """Extract useful fields from a PDB line."""
    record = clean(line[:6])
    atom = clean(line[12:16])
    chain = line[21]
    res = clean(line[17:21])

    return record, atom, chain, res


def get_coordinates(line):
    return [
        float(clean(line[30:38])),
        float(clean(line[38:46])),
        float(clean(line[46:54]))
    ]


# -----------------------
# Sequence extraction
# -----------------------
def retrieve_sequence(file_path):
    """Extract sequences from PDB file (per chain)."""
    sequences = []
    seq = ""

    with open(file_path, "r") as f:
        for line in f:
            record, atom, _, res = parse_atom_line(line)

            if record == "HETATOM":
                break

            if record == "TER":
                if seq:
                    sequences.append(seq)
                    seq = ""
                continue

            if record != "ATOM" or atom != "CA":
                continue

            if res in AMINO_ACID_MAP:
                seq += AMINO_ACID_MAP[res]

    if seq:
        sequences.append(seq)

    return sequences


def retrieve_sequence_from_gz(file_path):
    """Extract sequences from gzipped PDB."""
    sequences = []
    seq = ""
    prev_chain = None

    with gzip.open(file_path, "rt") as f:
        for line in f:
            record, atom, chain, res = parse_atom_line(line)

            if record == "HETATOM":
                break

            if chain != prev_chain and seq:
                sequences.append(seq)
                seq = ""

            if record == "TER":
                if seq:
                    sequences.append(seq)
                    seq = ""
                continue

            if record != "ATOM" or atom != "CA":
                continue

            if res in AMINO_ACID_MAP:
                seq += AMINO_ACID_MAP[res]

            prev_chain = chain

    if seq:
        sequences.append(seq)

    return sequences


# -----------------------
# Nanobody / antigen extraction
# -----------------------
def extract_chain_sequence(gzip_file, chain_idx):
    """Generic function for extracting a specific chain."""
    file_name = os.path.basename(gzip_file)
    target_chain = file_name[chain_idx]

    seq = ""

    with gzip.open(gzip_file, "rt") as f:
        for line in f:
            record, atom, chain, res = parse_atom_line(line)

            if record == "HETATOM":
                break

            if record == "ATOM" and atom == "CA" and chain == target_chain:
                if res in AMINO_ACID_MAP:
                    seq += AMINO_ACID_MAP[res]

    return seq


def retrieve_nanobody_seq(gzip_file):
    return extract_chain_sequence(gzip_file, chain_idx=5)


def retrieve_antigen_seq(gzip_file):
    return extract_chain_sequence(gzip_file, chain_idx=6)


def retrieve_all_sequences(gzip_folder, mode="nanobody"):
    """Batch extraction."""
    files = [f for f in os.listdir(gzip_folder) if f.endswith(".pdb.gz")]
    result = {}

    for f in files:
        path = os.path.join(gzip_folder, f)
        pdb_id = f[:4]

        if mode == "nanobody":
            seq = retrieve_nanobody_seq(path)
        elif mode == "antigen":
            seq = retrieve_antigen_seq(path)
        else:
            raise ValueError("mode must be 'nanobody' or 'antigen'")

        result[pdb_id] = seq

    return result


# -----------------------
# Coordinate + distance
# -----------------------
def coordinate_distance(c1, c2):
    """Squared Euclidean distance."""
    return sum((float(c1[i]) - float(c2[i])) ** 2 for i in range(3))


def get_sequence_with_coords(file_path):
    """
    Return:
    {chain: [{AA: [x,y,z]}, ...]}
    """
    result = {}
    current_chain = None
    chain_data = []

    with open(file_path, "r") as f:
        for line in f:
            record, atom, chain, res = parse_atom_line(line)

            if record == "HETATOM":
                break

            if record == "TER":
                if current_chain is not None:
                    result[current_chain] = chain_data
                chain_data = []
                continue

            if record != "ATOM" or atom != "CA":
                continue

            if res not in AMINO_ACID_MAP:
                continue

            coords = get_coordinates(line)
            aa = AMINO_ACID_MAP[res]

            chain_data.append({aa: coords})
            current_chain = chain

    if current_chain and chain_data:
        result[current_chain] = chain_data

    return result


# -----------------------
# Chain interaction
# -----------------------
def determine_chain(chain, file_path, threshold):
    big_dict = get_sequence_with_coords(file_path)

    if chain not in big_dict:
        raise ValueError(f"Chain {chain} not found")

    aminolst = big_dict[chain]
    sequence = "".join(next(iter(d)) for d in aminolst)

    cdr = FindingCDRRegion.find_CDR(sequence)
    aminolst = (
        aminolst[cdr[0]:cdr[1]+1] +
        aminolst[cdr[2]:cdr[3]+1] +
        aminolst[cdr[4]:cdr[5]+1]
    )

    results = {}

    for other_chain, other_list in big_dict.items():
        if other_chain == chain:
            continue

        count = 0

        for aa1_dict in aminolst:
            coord1 = next(iter(aa1_dict.values()))

            for aa2_dict in other_list:
                coord2 = next(iter(aa2_dict.values()))

                if coordinate_distance(coord1, coord2) < threshold ** 2:
                    count += 1

        results[other_chain] = count

    return results


# -----------------------
# Main
# -----------------------
def main():
    sequences = retrieve_sequence("data/PDB/4LSP.pdb")
    print(sequences)


if __name__ == "__main__":
    main()