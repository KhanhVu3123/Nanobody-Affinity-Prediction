#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Extract specific chains (nanobody + antigen) from PDB files
and save them as gzipped PDB files.
"""

import gzip
import os

import ProcessPdbFile


# -----------------------
# Core function
# -----------------------
def write_filtered_pdb(pdb_id, nanobody_chain, antigen_chain, base_dir, output_dir):
    """Write selected chains from a PDB file into a gzipped file."""
    
    input_path = os.path.join(base_dir, f"{pdb_id}.pdb")
    output_path = os.path.join(
        output_dir, f"{pdb_id}_{nanobody_chain}{antigen_chain}.pdb.gz"
    )

    with open(input_path, "r") as infile, gzip.open(output_path, "wt") as outfile:
        for line in infile:
            record_type = ProcessPdbFile.clean(line[:6])

            if record_type != "ATOM":
                continue

            chain = line[21]

            if chain in (nanobody_chain, antigen_chain):
                outfile.write(line)


# -----------------------
# Helper
# -----------------------
def load_chain_mapping(file_path):
    """
    Reads mapping file and returns:
    {pdb_id: (nanobody_chain, antigen_chain)}
    """
    mapping = {}

    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()

            pdb_id = line[:4]
            nanobody_chain = line[22]
            antigen_chain = line[35]

            mapping[pdb_id] = (nanobody_chain, antigen_chain)

    return mapping


# -----------------------
# Main
# -----------------------
def main():
    mapping_file = "NearestChain.txt"
    base_dir = "data/PDB"
    output_dir = "data/PDB"

    chain_map = load_chain_mapping(mapping_file)

    print(f"Processing {len(chain_map)} PDB files...")

    for i, (pdb_id, (nano_chain, anti_chain)) in enumerate(chain_map.items()):
        print(f"{i+1}/{len(chain_map)}: {pdb_id}")

        write_filtered_pdb(
            pdb_id,
            nano_chain,
            anti_chain,
            base_dir,
            output_dir
        )

    print("Done.")


if __name__ == "__main__":
    main()