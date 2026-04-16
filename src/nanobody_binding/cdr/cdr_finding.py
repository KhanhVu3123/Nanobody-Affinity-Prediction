#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CDR detection using alignment against nanobody framework regions.
"""

from Bio import pairwise2
from Bio.Seq import Seq


# -----------------------
# Utilities
# -----------------------
def edit_distance(x, y):
    """Compute Levenshtein distance."""
    m, n = len(x), len(y)
    D = [[0] * (n + 1) for _ in range(m + 1)]

    for i in range(m + 1):
        D[i][0] = i
    for j in range(n + 1):
        D[0][j] = j

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            cost = 0 if x[i - 1] == y[j - 1] else 1
            D[i][j] = min(
                D[i - 1][j - 1] + cost,  # substitution
                D[i - 1][j] + 1,         # deletion
                D[i][j - 1] + 1          # insertion
            )

    return D[m][n]


def remove_dash(s):
    return s.replace("-", "")


def split_blocks(s):
    """Split string by '-' and remove empty parts."""
    return [block for block in s.split("-") if block]


# -----------------------
# Alignment helper
# -----------------------
def select_best_alignment(seq, framework_seq):
    """Pick alignment with minimum fragmentation."""
    alignments = pairwise2.align.globalxx(seq, framework_seq)

    if not alignments:
        return None

    min_blocks = min(len(split_blocks(a[1])) for a in alignments)

    for align in alignments:
        if len(split_blocks(align[1])) == min_blocks:
            return align

    return None


# -----------------------
# CDR detection
# -----------------------
def find_CDR(sequence):
    """
    Return indices of 3 CDR regions:
    [start1, end1, start2, end2, start3, end3]
    """

    # Framework regions
    FR = [
        "QVQLVESGGGLVQAGGSLRLSCAASG",
        "WYRQAPGKQRELVA",
        "DSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYC",
        "WGQGTQVTVSS"
    ]

    framework_concat = "".join(FR)
    sequence = Seq(sequence)

    align = select_best_alignment(sequence, framework_concat)
    if align is None:
        raise ValueError("Alignment failed. Check input sequence.")

    seqA, seqB = align[0], align[1]

    index_list = []

    # Locate framework regions
    for fr in FR:
        best_idx = 0
        best_score = float("inf")

        window_len = len(fr) + 5

        # Sliding window search
        for i in range(len(seqB) - window_len + 1):
            kmer = remove_dash(seqB[i:i + window_len])
            score = edit_distance(kmer, fr)

            if score < best_score:
                best_score = score
                best_idx = i

        # Fine-tune alignment
        best_extension = 0
        best_score = float("inf")

        for j in range(min(6, len(seqB) - best_idx - len(fr))):
            segment = remove_dash(seqB[best_idx:best_idx + len(fr) + j])
            score = edit_distance(segment, fr)

            if score < best_score:
                best_score = score
                best_extension = j

        # Adjust for gaps
        start = best_idx
        while seqB[start] == "-":
            start += 1

        end = best_idx + len(fr) + best_extension
        while seqB[end] == "-":
            end -= 1

        index_list.extend([start, end])

    # Extract CDR regions from aligned sequence
    cdrs = [
        seqA[index_list[1] + 1:index_list[2]],
        seqA[index_list[3] + 1:index_list[4]],
        seqA[index_list[5] + 1:index_list[6]],
    ]

    # Map back to original sequence
    clean_seq = str(sequence)
    result = []

    for cdr in cdrs:
        cdr_clean = remove_dash(cdr)
        start = clean_seq.index(cdr_clean)
        end = start + len(cdr_clean)
        result.extend([start, end])

    return result


# -----------------------
# Main
# -----------------------
def main():
    seq = input("Enter sequence: ").strip()
    indices = find_CDR(seq)
    print(indices)


if __name__ == "__main__":
    main()