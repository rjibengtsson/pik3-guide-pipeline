#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys, os
import subprocess

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from conservation import align_with_mafft



def main():
    parser = argparse.ArgumentParser(description="Identify regions of conservation in a set of sequences.")
    parser.add_argument("-f", "--fasta", type=str, help="FASTA file containing sequences to analyze.")
    # parser.add_argument("-o", "--out_file", type=str, help="File where the downloaded sequences are written.")
    args = parser.parse_args()

    # Align the sequences using MAFFT
    # aligned_fasta = align_with_mafft(args.fasta)


    aligned_fasta = ""

if __name__ == "__main__":
    main()