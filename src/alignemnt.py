#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys, os
import subprocess
from pathlib import Path


def align_with_mafft(input_fasta: str) -> None:
    """
    Aligns sequences in the input FASTA file using MAFFT and writes the aligned sequences to the output FASTA file.
    
    Args:
        input_fasta (str): Path to the input FASTA file containing unaligned sequences.
    """
    
    file_path = Path(input_fasta)
    if not file_path.exists():
        raise FileNotFoundError(f"Input FASTA file '{input_fasta}' does not exist.")
    folder_path = file_path.parent
    output_fasta = folder_path / f"{file_path.stem}_aligned.fasta"
    
    cmd = ["mafft", "--auto", input_fasta]

    # MAFFT writes the alignment to stdout, so capture it and save to a file
    with open(output_fasta, "w") as out:
        subprocess.run(cmd, stdout=out, check=True)

    return output_fasta


def main():
    parser = argparse.ArgumentParser(description="Identify regions of conservation in a set of sequences.")
    parser.add_argument("-f", "--fasta", type=str, help="FASTA file containing sequences to analyze.")
    parser.add_argument("-o", "--out_file", type=str, help="File where the downloaded sequences are written.")
    args = parser.parse_args()

    # Align the sequences using MAFFT
    aligned_fasta = align_with_mafft(args.fasta)

if __name__ == "__main__":
    main()