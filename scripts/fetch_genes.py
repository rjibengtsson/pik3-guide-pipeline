#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from fetch import download_gene_id

def main():
    parser = argparse.ArgumentParser(description="Download RNA and protein sequences for an NCBI Gene ID.")
    parser.add_argument("-id", "--gene_id", type=str, help="NCBI Gene ID (e.g. 672 for human BRCA1).")
    parser.add_argument("-o", "--out_file", type=str, help="File where the downloaded sequences are written.")
    args = parser.parse_args()

    download_gene_id(args.gene_id, args.out_file)

if __name__ == "__main__":
    main()