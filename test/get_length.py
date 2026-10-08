import os, sys
from Bio import SeqIO


def main():
    path = sys.argv[1]
    for record in SeqIO.parse(path, "fasta"):
        print(record.id, len(record.seq))

if __name__ == "__main__":
    main()
