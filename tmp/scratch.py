import os, sys
from typing import Iterable, Iterator, NamedTuple
import subprocess
from Bio import SeqIO
from Bio.Data import IUPACData
from Bio.Seq import (
    back_transcribe,
    reverse_complement,
    reverse_complement_rna,
    transcribe,
)



def _run_command(command: list[str], capture: bool = False):
    # Not `run`: that would shadow TestCase.run and unittest would call it
    # with a TestResult.
    completed = subprocess.run(
        command,
        stdout=subprocess.PIPE if capture else subprocess.DEVNULL,
        stderr=subprocess.PIPE,
    )
    if completed.returncode != 0:
        raise AssertionError(
            f"{command[0]} failed (exit {completed.returncode}): "
            f"{completed.stderr.decode('utf-8', 'replace').strip()}"
        )
    return completed


def back_transcribe_fasta(in_file: str, out_file: str) -> None:

    with open(out_file, "w") as out_f:
        for record in SeqIO.parse(in_file, "fasta"):
            seq = record.seq
            target_seq = seq.back_transcribe()
            SeqRecord = SeqIO.SeqRecord(target_seq, id=record.id, description=record.description)
            SeqIO.write(SeqRecord, out_f, "fasta")




FIELDS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
          "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

def run_blastn(query_file: str, db_file: str, output_file: str) -> None:
    command = [
        "blastn",
        "-task", "blastn-short",
        "-query", query_file,
        "-db", db_file,
        "-dust", "no",
        "-soft_masking", "false",
        "-ungapped",
        "-perc_identity", "100",
        "-evalue", "1000",
        "-max_target_seqs", "100",
        "-outfmt", " ".join(["6", *FIELDS]),
    ]
    with open(output_file, "w", encoding="utf-8") as handle:
        handle.write("\t".join(FIELDS) + "\n")
        handle.flush()                      # before the child appends to the same fd
        completed = subprocess.run(command, stdout=handle, stderr=subprocess.PIPE)
    if completed.returncode != 0:
        raise AssertionError(f"blastn failed: {completed.stderr.decode().strip()}")


def get_length_of_fasta(fasta_file: str) -> int:

    for record in SeqIO.parse(fasta_file, "fasta"):
        print(f"Record ID: {record.id}, Length: {len(record.seq)}")




def main():

    in_file = sys.argv[1]
    # out_file = sys.argv[2]
    # db_file = sys.argv[3]

    # back_transcribe_fasta(in_file, out_file)
    # run_blastn(in_file, db_file, out_file)
    get_length_of_fasta(in_file)



if __name__ == "__main__":
    main()