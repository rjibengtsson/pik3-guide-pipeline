#!/usr/bin/env python3
"""Map alignment coordinates onto the reference sequence's own coordinates.

Reads a FASTA file containing a multiple sequence alignment, takes one or more
regions given as **alignment** coordinates, and reports where each region lands
on the reference sequence once its gaps are removed.

The reference sequence must be present in the alignment, named by ``--reference``
(or, with ``--bed``, by each BED record's first column). Its gapped form is the
query and its ungapped form is the BLAST database, so the hit coordinates are
positions in the ungapped reference.

Regions come in one of two ways:

    # one region, alignment coordinates, 1-based inclusive
    python get_conserved_regions.py -f alignment.fasta -r NM_006218.4 \
        -i 1200,1320 -o conserved.tsv

    # many regions, from a BED file (0-based half-open, as BED always is)
    python get_conserved_regions.py -f alignment.fasta -b regions.bed \
        -o conserved.tsv

Output is a tab-separated file, one row per region: the region as given, then
the best blastn hit for it (highest bitscore), with `ref_start`/`ref_end`
1-based inclusive and always low -> high. `n_hits` is how many HSPs blastn
reported for that region, so a region sitting in a repeat is visible rather
than silently collapsed to one position. A region with no hit is written with
empty hit columns.

A relative ``--out_file`` is written under ``results/blastn/``.
"""
from __future__ import annotations

import argparse
import os
import sys
import tempfile

import pandas as pd
from pathlib import Path
from Bio import SeqIO

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from blast_utils import (
    DEFAULT_EVALUE,
    DEFAULT_MAX_TARGET_SEQS,
    DEFAULT_TASK,
    DEFAULT_WORD_SIZE,
    BlastError,
    best_hit_per_query,
    count_hits_per_query,
    make_blast_db,
    run_blastn,
    write_fasta,
)

ROOT_DIR = Path(__file__).resolve().parents[1]
BLAST_RESULTS_DIR = ROOT_DIR / "results" / "blastn"

# megablast's 28 nt word cannot find a region shorter than itself, so short
# queries switch to the same settings the stage-2 blast cross-check uses.
SHORT_QUERY_LEN = 30
SHORT_TASK = "blastn-short"
SHORT_WORD_SIZE = 7

HIT_COLUMNS = (
    "ref_start",
    "ref_end",
    "strand",
    "pident",
    "aln_length",
    "mismatch",
    "gapopen",
    "evalue",
    "bitscore",
    "n_hits",
)

INT_COLUMNS = ("ref_start", "ref_end", "aln_length", "mismatch", "gapopen", "n_hits")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Map alignment coordinates onto reference sequence coordinates."
    )
    parser.add_argument("-f", "--fasta", required=True,
                        help="Multiple sequence alignment FASTA file.")
    parser.add_argument("-o", "--out_file", required=True,
                        help="Output TSV. A relative path is written under results/blastn/.")
    # Option A: a BED file
    parser.add_argument("-b", "--bed",
                        help="BED file of regions to map (alternative to -r and -i). "
                             "Column 1 names the reference sequence in the alignment.")
    # Option B: reference accession and a single region
    parser.add_argument("-r", "--reference",
                        help="Reference sequence id in the alignment (required unless --bed is given).")
    parser.add_argument("-i", "--index",
                        help="Alignment positions of interest start,end (1-based inclusive). "
                             "Required unless --bed is given.")
    parser.add_argument("--threads", type=int, default=1,
                        help="blastn -num_threads (default: 1).")
    return parser


def load_aligned_sequences(fasta_file: str) -> dict[str, str]:
    """Read the alignment once, upper-cased, keyed by record id."""
    sequences: dict[str, str] = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id in sequences:
            raise ValueError(f"duplicate sequence id in {fasta_file}: {record.id!r}")
        sequences[record.id] = str(record.seq).upper()
    if not sequences:
        raise ValueError(f"no FASTA records found in {fasta_file}")
    return sequences


def map_regions(aligned_seq: str,
                ref_acc: str,
                regions: list[tuple[int, int]],
                tmpdir: str,
                threads: int = 1) -> list[dict]:
    """Blast each alignment region of *aligned_seq* against its ungapped self.

    *regions* are 1-based inclusive alignment coordinates. Returns one record
    per region, in the order given.
    """
    ungapped = aligned_seq.replace("-", "")

    os.makedirs(tmpdir, exist_ok=True)
    db_fasta = os.path.join(tmpdir, "reference.fasta")
    with open(db_fasta, "w") as fh:
        if write_fasta([(ref_acc, ungapped)], fh) == 0:
            raise ValueError(f"reference {ref_acc!r} is empty once gaps are removed")
    db = make_blast_db(db_fasta, os.path.join(tmpdir, "db"), force=False)

    # Regions are blasted one at a time: each needs its own task/word_size, and
    # a region can repeat, so hits must stay attributable to one region.
    records: list[dict] = []
    for start, end in regions:
        if not 1 <= start <= end <= len(aligned_seq):
            raise ValueError(
                f"region {start}-{end} is outside the alignment "
                f"(1-{len(aligned_seq)}) for {ref_acc!r}"
            )

        # Gaps are alignment artefacts and not nucleotides; blastn would treat
        # them as mismatches, so the query is the region's real bases only.
        query = aligned_seq[start - 1:end].replace("-", "")
        row = {
            "reference": ref_acc,
            "aln_start": start,
            "aln_end": end,
            "query_len": len(query),
        }

        if not query:
            print(f"warning: alignment region {start}-{end} of {ref_acc} is all gaps",
                  file=sys.stderr)
            records.append(row | {col: None for col in HIT_COLUMNS})
            continue

        short = len(query) < SHORT_QUERY_LEN
        hits = run_blastn(
            [(f"{ref_acc}:{start}-{end}", query)],
            db,
            task=SHORT_TASK if short else DEFAULT_TASK,
            word_size=SHORT_WORD_SIZE if short else DEFAULT_WORD_SIZE,
            evalue=DEFAULT_EVALUE,
            max_target_seqs=DEFAULT_MAX_TARGET_SEQS,
            threads=threads,
        )

        best = best_hit_per_query(hits)
        counts = count_hits_per_query(hits)
        if not best:
            print(f"warning: no blastn hit for {ref_acc} alignment region {start}-{end}",
                  file=sys.stderr)
            records.append(row | {col: None for col in HIT_COLUMNS})
            continue

        query_id, hit = next(iter(best.items()))
        hit_record = hit.as_record()
        row.update({col: hit_record[col] for col in HIT_COLUMNS if col != "n_hits"})
        row["n_hits"] = counts[query_id]
        records.append(row)

    return records


def read_bed(path: str) -> list[tuple[str, int, int]]:
    """BED is 0-based half-open; return 1-based inclusive regions."""
    bed_df = pd.read_csv(path, sep="\t", header=None,
                         usecols=[0, 1, 2], names=["chrom", "start", "end"])
    print(bed_df)
    return [(str(r.chrom), int(r.start) + 1, int(r.end)) for r in bed_df.itertuples()]


def resolve_out_path(out_file: str) -> Path:
    path = Path(out_file)
    if not path.is_absolute():
        path = BLAST_RESULTS_DIR / path
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    # Either --bed on its own, or both -r and -i
    if args.bed:
        if args.reference or args.index:
            parser.error("--bed cannot be combined with -r or -i")
    elif not (args.reference and args.index):
        parser.error("provide either --bed, or both -r and -i")

    if args.bed:
        regions = read_bed(args.bed)
    else:
        try:
            start, end = (int(v) for v in args.index.split(","))
        except ValueError:
            parser.error("-i/--index must be two integers, 'start,end'")
        regions = [(args.reference, start, end)]

    if not regions:
        parser.error("no regions to map")

    try:
        sequences = load_aligned_sequences(args.fasta)

        missing = sorted({acc for acc, _, _ in regions} - set(sequences))
        if missing:
            raise ValueError(
                "reference sequence(s) not found in "
                f"{args.fasta}: {', '.join(repr(m) for m in missing)}"
            )

        # Group by reference so each one's database is built exactly once.
        by_ref: dict[str, list[tuple[int, int]]] = {}
        for acc, start, end in regions:
            by_ref.setdefault(acc, []).append((start, end))

        records: list[dict] = []
        with tempfile.TemporaryDirectory(prefix="blast_") as tmpdir:
            for acc, acc_regions in by_ref.items():
                records.extend(
                    map_regions(sequences[acc], acc, acc_regions,
                                os.path.join(tmpdir, acc.replace(os.sep, "_")),
                                threads=args.threads)
                )
    except (BlastError, FileNotFoundError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    out_path = resolve_out_path(args.out_file)
    df = pd.DataFrame(records, columns=["reference", "aln_start", "aln_end", "query_len",
                                        *HIT_COLUMNS])
    # A region with no hit leaves NaN behind, which would otherwise float-ify
    # the coordinates ("247.0"); nullable Int64 keeps them as integers.
    for col in INT_COLUMNS:
        df[col] = df[col].astype("Int64")
    df.to_csv(out_path, index=False, sep="\t")

    mapped = int(df["ref_start"].notna().sum())
    print(f"{mapped}/{len(df)} region(s) mapped -> {out_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
