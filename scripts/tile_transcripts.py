#!/usr/bin/env python3
"""Tile every input sequence in a folder into fixed-length k-mers and derive guides.

    python scripts/tile_transcripts.py --input-dir data/pik3-iso \
        --output results/pik3/tiles.tsv

Stages:

1. Take the molecule type from ``--input-type dna|mrna`` (required).
2. Tile each transcript into k-mers, one per starting position. ``target_seq``
   holds each tile verbatim, in the input's own alphabet.
3. Transcribe to RNA (T -> U) and reverse complement to give the Cas13 guide.
   DNA input needs the transcription step; mRNA input written with U skips it.

Target and guide are both written 5' -> 3'. The guide is antiparallel to its
target, so the guide's 5' end pairs with the target's 3' end.

Each output row is one tile, tagged with the transcript it came from and its
1-based inclusive start/end on that transcript. ``tile_id`` is the join key for
downstream stages (GC/MFE filtering, off-target mapping, ...).
"""

from __future__ import annotations

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from tile import (  # noqa: E402
    Record,
    Tile,
    back_transcribe,
    find_fasta_files,
    guide_rna,
    read_fasta,
    tile_record,
)

COLUMNS = [
    "tile_id",
    "transcript_id",
    "source_file",
    "start",
    "end",
    "strand",
    "length",
    "target_seq",
    "guide_rna",
    "has_ambiguous",
]

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input-dir",
        default="data/pik3-iso",
        help="folder of input FASTA files (default: %(default)s). Inputs are "
        "organised per dataset under data/; pass another folder for another set.",
    )
    parser.add_argument(
        "--output",
        default="results/pik3/tiles.tsv",
        help="combined TSV of all tiles (default: %(default)s). Results are "
        "organised per dataset under results/; point this at the matching folder "
        "when tiling another set, or the default will overwrite the PI3K run.",
    )
    parser.add_argument(
        "--input-type",
        choices=("dna", "mrna"),
        required=True,
        help="molecule type of the input sequences: 'dna' (coding/sense strand) "
        "or 'mrna'. Taken on trust and recorded in the run summary -- the "
        "sequences themselves are never inspected.",
    )
    parser.add_argument("-k", type=int, default=30, help="tile length (default: %(default)s)")
    parser.add_argument(
        "--step", type=int, default=1, help="bases to advance between tiles (default: %(default)s)"
    )
    parser.add_argument(
        "--guide-alphabet",
        choices=("rna", "dna"),
        default="rna",
        help="alphabet for the guide_rna column. 'rna' (default) is the real "
        "molecule; 'dna' writes T instead of U, which is what aligners like "
        "bowtie and blast expect for off-target mapping. target_seq is always "
        "verbatim from the input and is unaffected.",
    )
    parser.add_argument(
        "--per-transcript",
        action="store_true",
        help="also write {transcript_id}.tsv per transcript, into a tiles/ folder "
        "beside --output",
    )
    parser.add_argument(
        "--fasta-out",
        metavar="PATH",
        help="also write the guide sequences as FASTA, headers set to tile_id",
    )
    parser.add_argument(
        "--fasta-field",
        choices=("guide_rna", "target_seq"),
        default="guide_rna",
        help="which sequence --fasta-out writes (default: %(default)s)",
    )
    parser.add_argument(
        "--drop-ambiguous",
        action="store_true",
        help="omit tiles containing non-ACGTU bases (default: keep them, flagged in "
        "the has_ambiguous column)",
    )
    return parser.parse_args(argv)


def tile_row(tile: Tile, *, guide_alphabet: str) -> dict[str, object]:
    """One output row: the target tile and the guide that pairs with it.

    ``target_seq`` is the input sequence verbatim, in whatever alphabet the
    FASTA used. ``guide_rna`` is transcribed to RNA so it base-pairs correctly.
    Both are written 5' -> 3'.
    """
    target = tile.sequence  # verbatim from the input
    guide = guide_rna(target)
    if guide_alphabet == "dna":
        guide = back_transcribe(guide)
    return {
        "tile_id": tile.tile_id,
        "transcript_id": tile.transcript_id,
        "source_file": tile.source_file,
        "start": tile.start,
        "end": tile.end,
        "strand": "+",
        "length": tile.length,
        "target_seq": target,
        "guide_rna": guide,
        "has_ambiguous": int(tile.has_ambiguous),
    }


def writer_for(path: str, columns: list[str]) -> tuple[object, csv.DictWriter]:
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
    handle = open(path, "w", newline="", encoding="utf-8")
    writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    return handle, writer


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    if args.k < 1:
        sys.exit(f"error: -k must be >= 1, got {args.k}")
    if args.step < 1:
        sys.exit(f"error: --step must be >= 1, got {args.step}")

    try:
        paths = find_fasta_files(args.input_dir)
    except NotADirectoryError as exc:
        sys.exit(f"error: {exc}")
    if not paths:
        sys.exit(f"error: no FASTA files (*.fasta, *.fa, *.fna) in {args.input_dir}")

    combined_handle, combined = writer_for(args.output, COLUMNS)
    fasta_handle = None
    if args.fasta_out:
        os.makedirs(os.path.dirname(os.path.abspath(args.fasta_out)), exist_ok=True)
        fasta_handle = open(args.fasta_out, "w", encoding="utf-8")

    per_transcript_dir = os.path.join(os.path.dirname(os.path.abspath(args.output)), "tiles")
    total_written = 0
    total_ambiguous = 0
    seen_ids: dict[str, str] = {}
    summary: list[tuple[str, int, int, int]] = []

    try:
        for path in paths:
            for record in read_fasta(path):
                if not record.sequence:
                    sys.exit(f"error: {record.source_file}: record {record.transcript_id} has no sequence")
                if len(record.sequence) < args.k:
                    sys.exit(
                        f"error: {record.source_file}: {record.transcript_id} is "
                        f"{len(record.sequence)} nt, shorter than k={args.k}"
                    )
                if record.transcript_id in seen_ids:
                    sys.exit(
                        f"error: transcript id {record.transcript_id} appears in both "
                        f"{seen_ids[record.transcript_id]} and {record.source_file}; "
                        "tile_ids would collide"
                    )
                seen_ids[record.transcript_id] = record.source_file

                written, ambiguous = write_record(
                    record,
                    args,
                    combined,
                    fasta_handle,
                    per_transcript_dir,
                )
                total_written += written
                total_ambiguous += ambiguous
                summary.append((record.transcript_id, len(record.sequence), written, ambiguous))
    finally:
        combined_handle.close()
        if fasta_handle is not None:
            fasta_handle.close()

    report(summary, args, total_written, total_ambiguous)
    return 0


def write_record(
    record: Record,
    args: argparse.Namespace,
    combined: csv.DictWriter,
    fasta_handle,
    per_transcript_dir: str,
) -> tuple[int, int]:
    """Write every tile of one record; return (written, ambiguous) counts."""
    per_handle = per_writer = None
    if args.per_transcript:
        safe = record.transcript_id.replace(os.sep, "_")
        per_handle, per_writer = writer_for(
            os.path.join(per_transcript_dir, f"{safe}.tsv"), COLUMNS
        )

    written = ambiguous = 0
    try:
        for tile in tile_record(record, k=args.k, step=args.step):
            if tile.has_ambiguous:
                ambiguous += 1
                if args.drop_ambiguous:
                    continue
            row = tile_row(tile, guide_alphabet=args.guide_alphabet)
            combined.writerow(row)
            if per_writer is not None:
                per_writer.writerow(row)
            if fasta_handle is not None:
                fasta_handle.write(f">{row['tile_id']}\n{row[args.fasta_field]}\n")
            written += 1
    finally:
        if per_handle is not None:
            per_handle.close()

    return written, ambiguous


def report(
    summary: list[tuple[str, int, int, int]],
    args: argparse.Namespace,
    total_written: int,
    total_ambiguous: int,
) -> None:
    step_taken = (
        "DNA input; transcribed to RNA when building guides"
        if args.input_type == "dna"
        else "mRNA input"
    )
    ambiguous_policy = "omitted" if args.drop_ambiguous else "kept and flagged"
    print(f"input type: {args.input_type} (as declared) - {step_taken}", file=sys.stderr)
    print(
        f"k={args.k}  step={args.step}  guide alphabet: {args.guide_alphabet}  "
        f"ambiguous tiles: {ambiguous_policy}",
        file=sys.stderr,
    )
    print(
        "target_seq is verbatim from the input; guide_rna is RNA. "
        "Both are written 5' -> 3'",
        file=sys.stderr,
    )
    print(f"{'transcript':<16}{'length':>10}{'tiles':>10}{'ambiguous':>12}", file=sys.stderr)
    for transcript_id, length, written, ambiguous in summary:
        print(f"{transcript_id:<16}{length:>10}{written:>10}{ambiguous:>12}", file=sys.stderr)
    print(f"{'TOTAL':<16}{'':>10}{total_written:>10}{total_ambiguous:>12}", file=sys.stderr)
    print(f"wrote {args.output}", file=sys.stderr)
    if args.per_transcript:
        per_dir = os.path.join(os.path.dirname(os.path.abspath(args.output)), "tiles")
        print(f"wrote per-transcript TSVs under {per_dir}", file=sys.stderr)
    if args.fasta_out:
        print(f"wrote {args.fasta_out} ({args.fasta_field})", file=sys.stderr)


if __name__ == "__main__":
    raise SystemExit(main())
