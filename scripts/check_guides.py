#!/usr/bin/env python3
"""Flag guides on homopolymer runs and on GC outside an acceptable band.

    python scripts/check_guides.py \
        --guides results/pik3/guides.fasta \
        --output results/pik3/guide_qc.tsv

Reads a guide FASTA -- the one stage 1 writes, headers set to ``tile_id`` --
and writes one row per guide carrying its GC fraction, its homopolymer runs,
and the flags those raise:

``gc_low`` / ``gc_high``
    GC outside the band, by default 40-60%. The two are separate flags rather
    than one ``gc_out_of_band`` because they have opposite remedies, and the
    split is what tells a badly set band from a GC-skewed transcript.
``homopolymer``
    A run of 4 or more identical bases anywhere in the guide -- GGGG, UUUU,
    AAAA, CCCC. All four bases count: unlike the insert-scale check in
    check_target_sites.py there is no separate poly-T scan to hand T off to,
    because a run inside a spacer is a property of that spacer rather than a
    U6 termination signal in the construct.

A guide with no flags has an empty ``flags`` cell and ``n_flags`` 0. There is
no separate pass/fail column: the verdict is derivable from the row, and both
thresholds are echoed into every row (``gc_band``, ``homopolymer_min``) so a
flag can be audited without knowing how this was invoked. ``--flagged-only``
drops the passing rows if the flagged ones are all that is wanted; the default
writes every guide, so the file can be joined back onto ``tiles.tsv`` on
``tile_id``.

Both alphabets are accepted. GC counts U with T, and a run is named in
whichever alphabet the input used, so a run of four U reports base ``U`` and
detail ``UUUU`` rather than the ``T`` the scanner works in internally. That
matters because ``guides.fasta`` is RNA while ``--guide-alphabet dna`` output
is not, and a detail string naming a base the guide does not contain would be
unsearchable.

Offsets are 1-based inclusive and relative to the guide, so ``guide[start - 1:
end]`` is the run. Guide coordinates are not transcript coordinates: the guide
is the reverse complement of its target, so offset 1 in the guide is the
*last* base of the target span named in its ``tile_id``.

Bases outside ACGTU are excluded from the GC denominator and counted in
``n_ambiguous`` rather than scored as not-GC. They are reported, not flagged.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from typing import Sequence

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

import seqqc  # noqa: E402
import tile as tile_mod  # noqa: E402

COLUMNS = (
    "guide_id", "length", "sequence",
    "gc_frac", "gc_percent", "gc_band",
    "n_ambiguous",
    "n_homopolymer_runs", "longest_homopolymer", "longest_homopolymer_base",
    "longest_homopolymer_at", "homopolymer_min", "homopolymer_runs",
    "flags", "n_flags",
)

FLAG_NAMES = ("gc_low", "gc_high", "homopolymer")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--guides",
        default=os.path.join("results", "pik3", "guides.fasta"),
        metavar="PATH",
        help="guide FASTA, headers used as guide_id (default: %(default)s)",
    )
    parser.add_argument(
        "--output",
        default=os.path.join("results", "pik3", "guide_qc.tsv"),
        metavar="PATH",
        help="output table (default: %(default)s). Parent directory is "
        "created if absent.",
    )
    parser.add_argument(
        "--gc-min", type=float, default=seqqc.DEFAULT_GUIDE_GC_MIN, metavar="X",
        help="lower edge of the acceptable GC band for one guide, as a "
        "fraction. Inclusive, so a guide at exactly this GC passes "
        "(default: %(default)s)",
    )
    parser.add_argument(
        "--gc-max", type=float, default=seqqc.DEFAULT_GUIDE_GC_MAX, metavar="X",
        help="upper edge of the acceptable GC band, inclusive "
        "(default: %(default)s)",
    )
    parser.add_argument(
        "--homopolymer-min", type=int,
        default=seqqc.DEFAULT_GUIDE_HOMOPOLYMER_MIN, metavar="N",
        help="shortest run of identical bases that raises the homopolymer "
        "flag, counted over all four bases (default: %(default)s)",
    )
    parser.add_argument(
        "--flagged-only", action="store_true",
        help="write only the guides that raised at least one flag. The "
        "default writes every guide, which is what lets the output be joined "
        "onto tiles.tsv on tile_id.",
    )
    args = parser.parse_args(argv)
    if args.gc_min > args.gc_max:
        parser.error(
            f"--gc-min ({args.gc_min}) is above --gc-max ({args.gc_max}), so "
            "no GC can pass"
        )
    if args.homopolymer_min < 2:
        parser.error(
            f"--homopolymer-min must be at least 2, got {args.homopolymer_min}: "
            "a run of 1 is every base, so every guide would be flagged"
        )
    return args


def read_guides(path: str) -> list[tuple[str, str]]:
    """``(guide_id, sequence)`` in file order.

    Duplicate ids are a hard error, as they are in stage 1: the output is keyed
    on ``guide_id`` and a duplicate would make two different rows claim the
    same key, so a downstream join would silently pick one.
    """
    guides: list[tuple[str, str]] = []
    seen: dict[str, int] = {}
    for position, record in enumerate(tile_mod.read_fasta(path), start=1):
        if record.transcript_id in seen:
            raise ValueError(
                f"{path}: duplicate guide id {record.transcript_id!r} at "
                f"record {position}, first seen at record "
                f"{seen[record.transcript_id]}"
            )
        seen[record.transcript_id] = position
        guides.append((record.transcript_id, record.sequence))
    return guides


def row_for(result: seqqc.GuideQC) -> dict:
    """One output row. ``result`` carries its own thresholds, so this adds none."""
    longest = result.longest_run
    return {
        "guide_id": result.guide_id,
        "length": result.length,
        "sequence": result.sequence,
        "gc_frac": round(result.gc, 4),
        # Both forms, because the band is quoted in percent in every
        # conversation about it and in fractions in the code.
        "gc_percent": round(100 * result.gc, 2),
        "gc_band": f"{result.gc_min}-{result.gc_max}",
        "n_ambiguous": result.n_ambiguous,
        "n_homopolymer_runs": len(result.runs),
        "longest_homopolymer": longest.length if longest else 0,
        "longest_homopolymer_base": longest.base if longest else "",
        "longest_homopolymer_at": longest.start if longest else "",
        "homopolymer_min": result.homopolymer_min,
        # Every run, not just the longest: a guide with three separate 4-runs
        # is a different proposition from one with a single 4-run, and the
        # longest-only columns cannot tell them apart.
        "homopolymer_runs": ";".join(
            f"{run.base}:{run.start}-{run.end}(x{run.length})"
            for run in result.runs
        ),
        "flags": ";".join(result.flags),
        "n_flags": len(result.flags),
    }


def write_tsv(path: str, columns: Sequence[str], rows: list[dict]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as handle:
        # lineterminator="\n" to match every other table this pipeline writes;
        # csv defaults to "\r\n", which leaves a stray CR on the last column.
        writer = csv.DictWriter(
            handle,
            fieldnames=list(columns),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    if not os.path.isfile(args.guides):
        print(f"error: no such guide FASTA: {args.guides}", file=sys.stderr)
        return 1
    try:
        guides = read_guides(args.guides)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    if not guides:
        print(f"error: no FASTA records in {args.guides}", file=sys.stderr)
        return 1

    results = [
        seqqc.check_guide(
            guide_id,
            sequence,
            gc_min=args.gc_min,
            gc_max=args.gc_max,
            homopolymer_min=args.homopolymer_min,
        )
        for guide_id, sequence in guides
    ]

    rows = [
        row_for(result)
        for result in results
        if result.flags or not args.flagged_only
    ]

    parent = os.path.dirname(args.output)
    if parent:
        os.makedirs(parent, exist_ok=True)
    write_tsv(args.output, COLUMNS, rows)

    # Counts over every guide, not over the rows written: with
    # --flagged-only the two differ, and a summary that only described the
    # written subset could not state what fraction passed.
    flagged = [r for r in results if r.flags]
    print(f"guides read:  {len(results)} from {args.guides}", file=sys.stderr)
    print(f"GC band:      {args.gc_min}-{args.gc_max} (inclusive)",
          file=sys.stderr)
    print(f"homopolymer:  runs of >= {args.homopolymer_min} identical bases",
          file=sys.stderr)
    print(f"passing:      {len(results) - len(flagged)}", file=sys.stderr)
    print(f"flagged:      {len(flagged)}", file=sys.stderr)
    for name in FLAG_NAMES:
        # Guides can raise more than one flag, so these do not sum to the
        # flagged total.
        count = sum(1 for r in results if name in r.flags)
        print(f"  {name:<12}{count}", file=sys.stderr)
    ambiguous = sum(1 for r in results if r.n_ambiguous)
    if ambiguous:
        print(
            f"note:         {ambiguous} guide(s) contain bases outside ACGTU; "
            "excluded from the GC denominator, not flagged",
            file=sys.stderr,
        )
    print(f"rows written: {len(rows)} -> {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
