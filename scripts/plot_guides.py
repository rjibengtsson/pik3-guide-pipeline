#!/usr/bin/env python3
"""Plot where the surviving guides sit on one transcript, as interactive HTML.

    python scripts/plot_guides.py --fasta data/ppib/NM_000942.5.fasta \
        --predictions results/ppib/ppib-predictions.csv \
        --min-efficacy 0.95 --output results/ppib/guides.html

Stages:

1. Read the transcript FASTA and take its **accession** and **length** from the
   record itself (Biopython). One transcript per file, as under ``data/``.
2. Read the predictions CSV and keep only the rows whose ``transcript_id``
   matches that accession, so the guides plotted and the ruler they are measured
   against provably come from the same sequence.
3. Keep the rows with ``predicted_efficacy >= --min-efficacy``. The threshold is
   required: there is no default worth guessing, and a silent one would make two
   runs of this script incomparable.
4. Draw the kept guides along the full length of the transcript and write a
   standalone HTML file.

The plot is a view of one particular efficacy cut-off. It shows where the guides
are, not whether they are any good.

Stage 2's ``classification.tsv`` is **optional**: pass ``--classification`` and
each hover box also carries the guide's family classification, joined on
``tile_id``. Without it the plot is built from the predictions file alone.
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Iterable, Sequence

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

import pandas as pd  # noqa: E402

from tile import SeqIO  # noqa: E402  - re-exported so one import covers the stage
from visualize import guide_track, guides_from_rows, pack_lanes, write_html  # noqa: E402

EFFICACY_COLUMN = "predicted_efficacy"
REQUIRED_COLUMNS = ("tile_id", "transcript_id", "start", "end", EFFICACY_COLUMN)

#: Column names holding the guide (spacer) sequence, in the order they are tried.
#: A predictions file calls it ``spacer_seq``; a stage-1 tiles file calls it
#: ``guide_rna``. ``target_seq`` is deliberately **not** here -- it is the target,
#: not the guide, and showing it in a box labelled "guide" would be wrong.
SEQUENCE_COLUMNS = ("spacer_seq", "guide_rna")

#: Stage 2 columns carried into the hover box when --classification is given, in
#: the order they appear there. All three are the **mismatch-tolerant** (NM<=3)
#: view, so the label and the hits it was computed from are read from the same
#: cutoff and cannot be mistaken for each other.
#:
#: That whole view is **advisory**, not a complete cross-reactivity census: only
#: the NM=0 columns have an oracle behind them, so an absent NM<=3 hit means "not
#: found", never "does not exist" (see the stage 2 notes in CLAUDE.md). The
#: exact-match label lives in ``class_exact``, which is not shown -- pass it to
#: ``--label-column`` to get both.
CLASSIFICATION_COLUMNS = ("class_nm", "hit_genes_nm", "hit_locs_nm")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--fasta",
        required=True,
        metavar="PATH",
        help="the original mRNA transcript FASTA, one record. Its accession is "
        "matched against the transcript_id column of --predictions, and its "
        "length sets the extent of the x axis.",
    )
    parser.add_argument(
        "--predictions",
        required=True,
        metavar="PATH",
        help=f"predictions CSV, with at least {', '.join(REQUIRED_COLUMNS)} columns "
        "(e.g. results/ppib/ppib-predictions.csv)",
    )
    parser.add_argument(
        "--min-efficacy",
        type=float,
        required=True,
        metavar="FLOAT",
        help=f"keep guides with {EFFICACY_COLUMN} >= this value. Required -- the "
        "cut-off is the whole content of the plot, so it is never defaulted.",
    )
    parser.add_argument(
        "--output",
        metavar="PATH",
        help="HTML file to write (default: guides-{accession}-eff{threshold}.html "
        "beside --predictions, which keeps one run per threshold rather than "
        "overwriting the last)",
    )
    parser.add_argument(
        "--y-mode",
        choices=("lane", "score"),
        default="lane",
        help="'lane' (default) stacks overlapping guides onto separate rows, so "
        "extent and local density are readable; 'score' puts each guide at the "
        f"height of its {EFFICACY_COLUMN}, which shows quality and position at "
        "the cost of overlapping bars.",
    )
    parser.add_argument(
        "--label-column",
        metavar="NAME",
        help="extra column to show in the hover box, e.g. class_exact after a join "
        "against classification.tsv",
    )
    parser.add_argument(
        "--classification",
        metavar="PATH",
        help="optional stage 2 classification.tsv. When given, each guide's hover "
        f"box also shows {', '.join(CLASSIFICATION_COLUMNS)}, joined on tile_id. "
        "Omit it and the plot is built from --predictions alone.",
    )
    parser.add_argument(
        "--sequence-column",
        metavar="NAME",
        help="column holding the guide (spacer) sequence to show in the hover box. "
        f"Default: the first of {', '.join(SEQUENCE_COLUMNS)} present in the file, "
        "or no sequence if neither is. Pass 'none' to leave it out.",
    )
    parser.add_argument(
        "--inline-plotlyjs",
        action="store_true",
        help="embed plotly.js in the HTML (~3 MB) so it opens without a network "
        "connection; the default links it from a CDN (~50 kB file)",
    )
    return parser.parse_args(argv)


def read_transcript(path: str) -> tuple[str, int]:
    """The accession and length of the single record in ``path``.

    More than one record is a hard error rather than a silent "use the first":
    the whole plot is measured against this one length, so picking a record for
    the user would risk drawing guides against the wrong ruler.
    """
    try:
        # Opened here rather than handing SeqIO the path, so the handle is closed
        # even though parse() is lazy.
        with open(path, encoding="utf-8") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
    except FileNotFoundError:
        sys.exit(f"error: no such FASTA file: {path}")
    if not records:
        sys.exit(f"error: no FASTA records in {path}")
    if len(records) > 1:
        found = ", ".join(record.id for record in records[:5])
        sys.exit(
            f"error: {path} holds {len(records)} records ({found}...); --fasta takes "
            "the single transcript to plot. Inputs under data/ are one transcript "
            "per file."
        )
    record = records[0]
    return record.id, len(record.seq)


def read_predictions(path: str, accession: str) -> pd.DataFrame:
    """The prediction rows belonging to ``accession``.

    No matching rows is a hard error: it almost always means the FASTA and the
    predictions file are from different datasets, which would otherwise render as
    a plausible-looking empty plot.
    """
    try:
        frame = pd.read_csv(path)
    except FileNotFoundError:
        sys.exit(f"error: no such predictions file: {path}")
    except pd.errors.ParserError as exc:
        sys.exit(f"error: cannot parse {path} as CSV: {exc}")

    missing = [column for column in REQUIRED_COLUMNS if column not in frame.columns]
    if missing:
        sys.exit(
            f"error: {path} is missing required column(s): {', '.join(missing)}\n"
            f"       found: {', '.join(map(str, frame.columns))}"
        )

    for_transcript = frame[frame["transcript_id"].astype(str) == accession]
    if for_transcript.empty:
        present = sorted(frame["transcript_id"].astype(str).unique())
        sys.exit(
            f"error: no rows in {path} have transcript_id == {accession!r} (the "
            f"accession in the FASTA).\n       transcript_id values present: "
            f"{', '.join(present[:10])}{' ...' if len(present) > 10 else ''}\n"
            "       Are the FASTA and the predictions file from the same dataset?"
        )
    return for_transcript


def join_classification(
    frame: pd.DataFrame, path: str, *, extra_columns: Sequence[str] = ()
) -> tuple[pd.DataFrame, int]:
    """``frame`` with the stage 2 columns joined on ``tile_id``, and how many matched.

    A left join, so a guide missing from the classification simply carries fewer
    hover lines rather than vanishing from the plot -- the plot's job is to show
    the filtered guides, and it should not quietly drop any because a second file
    is incomplete.

    Zero matches is a hard error, for the same reason a mismatched predictions
    file is: tile_ids are the join key across every stage, so no overlap at all
    means the two files describe different tilings (a different dataset, or a
    different k) and every hover box would be silently unannotated.
    """
    try:
        classification = pd.read_csv(path, sep="\t")
    except FileNotFoundError:
        sys.exit(f"error: no such classification file: {path}")
    except pd.errors.ParserError as exc:
        sys.exit(f"error: cannot parse {path} as TSV: {exc}")

    missing = [
        column
        for column in ("tile_id", *CLASSIFICATION_COLUMNS)
        if column not in classification.columns
    ]
    if missing:
        sys.exit(
            f"error: {path} is missing required column(s): {', '.join(missing)}\n"
            f"       found: {', '.join(map(str, classification.columns))}\n"
            "       Is this a stage 2 classification.tsv?"
        )

    # ``extra_columns`` lets --label-column name a classification column such as
    # class_exact; one that is not in this file is ignored here and reported by
    # the caller against both files.
    columns = [
        "tile_id",
        *CLASSIFICATION_COLUMNS,
        *(
            column
            for column in extra_columns
            if column in classification.columns and column not in CLASSIFICATION_COLUMNS
        ),
    ]
    # Dropping duplicate tile_ids keeps the join from multiplying rows; stage 2
    # asserts one row per tile, so this is a guard, not an expected case.
    annotations = classification[columns].drop_duplicates(subset="tile_id")
    # A predictions file carrying one of these names would otherwise win the
    # merge: pandas suffixes the incoming column, leaving the stale value under
    # the name everything downstream reads. The classification file is
    # authoritative for its own columns, so the predictions copy goes first.
    shadowed = [column for column in columns if column != "tile_id" and column in frame.columns]
    joined = frame.drop(columns=shadowed).merge(annotations, on="tile_id", how="left")
    matched = int(joined[CLASSIFICATION_COLUMNS[0]].notna().sum())
    # An empty frame matches nothing for a legitimate reason -- the efficacy
    # threshold excluded every guide -- so only a non-empty one can be mismatched.
    if matched == 0 and not joined.empty:
        sys.exit(
            f"error: no tile_id in {path} matches the rows plotted from the "
            "predictions file.\n       tile_id is the join key across all stages, "
            "so no overlap means these files are different tilings (another "
            "dataset, or another k)."
        )
    return joined, matched


def resolve_sequence_column(requested: str | None, columns: Iterable[str]) -> str | None:
    """Which column holds the guide sequence, or ``None`` to show none.

    An explicit name must exist -- asking for a column and silently getting no
    sequence would be worse than an error. Without one, the first of
    :data:`SEQUENCE_COLUMNS` present is used, and a file with neither simply gets
    no sequence line.
    """
    present = list(columns)
    if requested is None:
        return next((name for name in SEQUENCE_COLUMNS if name in present), None)
    if requested.lower() == "none":
        return None
    if requested not in present:
        sys.exit(
            f"error: --sequence-column {requested!r} is not a column in the "
            f"predictions file; found: {', '.join(map(str, present))}"
        )
    return requested


def default_output(predictions: str, accession: str, threshold: float) -> str:
    """``guides-{accession}-eff{threshold}.html`` beside the predictions file.

    The threshold is in the name on purpose: the plot is meaningless without it,
    and two cut-offs should not overwrite one another.
    """
    stem = f"guides-{accession}-eff{threshold:g}.html"
    return os.path.join(os.path.dirname(os.path.abspath(predictions)), stem)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    accession, length = read_transcript(args.fasta)
    for_transcript = read_predictions(args.predictions, accession)

    sequence_column = resolve_sequence_column(args.sequence_column, for_transcript.columns)

    efficacy = pd.to_numeric(for_transcript[EFFICACY_COLUMN], errors="coerce")
    unparsed = int(efficacy.isna().sum())
    kept = for_transcript[efficacy >= args.min_efficacy]

    # Joined after filtering: only the guides being plotted need annotating, and
    # `matched` then counts those, which is the number worth reporting.
    annotation_columns: tuple[str, ...] = ()
    matched = 0
    if args.classification:
        kept, matched = join_classification(
            kept,
            args.classification,
            extra_columns=(args.label_column,) if args.label_column else (),
        )
        annotation_columns = CLASSIFICATION_COLUMNS

    # Checked after the join, so --label-column can name a classification column
    # (class_exact, say) and not only one of the predictions file's own.
    if args.label_column and args.label_column not in kept.columns:
        sources = args.predictions
        if args.classification:
            sources += f" or {args.classification}"
        sys.exit(f"error: --label-column {args.label_column!r} is not a column in {sources}")

    guides = guides_from_rows(
        kept.to_dict("records"),
        transcript_id=accession,
        score_column=EFFICACY_COLUMN,
        label_column=args.label_column,
        sequence_column=sequence_column,
        annotation_columns=annotation_columns,
    )

    try:
        figure = guide_track(
            length,
            guides,
            transcript_id=accession,
            y_mode=args.y_mode,
            score_label=EFFICACY_COLUMN.replace("_", " "),
        )
    except ValueError as exc:
        sys.exit(f"error: {exc}")

    output = args.output or default_output(args.predictions, accession, args.min_efficacy)
    directory = os.path.dirname(os.path.abspath(output))
    os.makedirs(directory, exist_ok=True)
    write_html(figure, output, include_plotlyjs=True if args.inline_plotlyjs else "cdn")

    report(
        args, accession, length, for_transcript, efficacy, guides, unparsed, output,
        sequence_column, matched,
    )
    return 0


def report(
    args: argparse.Namespace,
    accession: str,
    length: int,
    for_transcript: pd.DataFrame,
    efficacy: pd.Series,
    guides: list,
    unparsed: int,
    output: str,
    sequence_column: str | None,
    matched: int,
) -> None:
    # read_predictions exits when no rows match the accession, so total >= 1 here.
    total = len(for_transcript)
    kept = len(guides)
    observed = efficacy.dropna()
    print(f"transcript: {accession}  length: {length:,} nt  ({args.fasta})", file=sys.stderr)
    print(
        f"predictions: {total:,} row(s) for this transcript in {args.predictions}",
        file=sys.stderr,
    )
    if unparsed:
        print(
            f"warning: {unparsed:,} row(s) have a non-numeric {EFFICACY_COLUMN} and "
            "were treated as failing the threshold",
            file=sys.stderr,
        )
    print(
        f"filter: {EFFICACY_COLUMN} >= {args.min_efficacy:g}  ->  kept {kept:,} "
        f"({kept / total:.1%} of rows)",
        file=sys.stderr,
    )
    if kept:
        covered = len({position for guide in guides for position in range(guide.start, guide.end + 1)})
        print(
            f"guides span {min(g.start for g in guides):,}-{max(g.end for g in guides):,}; "
            f"{covered:,} of {length:,} nt covered ({covered / length:.1%}) in "
            f"{max(pack_lanes(guides)) + 1} lane(s)",
            file=sys.stderr,
        )
        print(
            f"{EFFICACY_COLUMN} of kept guides: {min(g.score for g in guides):.4g} "
            f"to {max(g.score for g in guides):.4g}  "
            f"(all rows: {observed.min():.4g} to {observed.max():.4g})",
            file=sys.stderr,
        )
    else:
        highest = f"{observed.max():.4g}" if not observed.empty else "n/a"
        print(
            f"note: nothing passed the threshold -- the highest {EFFICACY_COLUMN} on "
            f"this transcript is {highest}. The plot shows the transcript with no "
            "guides on it.",
            file=sys.stderr,
        )
    if args.classification:
        print(
            f"classification: {args.classification} -> {matched:,} of {kept:,} plotted "
            f"guide(s) annotated with {', '.join(CLASSIFICATION_COLUMNS)}",
            file=sys.stderr,
        )
        if matched < kept:
            print(
                f"warning: {kept - matched:,} plotted guide(s) have no row in the "
                "classification file and so carry no classification in the hover "
                "box; was it generated from the same tiling?",
                file=sys.stderr,
            )
        print(
            f"note: all three ({', '.join(CLASSIFICATION_COLUMNS)}) are the NM<=3 "
            "view, which rests on bwa's heuristic with no oracle behind it and is "
            "advisory -- an absent hit means 'not found', not 'does not exist'. The "
            "validated exact-match label is class_exact; --label-column shows it too.",
            file=sys.stderr,
        )
    print(
        f"y axis: {args.y_mode}  hover sequence: "
        + (f"{sequence_column} (shown 5' -> 3')" if sequence_column else "none"),
        file=sys.stderr,
    )
    if sequence_column is None and args.sequence_column is None:
        print(
            "note: no guide sequence in the hover box -- the predictions file has "
            f"none of {', '.join(SEQUENCE_COLUMNS)}. Name the column with "
            "--sequence-column.",
            file=sys.stderr,
        )
    print(f"wrote {output}", file=sys.stderr)
    if not args.inline_plotlyjs:
        print(
            "plotly.js is loaded from a CDN, so opening the file needs a network "
            "connection; pass --inline-plotlyjs for a self-contained copy.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    raise SystemExit(main())
