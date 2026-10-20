#!/usr/bin/env python3
"""Plot guide efficacy along each transcript, against known siRNA target sites.

    python scripts/plot_scatterlines.py \
        --blastn results/blastn/siRNA_targetSeq.blastn \
        --predictions results/pik3/pik3-predictions.csv \
        --output-dir results/plots

Stages:

1. Read the blastn table of siRNAs mapped against the transcripts and keep one
   site per siRNA -- its highest-scoring HSP. A minus-strand hit is reported
   with ``sstart > send``, so the span is ordered before anything uses it.
2. Read the predictions CSV and classify every guide as ``upstream``,
   ``downstream`` or ``overlap`` **relative to one siRNA**. A transcript can
   carry several siRNAs, so a guide is emitted once per siRNA and each copy
   carries its own classification; guides on a transcript with no siRNA are
   kept unclassified rather than dropped.
3. Write one HTML plot per transcript: every guide as a point at its start
   coordinate and at the height of its predicted efficacy, filled green-to-
   yellow by that efficacy, with the guides overlapping an siRNA ringed in
   charcoal and the siRNA footprints shaded light red.

**One plot is one transcript.** The transcripts come from the predictions file,
so a file covering four paralogs gives four plots. Every one is drawn across
the full transcript, 1 to its length, so an uncovered stretch reads as
genuinely uncovered.

The transcript length is taken as the largest tile ``end`` in the predictions
file. Stage 1 tiles at step 1 across the whole transcript, so that is exact for
a full tiling -- but it is a **floor**, not a measurement: pass an
already-filtered predictions file and the axis stops at the last surviving
guide rather than at the end of the transcript.

Nothing here filters by efficacy. Every row in the predictions file is plotted,
and the colour scale is pinned to 0-1 so two plots can be compared.
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

import pandas as pd  # noqa: E402

from visualise_scatterlines import (  # noqa: E402
    BLASTN_COLUMNS,
    CLASS_STYLE,
    PREDICTION_COLUMNS,
    build_guides,
    build_sirnas,
    plot_scatterlines,
    transcript_lengths,
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--blastn",
        required=True,
        metavar="PATH",
        help="blastn table (outfmt 6 with a header) of the siRNAs mapped against "
        f"the transcripts, with at least {', '.join(BLASTN_COLUMNS)} columns. The "
        "best HSP per qseqid is kept, so an unfiltered table is fine.",
    )
    parser.add_argument(
        "--predictions",
        required=True,
        metavar="PATH",
        help="predictions CSV, with at least "
        f"{', '.join(PREDICTION_COLUMNS)} columns "
        "(e.g. results/pik3/pik3-predictions.csv)",
    )
    parser.add_argument(
        "--output-dir",
        default=os.path.join("results", "plots"),
        metavar="DIR",
        help="directory for the HTML plots, one per transcript, named "
        "{transcript}-scatterlines.html (default: %(default)s). Created if "
        "absent.",
    )
    parser.add_argument(
        "--cdn-plotlyjs",
        action="store_true",
        help="link plotly.js from a CDN (~50 kB per file) instead of embedding "
        "it (~3 MB); the small file needs a network connection to open",
    )
    return parser.parse_args(argv)


def require_columns(df: pd.DataFrame, required, source: str) -> None:
    """Fail up front, naming the columns actually found.

    Checked here rather than in the library because this is where a path the
    user typed turns into a DataFrame -- and a wrong file is far likelier than
    a wrong call.
    """
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise SystemExit(
            f"error: {source} is missing required column(s) "
            f"{', '.join(missing)}.\n"
            f"       columns found: {', '.join(map(str, df.columns))}"
        )


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    # Checked before pandas opens them: a mistyped path is the likeliest thing
    # to go wrong here, and pandas answers it with a bare traceback.
    for flag, path in (("--blastn", args.blastn),
                       ("--predictions", args.predictions)):
        if not os.path.isfile(path):
            raise SystemExit(f"error: {flag} {path}: no such file")

    # Split on any run of whitespace rather than on tabs alone: no blastn
    # outfmt 6 field contains a space, and these tables are often hand-filtered
    # into a mix of tabs and spaces that a strict sep="\t" folds into one
    # column, giving a baffling "missing column" error on a file that looks fine.
    blastn_df = pd.read_csv(args.blastn, sep=r"\s+", header=0, engine="python")
    require_columns(blastn_df, BLASTN_COLUMNS, args.blastn)
    sirnas = build_sirnas(blastn_df)

    prediction_df = pd.read_csv(args.predictions, header=0)
    require_columns(prediction_df, PREDICTION_COLUMNS, args.predictions)

    lengths = transcript_lengths(prediction_df)
    guides = build_guides(prediction_df, sirnas)

    os.makedirs(args.output_dir, exist_ok=True)

    print(
        f"{len(prediction_df)} guides over {len(lengths)} transcript(s), "
        f"{len(sirnas)} siRNA site(s) -> {len(guides)} guide-siRNA pairings",
        file=sys.stderr,
    )

    for transcript in sorted(lengths):
        transcript_sirnas = [s for s in sirnas if s.transcript == transcript]
        transcript_guides = [g for g in guides if g.transcript == transcript]
        output = os.path.join(
            args.output_dir, f"{transcript}-scatterlines.html"
        )
        plot_scatterlines(
            transcript_guides,
            transcript_sirnas,
            {transcript: lengths[transcript]},
            output,
            inline_plotlyjs=not args.cdn_plotlyjs,
        )
        overlaps = sum(1 for g in transcript_guides if g.classification == "overlap")
        print(
            f"  {transcript}: {lengths[transcript]:,} nt, "
            f"{len(transcript_guides):,} points, "
            f"{len(transcript_sirnas)} siRNA site(s), {overlaps} overlapping "
            f"-> {output}",
            file=sys.stderr,
        )
        if not transcript_sirnas:
            print(
                "      no siRNA mapped to this transcript, so no band is drawn "
                "and every guide is unclassified",
                file=sys.stderr,
            )

    # An siRNA on a transcript absent from the predictions can never be
    # plotted, which usually means the two files describe different datasets.
    orphans = sorted(s.qseqid for s in sirnas if s.transcript not in lengths)
    if orphans:
        print(
            f"warning: {len(orphans)} siRNA(s) map to transcripts absent from "
            f"{args.predictions} and were not plotted: {', '.join(orphans)}",
            file=sys.stderr,
        )
    if not any(s.transcript in lengths for s in sirnas):
        print(
            "warning: not one siRNA maps to a transcript in the predictions "
            "file. No plot has an siRNA band, so the two files may describe "
            "different datasets.",
            file=sys.stderr,
        )

    # A class with no points is left out of the plot, so say what was drawn.
    # Accumulated rather than built by comprehension: upstream, downstream and
    # the no-siRNA class share one label, so keying on it has to sum them, not
    # overwrite.
    counts: dict[str, int] = {}
    for cls, (label, _, _) in CLASS_STYLE.items():
        n = sum(1 for g in guides if g.classification == cls)
        counts[label] = counts.get(label, 0) + n
    print(
        "classes drawn: "
        + ", ".join(f"{k}={v:,}" for k, v in counts.items() if v),
        file=sys.stderr,
    )
    if args.cdn_plotlyjs:
        print(
            "plotly.js is loaded from a CDN, so opening the files needs a "
            "network connection; omit --cdn-plotlyjs for self-contained copies.",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
