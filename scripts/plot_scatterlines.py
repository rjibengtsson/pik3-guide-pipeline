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

from annotation import (  # noqa: E402
    cds_from_genbank,
    exons_from_genbank,
    find_genbank,
    genbank_length,
    gene_symbol_from_genbank,
)
from visualise_scatterlines import (  # noqa: E402
    DEFAULT_THRESHOLD,
    BLASTN_COLUMNS,
    CLASS_STYLE,
    PREDICTION_COLUMNS,
    build_guides,
    build_sirnas,
    plot_scatterlines,
    highlight_category,
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
        "{accession}_{gene}_anno.html when a GenBank record supplies the "
        "gene symbol -- e.g. NM_006219.3_PIK3CB_anno.html -- and "
        "{accession}-scatterlines.html when none is available, since that "
        "plot carries no annotation (default: %(default)s). Created if "
        "absent.",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=DEFAULT_THRESHOLD,
        metavar="X",
        help="predicted-efficacy cutoff for highlighting (default: "
        "%(default)s). Guides at or below it are grey; those strictly above "
        "it are yellow, and the ones among those whose spacer starts with GG "
        "are orange. The scale is top-heavy, so a low threshold highlights "
        "nearly everything: at 0.7 about 79%% of guides clear it.",
    )
    parser.add_argument(
        "--genbank-dir",
        action="append",
        default=None,
        metavar="DIR",
        help="directory of RefSeq mRNA GenBank records named {accession}.gbk. "
        "Repeatable, and searched in the order given -- the records are "
        "organised per dataset, so the PI3K paralogs and PPIB live in "
        "different directories: --genbank-dir data/pik3-iso --genbank-dir "
        "data/ppib. Each transcript with a record gets an exon lane and a CDS "
        "lane above its scatter; a transcript without one is plotted bare and "
        "named in the summary.",
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

    genbank_dirs = args.genbank_dir or []
    for directory in genbank_dirs:
        if not os.path.isdir(directory):
            raise SystemExit(f"error: --genbank-dir {directory}: no such directory")

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

    missing_genbank: list[str] = []

    for transcript in sorted(lengths):
        transcript_sirnas = [s for s in sirnas if s.transcript == transcript]
        transcript_guides = [g for g in guides if g.transcript == transcript]

        # A missing record is a warning, not a failure: the scatter is still
        # correct without its lanes, and losing every figure because one
        # annotation file is absent is the wrong trade. A *malformed* record
        # does raise -- see the contiguity check in src/annotation.py.
        exons, cds, record_length, symbol = (), None, None, None
        if genbank_dirs:
            genbank = find_genbank(transcript, genbank_dirs)
            if genbank:
                exons = exons_from_genbank(genbank)
                cds = cds_from_genbank(genbank)
                record_length = genbank_length(genbank)
                symbol = gene_symbol_from_genbank(genbank)
            else:
                missing_genbank.append(transcript)

        # {accession}_{gene}_anno.html once there is annotation to put in the
        # plot, since that is what "anno" claims. Without a GenBank record
        # there is no gene symbol and no lane, so the name falls back to the
        # original form rather than inventing a symbol or writing a file whose
        # name promises annotation it does not contain.
        output = os.path.join(
            args.output_dir,
            f"{transcript}_{symbol}_anno.html" if symbol
            else f"{transcript}-scatterlines.html",
        )
        plot_scatterlines(
            transcript_guides,
            transcript_sirnas,
            {transcript: lengths[transcript]},
            output,
            inline_plotlyjs=not args.cdn_plotlyjs,
            exons=exons,
            cds=cds,
            threshold=args.threshold,
        )
        overlaps = sum(1 for g in transcript_guides if g.classification == "overlap")
        # Counted over distinct guides, not over the guide-siRNA pairs the plot
        # draws: a guide on a transcript with three siRNAs is plotted three
        # times, so counting pairs would treble these and they would not match
        # the predictions file.
        unique = {g.tile_id: g for g in transcript_guides}.values()
        states = [highlight_category(g, args.threshold) for g in unique]
        print(
            f"  {transcript}: {lengths[transcript]:,} nt, "
            f"{len(transcript_guides):,} points, "
            f"{len(transcript_sirnas)} siRNA site(s), {overlaps} overlapping"
            + (f", {len(exons)} exons" if exons else "")
            + (f", CDS {cds.start:,}-{cds.end:,}" if cds else "")
            + f" -> {output}",
            file=sys.stderr,
        )
        print(
            f"      {states.count('above') + states.count('gg'):,} of "
            f"{len(states):,} guides above {args.threshold:g}, of which "
            f"{states.count('gg'):,} start with GG",
            file=sys.stderr,
        )
        # Two different disagreements, with opposite consequences, so they get
        # separate messages rather than one guess at the cause.
        #
        # Record against predictions: the .gbk and the predictions are built
        # from different versions of the accession, so every feature in both
        # lanes is displaced by however much the transcript changed.
        if record_length is not None and record_length != lengths[transcript]:
            print(
                f"      warning: the GenBank record is {record_length:,} nt "
                f"but the predictions imply {lengths[transcript]:,} nt. These "
                "are different versions of the accession, so the lanes are "
                "misplaced.",
                file=sys.stderr,
            )
        # Exons against their own record: NCBI's exon table stops short of the
        # sequence it annotates. The exons that exist are placed correctly;
        # the uncovered tail simply has none.
        elif exons and record_length is not None and exons[-1].end != record_length:
            print(
                f"      note: the exon table covers 1..{exons[-1].end:,} of "
                f"{record_length:,} nt, so the last "
                f"{record_length - exons[-1].end} nt have no exon drawn. This "
                "is a gap in the record's annotation, not a coordinate error.",
                file=sys.stderr,
            )
        if not transcript_sirnas:
            print(
                "      no siRNA mapped to this transcript, so no band is drawn "
                "and every guide is unclassified",
                file=sys.stderr,
            )

    if missing_genbank:
        print(
            "warning: no {accession}.gbk found in "
            f"{', '.join(genbank_dirs)} for {', '.join(missing_genbank)}; "
            "those plots have no annotation lanes",
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
