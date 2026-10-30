#!/usr/bin/env python3
"""Plot only the selected guides along each transcript, each over its full length.

    python scripts/plot_selected_guides.py \
        --selected results/selected_guides.list \
        --target-sites results/pik3/target_site_selection.tsv \
        --predictions results/pik3/pik3-predictions.csv results/ppib/ppib-predictions.csv \
        --blastn results/blastn/siRNA_targetSeq.blastn \
        --conserved-dir data/pik3-iso \
        --genbank-dir data/pik3-iso --genbank-dir data/ppib \
        --output-dir results/plots-selected

One HTML plot per transcript that has selected guides. Exon and CDS lanes, siRNA
bands and hover content are those of ``plot_scatterlines.py``; the difference is
that only the guides named in ``--selected`` (by ``tile_id``) are drawn, each as
a horizontal segment spanning its whole length at the height of its predicted
efficacy. The target region from ``--target-sites`` and the conserved regions from the
``PIK3*_conserved_sites.bed`` files each get a lane above the plot.
"""

from __future__ import annotations

import argparse
import glob
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
    BLASTN_COLUMNS,
    DEFAULT_THRESHOLD,
    PREDICTION_COLUMNS,
    build_guides,
    build_sirnas,
    transcript_lengths,
)
from visualise_selected import plot_selected_guides, read_bed_regions, read_guide_list  # noqa: E402

TARGET_SITE_COLUMNS = ("transcript", "start", "end")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--selected", default=os.path.join("results", "selected_guides.list"),
                        metavar="PATH", help="one tile_id per line (default: %(default)s)")
    parser.add_argument("--target-sites",
                        default=os.path.join("results", "pik3", "target_site_selection.tsv"),
                        metavar="PATH",
                        help="TSV with transcript, start, end columns, shaded as the target "
                        "region (default: %(default)s)")
    parser.add_argument("--predictions", nargs="+", metavar="PATH",
                        default=[os.path.join("results", "pik3", "pik3-predictions.csv"),
                                 os.path.join("results", "ppib", "ppib-predictions.csv")],
                        help="one or more predictions CSVs; the selected guides are looked up "
                        "across all of them (default: %(default)s)")
    parser.add_argument("--blastn", default=os.path.join("results", "blastn", "siRNA_targetSeq.blastn"),
                        metavar="PATH", help="siRNA blastn table (default: %(default)s)")
    parser.add_argument("--genbank-dir", action="append", default=None, metavar="DIR",
                        help="directory of {accession}.gbk records; repeatable "
                        "(default: data/pik3-iso and data/ppib)")
    parser.add_argument("--conserved-dir", default=os.path.join("data", "pik3-iso"), metavar="DIR",
                        help="directory of PIK3*_conserved_sites.bed files (BED, 0-based half-open); "
                        "transcripts with no record there, such as PPIB, get no conserved lane "
                        "(default: %(default)s)")
    parser.add_argument("--table", default=os.path.join("results", "selected_guides.tsv"),
                        metavar="PATH", help="TSV of the selected guides: tile_id, transcript, "
                        "gene, guide_seq, target_seq, predicted_efficacy (default: %(default)s)")
    parser.add_argument("--output-dir", default=os.path.join("results", "plots-selected"),
                        metavar="DIR", help="default: %(default)s")
    parser.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD, metavar="X",
                        help="efficacy cutoff for the fill colours (default: %(default)s)")
    parser.add_argument("--cdn-plotlyjs", action="store_true",
                        help="link plotly.js from a CDN instead of embedding it")
    return parser.parse_args(argv)


def require_columns(df: pd.DataFrame, required, source: str) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise SystemExit(
            f"error: {source} is missing required column(s) {', '.join(missing)}.\n"
            f"       columns found: {', '.join(map(str, df.columns))}"
        )


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    genbank_dirs = args.genbank_dir or [os.path.join("data", "pik3-iso"), os.path.join("data", "ppib")]

    for flag, path in (("--selected", args.selected), ("--target-sites", args.target_sites),
                       ("--blastn", args.blastn), *(("--predictions", p) for p in args.predictions)):
        if not os.path.isfile(path):
            raise SystemExit(f"error: {flag} {path}: no such file")
    genbank_dirs = [d for d in genbank_dirs if os.path.isdir(d)]

    selected = read_guide_list(args.selected)
    if len(set(selected)) != len(selected):
        print("warning: duplicate ids in the selected list are plotted once", file=sys.stderr)
    selected_set = set(selected)

    frames = []
    for path in args.predictions:
        df = pd.read_csv(path, header=0)
        require_columns(df, PREDICTION_COLUMNS, path)
        frames.append(df)
    predictions = pd.concat(frames, ignore_index=True)
    # Lengths come from the full tiling, before selection, so the axis spans the transcript.
    lengths = transcript_lengths(predictions)

    chosen = predictions[predictions["tile_id"].isin(selected_set)].drop_duplicates("tile_id")
    absent = sorted(selected_set - set(chosen["tile_id"]))
    if absent:
        raise SystemExit(
            f"error: {len(absent)} selected guide(s) not found in the predictions: "
            f"{', '.join(absent)}"
        )

    sites = pd.read_csv(args.target_sites, sep=r"\s+", engine="python")
    require_columns(sites, TARGET_SITE_COLUMNS, args.target_sites)
    genes = dict(zip(sites["transcript"], sites["gene"])) if "gene" in sites.columns else {}
    regions = {r.transcript: (int(r.start), int(r.end)) for r in sites.itertuples()}

    # Written in the order of the selected list, so it reads like the list it came from.
    order = {tile: i for i, tile in enumerate(dict.fromkeys(selected))}
    table = chosen.assign(
        transcript=chosen["transcript_id"],
        gene=chosen["transcript_id"].map(genes).fillna("-"),
        guide_seq=chosen["spacer_seq"],
    ).sort_values("tile_id", key=lambda col: col.map(order))
    table[["tile_id", "transcript", "gene", "guide_seq", "target_seq", "predicted_efficacy"]].to_csv(
        args.table, sep="\t", index=False
    )
    print(f"wrote {len(table)} guides -> {args.table}", file=sys.stderr)

    conserved: dict[str, list[tuple[int, int]]] = {}
    for bed in sorted(glob.glob(os.path.join(args.conserved_dir, "PIK3*_conserved_sites.bed"))):
        for reference, spans in read_bed_regions(bed).items():
            conserved.setdefault(reference, []).extend(spans)

    blastn_df = pd.read_csv(args.blastn, sep=r"\s+", header=0, engine="python")
    require_columns(blastn_df, BLASTN_COLUMNS, args.blastn)
    sirnas = build_sirnas(blastn_df)

    guides = [g for g in build_guides(chosen, sirnas)]
    # build_guides emits one copy per siRNA; the plot picks the nearest itself.
    unique = list({g.tile_id: g._replace(sirna=None, classification=None, distance=None)
                   for g in guides}.values())

    os.makedirs(args.output_dir, exist_ok=True)
    print(f"{len(selected_set)} selected guides over "
          f"{chosen['transcript_id'].nunique()} transcript(s)", file=sys.stderr)

    for transcript in sorted(chosen["transcript_id"].unique()):
        exons, cds, symbol, record_length = (), None, None, None
        genbank = find_genbank(transcript, genbank_dirs) if genbank_dirs else None
        if genbank:
            exons = exons_from_genbank(genbank)
            cds = cds_from_genbank(genbank)
            symbol = gene_symbol_from_genbank(genbank)
            record_length = genbank_length(genbank)
        else:
            print(f"  warning: no {transcript}.gbk found; no annotation lanes", file=sys.stderr)
        if record_length is not None and record_length != lengths[transcript]:
            print(f"  warning: {transcript} GenBank record is {record_length:,} nt but the "
                  f"predictions imply {lengths[transcript]:,} nt; the lanes are misplaced",
                  file=sys.stderr)

        transcript_guides = [g for g in unique if g.transcript == transcript]
        transcript_sirnas = [s for s in sirnas if s.transcript == transcript]
        output = os.path.join(
            args.output_dir,
            f"{transcript}_{symbol}_selected.html" if symbol else f"{transcript}_selected.html",
        )
        plot_selected_guides(
            transcript_guides, transcript_sirnas, lengths[transcript], output,
            inline_plotlyjs=not args.cdn_plotlyjs, exons=exons, cds=cds,
            threshold=args.threshold, target_region=regions.get(transcript),
            conserved=conserved.get(transcript, ()), gene=symbol,
        )
        print(f"  {transcript}: {len(transcript_guides)} guides, "
              f"{len(transcript_sirnas)} siRNA site(s) -> {output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
