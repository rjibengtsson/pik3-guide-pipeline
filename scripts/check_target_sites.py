#!/usr/bin/env python3
"""Check the selected target regions for GC, poly-T, repeats and stop codons.

    python scripts/check_target_sites.py \
        --selection results/pik3/target_site_selection.tsv \
        --input-dir data/pik3-iso --input-dir data/ppib \
        --output-dir results/pik3

``--input-dir`` may be repeated and is searched in order, which is what a
selection table spanning datasets needs: the PI3K FASTAs live in
``data/pik3-iso`` and PPIB in ``data/ppib``, and no single directory holds both.

Reads the selection table, pulls each region out of its transcript FASTA, and
writes three files:

``target_site_qc.tsv``
    One row per region: every summary number, plus the thresholds that judged
    it, so a flag can be audited from its own row without rerunning anything.
``target_site_qc_features.tsv``
    One row per located feature. The per-region row collapses things; this is
    the file to debug from, and it carries each feature's coordinates **twice**
    -- as 1-based inclusive transcript coordinates, which read directly against
    ``tile_id`` and the FASTA, and as an offset within the region, so a feature
    can be found in the extracted sequence without arithmetic.
``target_site_qc.md``
    The same thing as prose, for the notebook.

Coordinates are 1-based inclusive throughout, following src/tile.py, so
``sequence[start - 1:end]`` is the region and ``end - start + 1`` its length.

The declared ``length`` column is checked against ``end - start + 1`` and the
mismatch is reported in ``length_actual``/``length_matches`` and as a warning.
The run then continues over the **actual** span: a selection table disagreeing
with itself is worth stopping a pipeline for, but not worth refusing to
describe, and the actual span is what the coordinates denote.

Stop codons are reported as counts and positions only -- no frame, no CDS, and
no flag. Chance alone puts ~94 stop triplets in 2 kb, so "has stop codons"
would be true of every region, correct ones included. The flags come from GC,
poly-T and the repeat/entropy checks.
"""

from __future__ import annotations

import argparse
import csv
import os
import statistics
import sys
from typing import Sequence

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

import seqqc  # noqa: E402
import tile as tile_mod  # noqa: E402

SELECTION_COLUMNS = ("transcript", "gene", "start", "end", "length")

SUMMARY_COLUMNS = (
    "transcript", "gene", "start", "end",
    "length_declared", "length_actual", "length_matches",
    "n_ambiguous",
    "gc_frac", "gc_min_window", "gc_max_window", "gc_mean_window", "gc_sd_window",
    "n_windows_out_of_band", "widest_out_of_band_run", "gc_band",
    "n_polyT_runs", "longest_polyT", "longest_polyT_at", "polyT_min",
    "n_homopolymers", "longest_homopolymer", "longest_homopolymer_base",
    "homopolymer_min",
    "n_tandem_repeats", "longest_tandem_repeat", "longest_tandem_motif",
    "min_window_entropy", "min_entropy_at", "n_windows_low_entropy", "entropy_min",
    "n_stop_TAA", "n_stop_TAG", "n_stop_TGA", "n_stops_total", "stops_per_100nt",
    "flags", "n_flags",
)

FEATURE_COLUMNS = (
    "transcript", "gene", "feature_type",
    "start", "end", "offset_in_region", "length", "detail", "severity",
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--selection",
        default=os.path.join("results", "pik3", "target_site_selection.tsv"),
        metavar="PATH",
        help="selection table, tab-separated, with columns "
        f"{', '.join(SELECTION_COLUMNS)} (default: %(default)s)",
    )
    parser.add_argument(
        "--input-dir",
        action="append",
        metavar="DIR",
        help="directory of transcript FASTAs named {accession}.fasta "
        f"(default: {os.path.join('data', 'pik3-iso')}). May be repeated, and "
        "is then searched in the order given -- a selection table spanning two "
        "datasets needs one of these per dataset directory.",
    )
    parser.add_argument(
        "--output-dir",
        default=os.path.join("results", "pik3"),
        metavar="DIR",
        help="directory for the three output files (default: %(default)s). "
        "Created if absent.",
    )
    parser.add_argument(
        "--window", type=int, default=seqqc.DEFAULT_WINDOW, metavar="N",
        help="width of the GC and entropy windows, in bases (default: "
        "%(default)s, a spacer length, so one window is one candidate guide)",
    )
    parser.add_argument(
        "--gc-min", type=float, default=seqqc.DEFAULT_GC_MIN, metavar="X",
        help="lower edge of the acceptable GC band, as a fraction "
        "(default: %(default)s)",
    )
    parser.add_argument(
        "--gc-max", type=float, default=seqqc.DEFAULT_GC_MAX, metavar="X",
        help="upper edge of the acceptable GC band (default: %(default)s)",
    )
    parser.add_argument(
        "--polyt-min", type=int, default=seqqc.DEFAULT_POLYT_MIN, metavar="N",
        help="shortest run of T to report (default: %(default)s, the "
        "premature-termination signal in a U6-driven construct)",
    )
    parser.add_argument(
        "--homopolymer-min", type=int, default=seqqc.DEFAULT_HOMOPOLYMER_MIN,
        metavar="N",
        help="shortest run of A, C or G to report (default: %(default)s). T is "
        "covered by --polyt-min and is not reported twice.",
    )
    parser.add_argument(
        "--repeat-copies", type=int, default=seqqc.DEFAULT_REPEAT_COPIES,
        metavar="N",
        help="fewest consecutive copies of a 2-4 nt motif to call a tandem "
        "repeat (default: %(default)s)",
    )
    parser.add_argument(
        "--min-entropy", type=float, default=seqqc.DEFAULT_MIN_ENTROPY,
        metavar="X",
        help="flag windows whose Shannon entropy falls below this, in bits out "
        "of a maximum of 2 (default: %(default)s)",
    )
    args = parser.parse_args(argv)
    # append leaves it None rather than the default when the flag is absent.
    if not args.input_dir:
        args.input_dir = [os.path.join("data", "pik3-iso")]
    return args


def read_selection(path: str) -> list[dict]:
    """Parse the selection table, failing on a header that does not line up.

    The header has been malformed before -- ``start`` and ``end`` joined by
    spaces instead of a tab, giving four header fields over five-field rows,
    which ``csv.DictReader`` turns into a column literally named
    ``"start    end"`` and a silently misread table. So the field count is
    checked per row rather than trusted.
    """
    with open(path, newline="", encoding="utf-8") as handle:
        rows = [line.rstrip("\n").split("\t") for line in handle if line.strip()]

    if not rows:
        raise SystemExit(f"error: {path}: empty")

    header = rows[0]
    missing = [c for c in SELECTION_COLUMNS if c not in header]
    if missing:
        raise SystemExit(
            f"error: {path}: header is missing column(s) {', '.join(missing)}.\n"
            f"       header fields found ({len(header)}): {header}\n"
            "       The columns must be separated by tabs, not spaces."
        )

    index = {name: header.index(name) for name in SELECTION_COLUMNS}
    records = []
    for number, row in enumerate(rows[1:], start=2):
        if len(row) != len(header):
            raise SystemExit(
                f"error: {path} line {number}: {len(row)} fields but the header "
                f"has {len(header)}. Mixed tabs and spaces are the usual cause."
            )
        records.append({name: row[position] for name, position in index.items()})
    return records


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    if not os.path.isfile(args.selection):
        raise SystemExit(f"error: --selection {args.selection}: no such file")
    for directory in args.input_dir:
        if not os.path.isdir(directory):
            raise SystemExit(f"error: --input-dir {directory}: no such directory")

    selection = read_selection(args.selection)
    os.makedirs(args.output_dir, exist_ok=True)

    summaries: list[dict] = []
    features: list[dict] = []
    mismatched: list[str] = []

    for record in selection:
        accession = record["transcript"]
        gene = record["gene"]
        start, end = int(record["start"]), int(record["end"])
        declared = int(record["length"])

        candidates = [
            os.path.join(directory, f"{accession}.fasta")
            for directory in args.input_dir
        ]
        fasta = next((path for path in candidates if os.path.isfile(path)), None)
        if fasta is None:
            raise SystemExit(
                f"error: no FASTA for {gene} / {accession}; looked for "
                + ", ".join(candidates)
            )
        records = list(tile_mod.read_fasta(fasta))
        if len(records) != 1:
            raise SystemExit(
                f"error: {fasta}: expected one record, found {len(records)}"
            )
        sequence = records[0].sequence

        if end > len(sequence):
            raise SystemExit(
                f"error: {gene}: region ends at {end:,} but {accession} is only "
                f"{len(sequence):,} nt"
            )

        # The whole convention in one line: 1-based inclusive, so this slice is
        # exactly the region the table names.
        region = sequence[start - 1:end]
        actual = end - start + 1
        if actual != declared:
            mismatched.append(
                f"{gene} {start:,}-{end:,}: declared {declared:,} nt, "
                f"actual {actual:,} nt"
            )

        summary, region_features = analyse(region, accession, gene, start, args)
        summary["length_declared"] = declared
        summary["length_actual"] = actual
        summary["length_matches"] = int(actual == declared)
        summaries.append(summary)
        features.extend(region_features)

    summary_path = os.path.join(args.output_dir, "target_site_qc.tsv")
    features_path = os.path.join(args.output_dir, "target_site_qc_features.tsv")
    report_path = os.path.join(args.output_dir, "target_site_qc.md")

    write_tsv(summary_path, SUMMARY_COLUMNS, summaries)
    write_tsv(features_path, FEATURE_COLUMNS, features)
    write_report(report_path, summaries, features, args)

    print(
        f"{len(summaries)} region(s) checked from {args.selection}",
        file=sys.stderr,
    )
    for summary in summaries:
        print(
            f"  {summary['gene']:7} {summary['start']:>6,}-{summary['end']:<7,} "
            f"{summary['length_actual']:>5,} nt  "
            f"GC {summary['gc_frac']:.3f}  "
            f"windows out of band {summary['n_windows_out_of_band']:>5}  "
            f"polyT {summary['n_polyT_runs']:>3}  "
            f"repeats {summary['n_tandem_repeats']:>3}  "
            f"stops {summary['n_stops_total']:>4}  "
            f"-> {summary['flags'] or 'no flags'}",
            file=sys.stderr,
        )
    if mismatched:
        print(
            "warning: the selection table disagrees with itself; the actual "
            "span was used:\n         " + "\n         ".join(mismatched),
            file=sys.stderr,
        )
    print(
        f"wrote {summary_path}, {features_path}, {report_path}",
        file=sys.stderr,
    )
    return 0


def analyse(
    region: str, accession: str, gene: str, offset: int, args: argparse.Namespace
) -> tuple[dict, list[dict]]:
    """Measure one region. ``offset`` is its 1-based start on the transcript."""
    features: list[dict] = []

    def emit(kind: str, start: int, end: int, detail: str, severity: str) -> None:
        """Record a feature at region offsets ``start``-``end``, 1-based."""
        features.append({
            "transcript": accession,
            "gene": gene,
            "feature_type": kind,
            # Both frames of reference, deliberately: see the module docstring.
            "start": offset + start - 1,
            "end": offset + end - 1,
            "offset_in_region": start,
            "length": end - start + 1,
            "detail": detail,
            "severity": severity,
        })

    gc, ambiguous = seqqc.gc_fraction(region)
    gc_windows = seqqc.sliding_gc(region, args.window)
    low = [w for w in gc_windows if w.value < args.gc_min]
    high = [w for w in gc_windows if w.value > args.gc_max]
    out_of_band = [w.value < args.gc_min or w.value > args.gc_max for w in gc_windows]
    for window in low:
        emit("gc_window_low", window.start, window.end,
             f"GC {window.value:.3f} < {args.gc_min}", "warn")
    for window in high:
        emit("gc_window_high", window.start, window.end,
             f"GC {window.value:.3f} > {args.gc_max}", "warn")

    polyt = seqqc.homopolymer_runs(region, args.polyt_min, bases="T")
    for run in polyt:
        emit("polyT", run.start, run.end, "T" * run.length, "warn")

    # T excluded: it is already reported above, under its own threshold.
    homopolymers = seqqc.homopolymer_runs(region, args.homopolymer_min, bases="ACG")
    for run in homopolymers:
        emit("homopolymer", run.start, run.end,
             f"{run.base} x{run.length}", "warn")

    repeats = seqqc.tandem_repeats(region, min_copies=args.repeat_copies)
    for repeat in repeats:
        emit("tandem_repeat", repeat.start, repeat.end,
             f"({repeat.motif}) x{repeat.copies}", "warn")

    entropy_windows = seqqc.sliding_entropy(region, args.window)
    low_entropy = [w for w in entropy_windows if w.value < args.min_entropy]
    for window in low_entropy:
        emit("low_entropy_window", window.start, window.end,
             f"{window.value:.3f} bits < {args.min_entropy}", "warn")

    stops = seqqc.stop_codons(region)
    for stop in stops:
        emit("stop_codon", stop.start, stop.end, stop.codon, "info")

    gc_values = [w.value for w in gc_windows]
    entropy_values = [w.value for w in entropy_windows]
    longest_polyt = max(polyt, key=lambda r: r.length, default=None)
    longest_homopolymer = max(homopolymers, key=lambda r: r.length, default=None)
    longest_repeat = max(repeats, key=lambda r: r.length, default=None)
    min_entropy_window = min(entropy_windows, key=lambda w: w.value, default=None)

    # Stop codons are deliberately absent: see the module docstring.
    flags = []
    if low or high:
        flags.append("gc_out_of_band")
    if polyt:
        flags.append("polyT")
    if homopolymers:
        flags.append("homopolymer")
    if repeats:
        flags.append("tandem_repeat")
    if low_entropy:
        flags.append("low_entropy")

    summary = {
        "transcript": accession,
        "gene": gene,
        "start": offset,
        "end": offset + len(region) - 1,
        "n_ambiguous": ambiguous,
        "gc_frac": round(gc, 4),
        "gc_min_window": round(min(gc_values), 4) if gc_values else "",
        "gc_max_window": round(max(gc_values), 4) if gc_values else "",
        "gc_mean_window": round(statistics.fmean(gc_values), 4) if gc_values else "",
        "gc_sd_window": round(statistics.pstdev(gc_values), 4)
        if len(gc_values) > 1 else 0.0,
        "n_windows_out_of_band": len(low) + len(high),
        "widest_out_of_band_run": seqqc.widest_run(out_of_band),
        "gc_band": f"{args.gc_min}-{args.gc_max}",
        "n_polyT_runs": len(polyt),
        "longest_polyT": longest_polyt.length if longest_polyt else 0,
        "longest_polyT_at": offset + longest_polyt.start - 1 if longest_polyt else "",
        "polyT_min": args.polyt_min,
        "n_homopolymers": len(homopolymers),
        "longest_homopolymer": longest_homopolymer.length if longest_homopolymer else 0,
        "longest_homopolymer_base": longest_homopolymer.base
        if longest_homopolymer else "",
        "homopolymer_min": args.homopolymer_min,
        "n_tandem_repeats": len(repeats),
        "longest_tandem_repeat": longest_repeat.length if longest_repeat else 0,
        "longest_tandem_motif": longest_repeat.motif if longest_repeat else "",
        "min_window_entropy": round(min(entropy_values), 4) if entropy_values else "",
        "min_entropy_at": offset + min_entropy_window.start - 1
        if min_entropy_window else "",
        "n_windows_low_entropy": len(low_entropy),
        "entropy_min": args.min_entropy,
        "n_stop_TAA": sum(1 for s in stops if s.codon == "TAA"),
        "n_stop_TAG": sum(1 for s in stops if s.codon == "TAG"),
        "n_stop_TGA": sum(1 for s in stops if s.codon == "TGA"),
        "n_stops_total": len(stops),
        "stops_per_100nt": round(100 * len(stops) / len(region), 2) if region else 0.0,
        "flags": ";".join(flags),
        "n_flags": len(flags),
    }
    return summary, features


def write_tsv(path: str, columns: Sequence[str], rows: list[dict]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as handle:
        # lineterminator="\n" to match every other table this pipeline writes.
        # csv defaults to "\r\n", which leaves a stray CR on the last column of
        # every row -- enough to turn a header field into "severity\r" and
        # break any consumer keying on the column name.
        writer = csv.DictWriter(
            handle,
            fieldnames=list(columns),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_report(
    path: str, summaries: list[dict], features: list[dict], args: argparse.Namespace
) -> None:
    lines = [
        "# Target site QC",
        "",
        f"Window {args.window} nt, GC band {args.gc_min}-{args.gc_max}, "
        f"poly-T >= {args.polyt_min}, homopolymer >= {args.homopolymer_min}, "
        f"tandem repeat >= {args.repeat_copies} copies, "
        f"entropy < {args.min_entropy} bits.",
        "",
        "Stop codons are counted, not judged: chance alone puts about 94 stop "
        "triplets in 2 kb, so their presence is reported and never flagged.",
        "",
        "| gene | region | nt | GC | windows out of band | poly-T | homopolymers "
        "| repeats | stops | flags |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for s in summaries:
        lines.append(
            f"| {s['gene']} | {s['start']:,}-{s['end']:,} | {s['length_actual']:,} "
            f"| {s['gc_frac']:.3f} | {s['n_windows_out_of_band']} "
            f"| {s['n_polyT_runs']} | {s['n_homopolymers']} "
            f"| {s['n_tandem_repeats']} | {s['n_stops_total']} "
            f"| {s['flags'] or '-'} |"
        )

    for s in summaries:
        gene = s["gene"]
        mine = [f for f in features if f["gene"] == gene]
        lines += [
            "",
            f"## {gene} ({s['transcript']}) {s['start']:,}-{s['end']:,}",
            "",
            f"- **Length**: {s['length_actual']:,} nt actual, "
            f"{s['length_declared']:,} declared"
            + ("" if s["length_matches"] else " -- **these disagree**"),
            f"- **GC**: {s['gc_frac']:.1%} overall; windows "
            f"{s['gc_min_window']:.1%}-{s['gc_max_window']:.1%}, "
            f"mean {s['gc_mean_window']:.1%}, SD {s['gc_sd_window']:.3f}. "
            f"{s['n_windows_out_of_band']} window(s) outside {s['gc_band']}, "
            f"widest consecutive run {s['widest_out_of_band_run']}.",
            f"- **Poly-T**: {s['n_polyT_runs']} run(s) of >= {s['polyT_min']} T"
            + (f", longest {s['longest_polyT']} nt at {s['longest_polyT_at']:,}."
               if s["n_polyT_runs"] else "."),
            f"- **Repeats**: {s['n_homopolymers']} homopolymer(s) of >= "
            f"{s['homopolymer_min']}"
            + (f" (longest {s['longest_homopolymer']} x "
               f"{s['longest_homopolymer_base']})" if s["n_homopolymers"] else "")
            + f"; {s['n_tandem_repeats']} tandem repeat(s)"
            + (f" (longest {s['longest_tandem_repeat']} nt, motif "
               f"{s['longest_tandem_motif']})" if s["n_tandem_repeats"] else "")
            + f"; lowest window entropy {s['min_window_entropy']:.3f} bits at "
            f"{s['min_entropy_at']:,}, {s['n_windows_low_entropy']} window(s) "
            f"below {s['entropy_min']}.",
            f"- **Stop codons**: {s['n_stops_total']} total "
            f"({s['n_stop_TAA']} TAA, {s['n_stop_TAG']} TAG, "
            f"{s['n_stop_TGA']} TGA), {s['stops_per_100nt']} per 100 nt. "
            "Presence only -- no frame, not a defect.",
            f"- **Flags**: {s['flags'] or 'none'} "
            f"({len([f for f in mine if f['severity'] == 'warn'])} flagged "
            f"feature(s) in the features file)",
        ]

    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
