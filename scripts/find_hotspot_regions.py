#!/usr/bin/env python3
"""Find the 1.5-2 kb sub-region of each transcript holding the most guide hotspots.

    python scripts/find_hotspot_regions.py \
        --predictions results/pik3/pik3-predictions.csv \
        --regions results/blastn/PIK3C*_location.blastn \
        --output results/pik3/hotspot_regions.tsv

Stages:

1. Read the predictions CSV and collapse every run of guides scoring at or above
   ``--min-efficacy`` into a **hotspot** -- one independent target site. At step-1
   tiling a single good site produces ~30 overlapping high-efficacy guides, so
   counting guides would count that site thirty times; see src/hotspots.py.
2. Read the conserved regions produced by scripts/get_conserved_regions.py and
   slide windows of ``--min-window``..``--max-window`` bases inside each one. A
   hotspot counts towards a window only if it lies wholly inside it.
3. Rank the windows per transcript -- most hotspots first, ties to the shorter
   and more evenly covered window -- suppress overlapping near-duplicates, and
   write the top ``--top`` candidates per transcript.

Windows are confined to the conserved regions. 

Two parameters dominate the result and are echoed into every row of the output,
because a hotspot count means nothing without them:

``--min-efficacy``
    Absolute, not a quantile, so counts stay comparable between transcripts --
    and so it must be set against the score distribution of the predictor that
    produced the file. The default 0.9 suits the committed PI3K predictions
    (median 0.823, so 0.8 would qualify over half of every transcript and
    hotspots would merge into kilobase-wide blobs).
``--merge-gap``
    How many low-scoring bases may sit between two clusters and still count as
    one site. The most sensitive knob here: on the PI3K run, 0 gives 47-87
    hotspots per transcript and 30 gives 8-22. The default 0 merges only
    clusters that overlap or abut.

Use ``--sweep-thresholds`` to check whether a winning window is stable or an
artefact of one cutoff, and ``--all-windows``/``--hotspots-out`` to inspect what
the ranking was built from.
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

import pandas as pd  # noqa: E402

from hotspots import (  # noqa: E402
    DEFAULT_MAX_WINDOW,
    DEFAULT_MERGE_GAP,
    DEFAULT_MIN_EFFICACY,
    DEFAULT_MIN_WINDOW,
    DEFAULT_SIZE_STEP,
    DEFAULT_STEP,
    PREDICTION_COLUMNS,
    REGION_COLUMNS,
    filter_by_core,
    find_hotspots,
    hotspots_to_frame,
    oversized_hotspots,
    rank_windows,
    regions_from_table,
    scan_regions,
    select_non_overlapping,
    sweep_thresholds,
    unplaceable_hotspots,
    whole_transcript_regions,
    windows_to_frame,
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--predictions",
        required=True,
        metavar="PATH",
        help="predictions CSV with at least "
        f"{', '.join(PREDICTION_COLUMNS)} columns "
        "(e.g. results/pik3/pik3-predictions.csv)",
    )
    parser.add_argument(
        "--regions",
        nargs="+",
        default=None,
        metavar="PATH",
        help="one or more region tables from scripts/get_conserved_regions.py, "
        f"with at least {', '.join(REGION_COLUMNS)} columns. Windows are "
        "confined to these spans. Omit to scan each transcript end to end, "
        "which gives up the conservation guarantee.",
    )
    parser.add_argument(
        "--output",
        required=True,
        metavar="PATH",
        help="TSV of the top --top non-overlapping windows per transcript, "
        "best first",
    )
    parser.add_argument(
        "--min-efficacy",
        type=float,
        default=DEFAULT_MIN_EFFICACY,
        metavar="X",
        help="a guide at or above this predicted efficacy can form a hotspot "
        "(default: %(default)s). Absolute, not a quantile: check it against "
        "your predictor's score distribution.",
    )
    parser.add_argument(
        "--merge-gap",
        type=int,
        default=DEFAULT_MERGE_GAP,
        metavar="N",
        help="low-scoring bases allowed between two clusters of the same "
        "hotspot (default: %(default)s, i.e. merge only when they overlap or "
        "abut). The most sensitive parameter here.",
    )
    parser.add_argument(
        "--min-window",
        type=int,
        default=DEFAULT_MIN_WINDOW,
        metavar="N",
        help="shortest window length to consider (default: %(default)s)",
    )
    parser.add_argument(
        "--max-window",
        type=int,
        default=DEFAULT_MAX_WINDOW,
        metavar="N",
        help="longest window length to consider (default: %(default)s)",
    )
    parser.add_argument(
        "--size-step",
        type=int,
        default=DEFAULT_SIZE_STEP,
        metavar="N",
        help="increment between window lengths in the sweep (default: "
        "%(default)s). The longest length that fits a region is always "
        "included regardless.",
    )
    parser.add_argument(
        "--step",
        type=int,
        default=DEFAULT_STEP,
        metavar="N",
        help="bases between successive window start positions (default: "
        "%(default)s)",
    )
    parser.add_argument(
        "--max-gap",
        type=int,
        default=None,
        metavar="N",
        help="reject any window containing a stretch of N+ bases with no "
        "hotspot in it (default: no limit). Hotspot count alone cannot tell "
        "many sites spread evenly from the same number clumped at one end, and "
        "a dead zone is where no guide can be placed at all; this makes even "
        "coverage a requirement rather than a tie-break.",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=5,
        metavar="N",
        help="non-overlapping windows to keep per transcript (default: "
        "%(default)s; 0 keeps all of them)",
    )
    parser.add_argument(
        "--all-windows",
        metavar="PATH",
        help="also write every scored window, ranked, before overlap "
        "suppression -- the table to reach for when a winner looks surprising",
    )
    parser.add_argument(
        "--hotspots-out",
        metavar="PATH",
        help="also write every hotspot found, with its guide count and best "
        "guide, so the windows can be audited against their inputs",
    )
    parser.add_argument(
        "--sweep-thresholds",
        metavar="X,Y,Z",
        help="comma-separated efficacy thresholds; writes the best window per "
        "transcript at each one, to show whether a winner is stable or an "
        "artefact of a single cutoff",
    )
    parser.add_argument(
        "--sweep-out",
        metavar="PATH",
        help="where to write the sweep table (default: --output with a "
        ".sweep.tsv suffix)",
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


def read_regions(paths: list[str]):
    """Read and concatenate region tables, keeping their source for the summary."""
    frames = []
    for path in paths:
        # Split on any run of whitespace, as plot_scatterlines.py does: these
        # tables are often hand-filtered into a mix of tabs and spaces that a
        # strict sep="\t" folds into one column, giving a baffling "missing
        # column" error on a file that looks fine.
        frame = pd.read_csv(path, sep=r"\s+", header=0, engine="python")
        require_columns(frame, REGION_COLUMNS, path)
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def write_tsv(frame: pd.DataFrame, path: str) -> None:
    parent = os.path.dirname(os.path.abspath(path))
    os.makedirs(parent, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    if args.min_window > args.max_window:
        raise SystemExit(
            f"error: --min-window {args.min_window} exceeds --max-window "
            f"{args.max_window}"
        )
    for flag, value in (("--size-step", args.size_step), ("--step", args.step)):
        if value < 1:
            raise SystemExit(f"error: {flag} must be at least 1, got {value}")

    # Parsed up here rather than at the point of use: the sweep runs last, and
    # rejecting a typo only after the scan has finished and written its output
    # wastes the run and leaves a half-finished set of files behind.
    sweep_list: list[float] = []
    if args.sweep_thresholds:
        try:
            sweep_list = [
                float(x) for x in args.sweep_thresholds.split(",") if x.strip()
            ]
        except ValueError:
            raise SystemExit(
                f"error: --sweep-thresholds {args.sweep_thresholds!r} is not a "
                "comma-separated list of numbers"
            )
        if not sweep_list:
            raise SystemExit("error: --sweep-thresholds is empty")

    # Checked before pandas opens them: a mistyped path is the likeliest thing
    # to go wrong here, and pandas answers it with a bare traceback.
    paths = [("--predictions", args.predictions)]
    paths += [("--regions", p) for p in (args.regions or [])]
    for flag, path in paths:
        if not os.path.isfile(path):
            raise SystemExit(f"error: {flag} {path}: no such file")

    predictions = pd.read_csv(args.predictions, header=0)
    require_columns(predictions, PREDICTION_COLUMNS, args.predictions)

    hotspots = find_hotspots(predictions, args.min_efficacy, args.merge_gap)
    if not hotspots:
        raise SystemExit(
            f"error: no guide in {args.predictions} scores at or above "
            f"--min-efficacy {args.min_efficacy}, so there are no hotspots to "
            f"place. The highest predicted efficacy in the file is "
            f"{predictions['predicted_efficacy'].max():.4f}."
        )

    if args.regions:
        regions = regions_from_table(read_regions(args.regions))
        if not regions:
            raise SystemExit(
                "error: the --regions table(s) contain no placed region "
                "(every row has empty ref_start/ref_end)"
            )
    else:
        regions = whole_transcript_regions(predictions)

    windows = scan_regions(
        hotspots,
        regions,
        min_window=args.min_window,
        max_window=args.max_window,
        size_step=args.size_step,
        step=args.step,
    )
    if not windows:
        raise SystemExit(
            "error: no candidate window could be placed. Every region is "
            "shorter than it should be, or the regions and the predictions use "
            "different coordinates."
        )

    # Applied before ranking, so a gappy window cannot win on count alone.
    # Reported per transcript rather than as one total: a limit that empties one
    # transcript while leaving the others untouched is the useful thing to see.
    kept = filter_by_core(windows, args.max_gap)
    if args.max_gap is not None and not kept:
        raise SystemExit(
            f"error: --max-gap {args.max_gap} rejected every one of the "
            f"{len(windows):,} candidate window(s). The most evenly covered "
            f"window anywhere has a "
            f"{min(w.largest_gap for w in windows):,} bp dead zone."
        )
    # `windows` stays the unfiltered scan: --all-windows exists to explain a
    # surprising winner, and a window rejected by --max-gap is exactly the kind
    # of surprise worth seeing. Only the ranking works from the filtered set.
    best = select_non_overlapping(kept, top=args.top)
    write_tsv(windows_to_frame(best), args.output)

    # ---- summary -------------------------------------------------------
    print(
        f"{len(predictions):,} guides -> {len(hotspots):,} hotspot(s) at "
        f"--min-efficacy {args.min_efficacy} --merge-gap {args.merge_gap}",
        file=sys.stderr,
    )
    if not args.regions:
        print(
            "no --regions given: scanning each transcript end to end, so a "
            "winning window is NOT guaranteed to be conserved across variants",
            file=sys.stderr,
        )
    print(
        f"{len(regions)} region(s) over "
        f"{len({r.transcript for r in regions})} transcript(s), "
        f"{len(windows):,} candidate window(s) scored",
        file=sys.stderr,
    )

    # Reported per transcript as well as in total: a limit that empties one
    # transcript while leaving the others untouched is the useful thing to see.
    if args.max_gap is not None:
        print(
            f"--max-gap {args.max_gap}: {len(windows) - len(kept):,} of "
            f"{len(windows):,} window(s) rejected for a longer dead zone, "
            f"{len(kept):,} ranked",
            file=sys.stderr,
        )
        for transcript in sorted(
            {w.transcript for w in windows} - {w.transcript for w in kept}
        ):
            tightest = min(
                w.largest_gap for w in windows if w.transcript == transcript
            )
            print(
                f"warning: no window on {transcript} satisfies --max-gap "
                f"{args.max_gap}; its most evenly covered window has a "
                f"{tightest:,} bp dead zone, so that transcript is absent from "
                f"the output",
                file=sys.stderr,
            )

    # A region too short for --min-window still yields one window, flagged
    # undersized; say so, because its hotspot count is not comparable with a
    # full-length window's.
    undersized = sorted(
        {(w.transcript, w.region_start, w.region_end, w.length)
         for w in windows if w.undersized}
    )
    if undersized:
        print(
            f"{len(undersized)} region(s) shorter than --min-window "
            f"{args.min_window}, each scanned as a single short window -- their "
            f"counts are not comparable with full-length windows:",
            file=sys.stderr,
        )
        for transcript, rs, re_, length in undersized:
            print(f"  {transcript}:{rs}-{re_}  {length:,} bp", file=sys.stderr)

    # Hotspots no window can contain. Loud, because a region set that excludes
    # most of them usually means the two files disagree about coordinates.
    unplaceable = unplaceable_hotspots(hotspots, regions)
    if unplaceable:
        share = len(unplaceable) / len(hotspots)
        print(
            f"{len(unplaceable):,} of {len(hotspots):,} hotspot(s) "
            f"({share:.0%}) lie outside every region or straddle a boundary, so "
            f"no window contains them",
            file=sys.stderr,
        )
        if share > 0.5:
            print(
                "warning: that is more than half of them. Check that the "
                "region coordinates and the predictions describe the same "
                "transcripts in the same coordinate system.",
                file=sys.stderr,
            )
    oversized = oversized_hotspots(hotspots, args.max_window)
    if oversized:
        print(
            f"{len(oversized)} hotspot(s) are longer than --max-window "
            f"{args.max_window} and so can never be contained by any window; "
            f"the widest is "
            f"{max(h.length for h in oversized):,} bp. A hotspot this wide "
            f"usually means --min-efficacy is too low for this predictor.",
            file=sys.stderr,
        )

    print(f"top window(s) per transcript -> {args.output}", file=sys.stderr)
    shown: set[str] = set()
    for w in rank_windows(best):
        if w.transcript in shown:
            continue
        shown.add(w.transcript)
        total = sum(1 for h in hotspots if h.transcript == w.transcript)
        print(
            f"  {w.transcript}: {w.start:,}-{w.end:,} ({w.length:,} bp"
            f"{', undersized' if w.undersized else ''}) "
            f"{w.n_hotspots} hotspot(s) of {total} on the transcript, "
            f"{w.n_guides:,} guide(s), largest gap {w.largest_gap:,} bp, "
            f"best efficacy {w.max_efficacy:.4f}",
            file=sys.stderr,
        )

    if args.all_windows:
        write_tsv(windows_to_frame(windows), args.all_windows)
        print(
            f"{len(windows):,} scored window(s), ranked -> {args.all_windows}",
            file=sys.stderr,
        )
    if args.hotspots_out:
        write_tsv(hotspots_to_frame(hotspots), args.hotspots_out)
        print(
            f"{len(hotspots):,} hotspot(s) -> {args.hotspots_out}",
            file=sys.stderr,
        )

    if sweep_list:
        sweep_out = args.sweep_out or f"{os.path.splitext(args.output)[0]}.sweep.tsv"
        sweep = sweep_thresholds(
            predictions,
            sweep_list,
            merge_gap=args.merge_gap,
            regions=regions,
            max_gap=args.max_gap,
            min_window=args.min_window,
            max_window=args.max_window,
            size_step=args.size_step,
            step=args.step,
        )
        write_tsv(sweep, sweep_out)
        print(
            f"threshold sweep over {', '.join(str(t) for t in sweep_list)} "
            f"-> {sweep_out}",
            file=sys.stderr,
        )
        # A window that wins at one threshold only is not a robust choice --
        # but a 50 bp shift and a jump to the other end of the transcript are
        # very different findings, so say which happened. Disjoint winners mean
        # the threshold, not the data, is choosing the region.
        #
        # How many thresholds actually produced a window matters as much as
        # whether they agreed. --max-gap can reject every window at the higher
        # thresholds, leaving a transcript with a single row -- and one row
        # agrees with itself, so calling that "stable" would be a false
        # reassurance built on missing evidence.
        for transcript, group in sweep.groupby("transcript_id"):
            spans = sorted(
                {(int(r.window_start), int(r.window_end)) for r in group.itertuples()}
            )
            tested = sorted(float(t) for t in group["min_efficacy"].unique())
            missing = [t for t in sweep_list if t not in tested]
            if missing:
                why = (
                    f" (--max-gap {args.max_gap} left no window there)"
                    if args.max_gap is not None
                    else ""
                )
                print(
                    f"  {transcript}: only {len(tested)} of {len(sweep_list)} "
                    f"threshold(s) produced a window -- nothing at "
                    f"{', '.join(str(t) for t in missing)}{why}. Stability is "
                    f"UNTESTED here, not confirmed.",
                    file=sys.stderr,
                )
                continue
            if len(spans) < 2:
                print(
                    f"  {transcript}: the same window wins at all "
                    f"{len(sweep_list)} threshold(s) -- a stable choice",
                    file=sys.stderr,
                )
                continue
            disjoint = [
                (a, b)
                for i, a in enumerate(spans)
                for b in spans[i + 1:]
                if a[1] < b[0] or b[1] < a[0]
            ]
            if disjoint:
                print(
                    f"  {transcript}: the best window JUMPS between "
                    f"non-overlapping spans as the threshold changes "
                    f"({'; '.join(f'{s}-{e}' for s, e in spans)}). The cutoff, "
                    f"not the data, is picking the region -- compare the "
                    f"candidates in {args.output} before committing.",
                    file=sys.stderr,
                )
            else:
                print(
                    f"  {transcript}: the best window shifts but stays on one "
                    f"overlapping span across all {len(sweep_list)} threshold(s) "
                    f"-- a stable choice",
                    file=sys.stderr,
                )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
