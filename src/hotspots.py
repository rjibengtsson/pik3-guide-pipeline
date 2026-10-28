"""Find the sub-region of a transcript holding the most guide hotspots.

Window is scored on guide count ranks the
*longest good stretch* first rather than the region offering the most
independent choices. This module collapses overlapping high-efficacy guides into
**hotspots** -- one hotspot is one site -- and scores a window by how many whole
hotspots it contains.

Definitions
-----------
hotspot
    A maximal run of guides scoring at or above ``min_efficacy``, merged
    whenever consecutive runs are separated by no more than ``merge_gap``
    unqualified bases. At the default ``merge_gap=0`` two clusters merge only
    when they overlap or abut, so guides sharing no base stay separate sites.
window
    A candidate construct: a span inside one conserved region. A hotspot counts
    towards a window only if it lies **wholly** inside it -- a hotspot clipped
    by the window edge is not a usable site in the cloned fragment.


Choosing ``min_efficacy``
-------------------------
The threshold is **absolute, not a per-transcript quantile**, so counts stay
comparable between transcripts.

Hotspots are found over the whole transcript and only then tested against the
regions, so the hotspot boundaries do not depend on which regions were supplied
and are stable as the region set is edited. The consequence is that a hotspot
straddling a region boundary belongs to no window; :func:`scan_regions` reports
those rather than dropping them silently.

Window length is swept, not fixed, because conserved regions are often shorter
than a nominal 2 kb. Longer windows trivially contain more hotspots, so raw
count alone would always rank the longest window first; :func:`rank_windows`
breaks ties towards the **shorter** window, which makes the sweep meaningful --
the answer to "most hotspots" becomes the smallest construct achieving that
count.

Conventions follow src/tile.py: coordinates are 1-based inclusive, so a span
covers ``end - start + 1`` bases and ``sequence[start - 1:end]`` is it.

This is the library half of the stage: it parses no arguments and reads no file
by a hard-coded path. Every function takes a DataFrame or a sequence of records,
so the merging and scoring rules are testable without any I/O.
``scripts/find_hotspot_regions.py`` is the CLI around it.
"""

from __future__ import annotations

from typing import Iterable, Iterator, NamedTuple, Sequence

import pandas as pd

#: Columns required of a predictions CSV (e.g. results/pik3/pik3-predictions.csv).
PREDICTION_COLUMNS = ("transcript_id", "start", "end", "predicted_efficacy")

#: Columns required of a region table -- the output of
#: scripts/get_conserved_regions.py, whose ``ref_start``/``ref_end`` are already
#: 1-based inclusive coordinates on the ungapped reference.
REGION_COLUMNS = ("reference", "ref_start", "ref_end")

DEFAULT_MIN_EFFICACY = 0.9
DEFAULT_MERGE_GAP = 0
DEFAULT_MIN_WINDOW = 1500
DEFAULT_MAX_WINDOW = 2000
DEFAULT_SIZE_STEP = 250
DEFAULT_STEP = 50


class Hotspot(NamedTuple):
    """One independent high-efficacy target site."""

    transcript: str
    start: int
    end: int
    n_guides: int
    max_efficacy: float
    mean_efficacy: float
    best_guide: str

    @property
    def length(self) -> int:
        return self.end - self.start + 1


class Region(NamedTuple):
    """A span a window is allowed to sit inside -- typically conserved."""

    transcript: str
    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start + 1


class Window(NamedTuple):
    """A scored candidate construct."""

    transcript: str
    region_start: int
    region_end: int
    start: int
    end: int
    length: int
    undersized: bool
    n_hotspots: int
    n_guides: int
    hotspots_per_kb: float
    max_efficacy: float
    mean_hotspot_efficacy: float
    covered_bp: int
    frac_covered: float
    largest_gap: int
    hotspot_loci: str


def merge_intervals(
    intervals: Iterable[tuple[int, int]], gap: int = DEFAULT_MERGE_GAP
) -> list[tuple[int, int]]:
    """Merge 1-based inclusive spans separated by at most ``gap`` bases.

    ``gap`` counts the bases strictly *between* two spans, so ``gap=0`` merges
    spans that overlap or abut and leaves spans sharing no base apart. Negative
    values are accepted and require genuine overlap.
    """
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if merged and start - merged[-1][1] - 1 <= gap:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def find_hotspots(
    predictions: pd.DataFrame,
    min_efficacy: float = DEFAULT_MIN_EFFICACY,
    merge_gap: int = DEFAULT_MERGE_GAP,
) -> list[Hotspot]:
    """Collapse high-efficacy guides into hotspots, one record per site.

    Returns hotspots sorted by ``(transcript, start)``. Guides below
    ``min_efficacy`` take no part: they neither form nor extend a hotspot, so
    ``n_guides`` on a hotspot counts qualifying guides only.
    """
    hotspots: list[Hotspot] = []
    qualifying = predictions[predictions["predicted_efficacy"] >= min_efficacy]

    for transcript, group in qualifying.groupby("transcript_id", sort=True):
        # Sorted by start so a single pass can close each cluster: the next row
        # can only ever extend the current one or begin a new one.
        rows = group.sort_values(["start", "end"]).itertuples(index=False)
        cluster: list = []
        cluster_end = 0

        def flush(cluster: list, cluster_end: int) -> None:
            efficacies = [r.predicted_efficacy for r in cluster]
            best = max(cluster, key=lambda r: r.predicted_efficacy)
            hotspots.append(
                Hotspot(
                    transcript=str(transcript),
                    start=int(cluster[0].start),
                    end=int(cluster_end),
                    n_guides=len(cluster),
                    max_efficacy=float(max(efficacies)),
                    mean_efficacy=float(sum(efficacies) / len(efficacies)),
                    # tile_id is optional in the input, so fall back to the
                    # coordinates rather than failing on a trimmed file.
                    best_guide=str(
                        getattr(best, "tile_id", f"{transcript}:{best.start}-{best.end}")
                    ),
                )
            )

        for row in rows:
            if cluster and int(row.start) - cluster_end - 1 <= merge_gap:
                cluster.append(row)
                cluster_end = max(cluster_end, int(row.end))
            else:
                if cluster:
                    flush(cluster, cluster_end)
                cluster = [row]
                cluster_end = int(row.end)
        if cluster:
            flush(cluster, cluster_end)

    return sorted(hotspots, key=lambda h: (h.transcript, h.start, h.end))


def regions_from_table(regions: pd.DataFrame) -> list[Region]:
    """Build regions from a get_conserved_regions.py table.

    A row whose hit was reported on the minus strand carries ``ref_start >
    ref_end``; the span is ordered here so nothing downstream has to care. Rows
    with no hit (empty coordinates) are skipped -- that is how
    get_conserved_regions.py records a region it could not place, and it has no
    position to scan.
    """
    out: list[Region] = []
    for row in regions.itertuples(index=False):
        if pd.isna(row.ref_start) or pd.isna(row.ref_end):
            continue
        start, end = int(row.ref_start), int(row.ref_end)
        out.append(
            Region(
                transcript=str(row.reference),
                start=min(start, end),
                end=max(start, end),
            )
        )
    return sorted(set(out), key=lambda r: (r.transcript, r.start, r.end))


def whole_transcript_regions(predictions: pd.DataFrame) -> list[Region]:
    """One region per transcript, spanning the full tiling.

    The fallback when no conserved regions are supplied. The end is the largest
    tile ``end``, which is exact for a step-1 tiling of the whole transcript but
    is a **floor**, not a measurement -- an already-filtered predictions file
    stops at its last surviving guide.
    """
    bounds = predictions.groupby("transcript_id")["end"].max()
    return [
        Region(transcript=str(t), start=1, end=int(e)) for t, e in sorted(bounds.items())
    ]


def window_bounds(
    region: Region,
    min_window: int = DEFAULT_MIN_WINDOW,
    max_window: int = DEFAULT_MAX_WINDOW,
    size_step: int = DEFAULT_SIZE_STEP,
    step: int = DEFAULT_STEP,
) -> Iterator[tuple[int, int, bool]]:
    """Yield ``(start, end, undersized)`` for every candidate window in a region.

    Window lengths are swept from ``min_window`` to ``max_window`` in
    ``size_step`` increments, capped at the region. The largest length that fits
    is always included even when ``size_step`` would step over it, so the
    widest usable construct is never missed.

    A region shorter than ``min_window`` yields exactly one window spanning the
    whole region, flagged ``undersized``: these are reported rather than skipped
    because a short conserved region is still a candidate -- it simply cannot
    reach the requested length, and silently dropping it would hide the region
    from the ranking entirely.
    """
    if region.length < min_window:
        yield region.start, region.end, True
        return

    largest = min(max_window, region.length)
    sizes = sorted(set(range(min_window, largest + 1, size_step)) | {largest})

    for size in sizes:
        last_start = region.end - size + 1
        starts = list(range(region.start, last_start + 1, step))
        # Flush right as well as left: otherwise a region whose length is not a
        # multiple of `step` never has its tail scanned.
        if starts and starts[-1] != last_start:
            starts.append(last_start)
        for start in starts:
            yield start, start + size - 1, False


def score_window(
    transcript: str,
    region: Region,
    start: int,
    end: int,
    hotspots: Sequence[Hotspot],
    undersized: bool = False,
) -> Window:
    """Score one window against the hotspots of its transcript.

    ``hotspots`` must already be restricted to ``transcript`` and sorted by
    start. Only wholly-contained hotspots count.
    """
    inside = [h for h in hotspots if h.start >= start and h.end <= end]
    length = end - start + 1
    covered = sum(h.length for h in inside)

    # Largest run of window with no hotspot in it, edges included: a window
    # whose hotspots all huddle at one end scores the same count as an evenly
    # covered one, and this is what separates them.
    largest_gap = 0
    cursor = start
    for h in inside:
        largest_gap = max(largest_gap, h.start - cursor)
        cursor = h.end + 1
    largest_gap = max(largest_gap, end - cursor + 1)

    return Window(
        transcript=transcript,
        region_start=region.start,
        region_end=region.end,
        start=start,
        end=end,
        length=length,
        undersized=undersized,
        n_hotspots=len(inside),
        n_guides=sum(h.n_guides for h in inside),
        hotspots_per_kb=round(len(inside) * 1000 / length, 3),
        max_efficacy=round(max((h.max_efficacy for h in inside), default=0.0), 6),
        mean_hotspot_efficacy=round(
            sum(h.max_efficacy for h in inside) / len(inside), 6
        )
        if inside
        else 0.0,
        covered_bp=covered,
        frac_covered=round(covered / length, 4),
        largest_gap=largest_gap,
        # Locus strings follow the convention of classification.tsv: ";" between
        # loci, coordinates 1-based inclusive, the per-locus detail in brackets.
        hotspot_loci=";".join(f"{h.start}-{h.end}({h.n_guides})" for h in inside),
    )


def scan_regions(
    hotspots: Sequence[Hotspot],
    regions: Sequence[Region],
    min_window: int = DEFAULT_MIN_WINDOW,
    max_window: int = DEFAULT_MAX_WINDOW,
    size_step: int = DEFAULT_SIZE_STEP,
    step: int = DEFAULT_STEP,
) -> list[Window]:
    """Score every candidate window in every region."""
    by_transcript: dict[str, list[Hotspot]] = {}
    for h in hotspots:
        by_transcript.setdefault(h.transcript, []).append(h)

    windows: list[Window] = []
    for region in regions:
        transcript_hotspots = by_transcript.get(region.transcript, [])
        # Pre-filter to the region so the per-window scan is over a short list.
        in_region = [
            h
            for h in transcript_hotspots
            if h.start >= region.start and h.end <= region.end
        ]
        for start, end, undersized in window_bounds(
            region, min_window, max_window, size_step, step
        ):
            windows.append(
                score_window(
                    region.transcript, region, start, end, in_region, undersized
                )
            )
    return windows


def unplaceable_hotspots(
    hotspots: Sequence[Hotspot], regions: Sequence[Region]
) -> list[Hotspot]:
    """Hotspots no window can ever contain, so the loss is visible.

    A hotspot is unplaceable when it lies wholly inside no region -- it is
    outside every one, or it straddles a boundary. Reported, never dropped
    quietly: a region set that excludes most hotspots is a sign the regions and
    the predictions describe different coordinates.
    """
    return [
        h
        for h in hotspots
        if not any(
            r.transcript == h.transcript and h.start >= r.start and h.end <= r.end
            for r in regions
        )
    ]


def oversized_hotspots(
    hotspots: Sequence[Hotspot], max_window: int = DEFAULT_MAX_WINDOW
) -> list[Hotspot]:
    """Hotspots longer than the widest window, which therefore never count."""
    return [h for h in hotspots if h.length > max_window]


def filter_by_core(windows: Iterable[Window], max_gap: int | None) -> list[Window]:
    """Drop windows with a dead zone longer than ``max_gap``.

    Hotspot count alone does not distinguish many sites spread evenly from the
    same number clumped at one end with a long stretch of nothing between --
    and a construct is only as useful as its worst region, since a dead zone is
    where no guide can be placed at all. ``largest_gap`` is a tie-break in
    :func:`rank_windows`, which lets a gappy window still win on count; this
    makes evenness a hard requirement instead.

    ``max_gap=None`` keeps everything, which is the default: a sensible bound
    depends on how far apart guides may sit for the experiment, so it has to be
    chosen rather than assumed.
    """
    if max_gap is None:
        return list(windows)
    return [w for w in windows if w.largest_gap <= max_gap]


def rank_windows(windows: Iterable[Window]) -> list[Window]:
    """Sort best-first within each transcript.

    Most hotspots wins. Ties break towards the **shorter** window: a sweep over
    window lengths otherwise always favours the longest, since a longer span
    cannot contain fewer hotspots. Then the smaller largest gap, so the more
    evenly covered of two equal windows wins, then the start coordinate to make
    the order total and the output deterministic.
    """
    return sorted(
        windows,
        key=lambda w: (
            w.transcript,
            -w.n_hotspots,
            w.length,
            w.largest_gap,
            w.start,
        ),
    )


def select_non_overlapping(windows: Iterable[Window], top: int = 0) -> list[Window]:
    """Greedily keep the best windows that do not overlap, per transcript.

    Adjacent windows at ``step=50`` overlap by 97%, so an unsuppressed top-10 is
    ten shifts of one region. ``top=0`` keeps every non-overlapping window.
    """
    chosen: dict[str, list[Window]] = {}
    for w in rank_windows(windows):
        kept = chosen.setdefault(w.transcript, [])
        if top and len(kept) >= top:
            continue
        if all(w.end < c.start or w.start > c.end for c in kept):
            kept.append(w)
    return rank_windows(w for group in chosen.values() for w in group)


def sweep_thresholds(
    predictions: pd.DataFrame,
    thresholds: Sequence[float],
    merge_gap: int = DEFAULT_MERGE_GAP,
    regions: Sequence[Region] | None = None,
    max_gap: int | None = None,
    **window_kwargs,
) -> pd.DataFrame:
    """Best window per transcript at each threshold, to expose instability.

    The ranking is threshold-sensitive -- on the PI3K run PIK3CG's best window
    moves from 1101-3100 at 0.9 to 1-2000 at 0.95 -- so a winner that holds
    across thresholds is worth more than one that does not. One row per
    (threshold, transcript).

    ``max_gap`` must be the same value the main ranking used, or the sweep
    reports on a ranking the caller never saw: a window the filter rejects could
    be named the best at some threshold and make a stable choice look unstable,
    or the reverse. A transcript left with no window at a threshold simply has
    no row there.
    """
    rows = []
    for threshold in thresholds:
        hotspots = find_hotspots(predictions, threshold, merge_gap)
        scope = regions if regions is not None else whole_transcript_regions(predictions)
        best = select_non_overlapping(
            filter_by_core(scan_regions(hotspots, scope, **window_kwargs), max_gap),
            top=1,
        )
        for w in best:
            rows.append(
                {
                    "min_efficacy": threshold,
                    "transcript_id": w.transcript,
                    "n_hotspots_total": sum(
                        1 for h in hotspots if h.transcript == w.transcript
                    ),
                    "window_start": w.start,
                    "window_end": w.end,
                    "window_length": w.length,
                    "n_hotspots": w.n_hotspots,
                    "largest_gap": w.largest_gap,
                }
            )
    return pd.DataFrame(rows)


#: Field -> output column. A NamedTuple field is named for its record, but a
#: column has to say which coordinate it is once the tables are joined.
WINDOW_COLUMNS = {
    "transcript": "transcript_id",
    "start": "window_start",
    "end": "window_end",
    "length": "window_length",
}

HOTSPOT_COLUMNS = {
    "transcript": "transcript_id",
    "start": "hotspot_start",
    "end": "hotspot_end",
}


def windows_to_frame(windows: Sequence[Window]) -> pd.DataFrame:
    """Ranked windows as a DataFrame, with a per-transcript ``rank`` column."""
    ranked = rank_windows(windows)
    if not ranked:
        # Renamed here too: an empty result is still read by the same code as a
        # full one, so the two must not differ in their column names.
        return pd.DataFrame(
            columns=["rank"] + [WINDOW_COLUMNS.get(f, f) for f in Window._fields]
        )
    frame = pd.DataFrame([w._asdict() for w in ranked])
    frame.insert(0, "rank", frame.groupby("transcript").cumcount() + 1)
    return frame.rename(columns=WINDOW_COLUMNS)


def hotspots_to_frame(hotspots: Sequence[Hotspot]) -> pd.DataFrame:
    """All hotspots as a DataFrame, for inspecting what the windows were built from."""
    if not hotspots:
        return pd.DataFrame(
            columns=[HOTSPOT_COLUMNS.get(f, f) for f in Hotspot._fields] + ["length"]
        )
    frame = pd.DataFrame([h._asdict() for h in hotspots])
    frame["length"] = frame["end"] - frame["start"] + 1
    return frame.rename(columns=HOTSPOT_COLUMNS)
