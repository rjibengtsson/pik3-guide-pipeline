"""Sequence QC for a selected target region: GC, poly-T, repeats, stop codons.

Four properties, measured over a region already extracted from its transcript:

GC content
    Over the whole region, and over a sliding window the length of a spacer --
    the window **is** a candidate guide, so a region averaging 50% that holds a
    20% window is a different object from a uniformly 50% one, and only the
    window view separates them.
poly-T stretches
    Runs of T at or above a threshold. Four is the default because a run of
    four or more T is the classic premature-termination signal in a U6-driven
    construct.
repeats and homopolymers
    Homopolymers of A, C and G (T is covered by the poly-T scan, so it is not
    reported twice), short tandem repeats of a 2-4 nt motif, and a per-window
    Shannon entropy that catches the diffuse low-complexity a motif scan
    misses and that DUST would mask.
stop codons
    Every occurrence of TAA, TAG or TGA at every offset. **Presence only** --
    no reading frame, no CDS. Reported as counts and positions rather than as a
    pass/fail: in 2 kb of sequence chance alone predicts about 94 of them
    (3 triplets over 1998 offsets, 1 in 64), so a flag on "has stop codons"
    would fire on every region including the correct ones.

Offsets are **1-based inclusive and relative to the region**, so offset 1 is
the region's first base. The caller adds ``region_start - 1`` to put a feature
back on the transcript; keeping the two apart is what lets these functions be
tested on short literal strings.

Sequence is taken verbatim, in whatever alphabet the FASTA used, matching the
``target_seq`` convention in src/tile.py. Comparisons upper-case first, and
anything outside ACGT is counted as ambiguous and excluded from the GC
denominator rather than silently scored as not-GC.

This is library code: it parses no command line, reads no file, and prints
nothing. scripts/check_target_sites.py is the CLI around it.
"""

from __future__ import annotations

import math
import re
from typing import NamedTuple, Sequence

#: A spacer length: the GC and entropy windows are this wide, so one window is
#: one candidate guide rather than an arbitrary smoothing span.
DEFAULT_WINDOW = 30

DEFAULT_GC_MIN = 0.30
DEFAULT_GC_MAX = 0.70
DEFAULT_POLYT_MIN = 4
DEFAULT_HOMOPOLYMER_MIN = 5
DEFAULT_MIN_ENTROPY = 1.5      # bits, out of a maximum of 2 over 4 bases
DEFAULT_REPEAT_MOTIFS = (2, 3, 4)
DEFAULT_REPEAT_COPIES = 3

STOP_CODONS = ("TAA", "TAG", "TGA")
UNAMBIGUOUS = frozenset("ACGT")


class Run(NamedTuple):
    """A homopolymer run: ``length`` copies of ``base``."""

    base: str
    start: int                  # 1-based inclusive, relative to the region
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start + 1


class Repeat(NamedTuple):
    """A short tandem repeat: ``copies`` consecutive copies of ``motif``."""

    motif: str
    start: int                  # 1-based inclusive, relative to the region
    end: int
    copies: int

    @property
    def length(self) -> int:
        return self.end - self.start + 1


class WindowStat(NamedTuple):
    """One sliding-window measurement."""

    start: int                  # 1-based inclusive, relative to the region
    end: int
    value: float


class Stop(NamedTuple):
    codon: str
    start: int                  # 1-based inclusive, relative to the region

    @property
    def end(self) -> int:
        return self.start + 2


def gc_fraction(sequence: str) -> tuple[float, int]:
    """``(gc_fraction, n_ambiguous)`` over the unambiguous bases only.

    Returns 0.0 for a sequence with no unambiguous base, rather than dividing
    by zero. U is counted with T, so an RNA-alphabet sequence is not reported
    as 100% ambiguous.
    """
    seq = sequence.upper().replace("U", "T")
    counts = {base: seq.count(base) for base in "ACGT"}
    usable = sum(counts.values())
    ambiguous = len(seq) - usable
    if not usable:
        return 0.0, ambiguous
    return (counts["G"] + counts["C"]) / usable, ambiguous


def sliding_gc(sequence: str, window: int = DEFAULT_WINDOW) -> list[WindowStat]:
    """GC fraction of every ``window``-wide window, step 1.

    A sequence shorter than ``window`` yields one window spanning all of it,
    so a short region is measured rather than skipped.
    """
    seq = sequence.upper().replace("U", "T")
    if len(seq) <= window:
        return [WindowStat(1, len(seq), gc_fraction(seq)[0])] if seq else []
    return [
        WindowStat(i + 1, i + window, gc_fraction(seq[i:i + window])[0])
        for i in range(len(seq) - window + 1)
    ]


def homopolymer_runs(
    sequence: str, min_length: int, bases: str = "ACGT"
) -> list[Run]:
    """Every run of one base in ``bases`` at least ``min_length`` long."""
    seq = sequence.upper().replace("U", "T")
    runs = []
    for match in re.finditer(r"(.)\1*", seq):
        base = match.group(1)
        if base in bases and len(match.group(0)) >= min_length:
            runs.append(Run(base, match.start() + 1, match.end()))
    return runs


def tandem_repeats(
    sequence: str,
    motif_sizes: Sequence[int] = DEFAULT_REPEAT_MOTIFS,
    min_copies: int = DEFAULT_REPEAT_COPIES,
) -> list[Repeat]:
    """Maximal short tandem repeats: ``min_copies`` consecutive copies of a motif.

    A motif made of one repeated base is skipped -- ``AA`` x3 is a 6 nt
    poly-A, which `homopolymer_runs` already reports, and reporting it twice
    under two names would double-count one feature.

    ``motif`` is **one rotation** of the repeat, the one beginning at the
    earliest qualifying offset: ``GGCAGCAGCAGGG`` is reported as ``(GCA) x3``
    from offset 2, not as ``(CAG) x3`` from offset 3, because the first is
    where the scan first finds three copies. Both describe the same feature.
    The span is the authoritative part; the motif names it.

    Overlapping reports are suppressed in favour of the longest: a 4-copy
    dinucleotide repeat also contains 3-copy ones, and a ``(AT)x6`` is also
    ``(ATAT)x3``. Scanning longer spans first and skipping anything already
    covered keeps one row per actual feature.
    """
    seq = sequence.upper().replace("U", "T")
    found: list[Repeat] = []
    claimed: list[tuple[int, int]] = []

    candidates: list[Repeat] = []
    for size in sorted(set(motif_sizes)):
        position = 0
        while position + size * min_copies <= len(seq):
            motif = seq[position:position + size]
            if len(set(motif)) == 1:
                position += 1
                continue
            copies = 1
            while seq[position + copies * size:position + (copies + 1) * size] == motif:
                copies += 1
            if copies >= min_copies:
                candidates.append(
                    Repeat(motif, position + 1, position + copies * size, copies)
                )
                # Past the whole repeat: a later start inside it describes the
                # same feature with a rotated motif.
                position += copies * size
            else:
                position += 1

    # Longest first, so the span that survives is the most informative one.
    for repeat in sorted(candidates, key=lambda r: (-r.length, r.start)):
        if any(repeat.start <= end and repeat.end >= start for start, end in claimed):
            continue
        claimed.append((repeat.start, repeat.end))
        found.append(repeat)
    return sorted(found, key=lambda r: r.start)


def shannon_entropy(sequence: str) -> float:
    """Base-composition entropy in bits, 0 for a homopolymer and 2 for uniform."""
    seq = sequence.upper().replace("U", "T")
    counts = [seq.count(base) for base in "ACGT"]
    total = sum(counts)
    if not total:
        return 0.0
    return -sum(
        (n / total) * math.log2(n / total) for n in counts if n
    )


def sliding_entropy(sequence: str, window: int = DEFAULT_WINDOW) -> list[WindowStat]:
    """Shannon entropy of every ``window``-wide window, step 1."""
    seq = sequence.upper().replace("U", "T")
    if len(seq) <= window:
        return [WindowStat(1, len(seq), shannon_entropy(seq))] if seq else []
    return [
        WindowStat(i + 1, i + window, shannon_entropy(seq[i:i + window]))
        for i in range(len(seq) - window + 1)
    ]


def stop_codons(sequence: str) -> list[Stop]:
    """Every TAA/TAG/TGA at every offset, in position order.

    Frame-agnostic by design: this answers "is a stop codon present here",
    which is all the check asks, so a hit at offset 2 counts as much as one at
    offset 1 or 3.

    Two stop codons can never overlap, so no deduplication is needed: none of
    the three contains a T except as its first base, and an overlap would
    require one codon's second or third base to start the next. The scan is
    still written with an overlap-tolerant search, because that is a property
    of these three codons rather than of the scan.
    """
    seq = sequence.upper().replace("U", "T")
    return sorted(
        (
            Stop(codon, index + 1)
            for codon in STOP_CODONS
            for index in _all_offsets(seq, codon)
        ),
        key=lambda s: (s.start, s.codon),
    )


def _all_offsets(sequence: str, needle: str) -> list[int]:
    """Every 0-based offset of ``needle``, overlaps included."""
    offsets = []
    start = sequence.find(needle)
    while start != -1:
        offsets.append(start)
        start = sequence.find(needle, start + 1)
    return offsets


def count_ambiguous(sequence: str) -> int:
    """Bases outside ACGT (U counted as T)."""
    seq = sequence.upper().replace("U", "T")
    return sum(1 for base in seq if base not in UNAMBIGUOUS)


def widest_run(flags: Sequence[bool]) -> int:
    """Longest run of True -- how many consecutive windows were out of band.

    Reported alongside the out-of-band count because the two say different
    things: forty scattered bad windows are a different region from forty
    consecutive ones, and only this separates them.
    """
    widest = current = 0
    for flag in flags:
        current = current + 1 if flag else 0
        widest = max(widest, current)
    return widest
