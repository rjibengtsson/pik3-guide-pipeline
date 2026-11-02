"""Sequence QC at two scales: the cloned insert, and the individual guide.

The module serves **two units of measurement**, and keeping them apart is the
whole organising idea here. Which one a function belongs to is stated in its
docstring; the guide-scale functions are gathered at the bottom of the file
under their own heading, and take a single guide as their argument.

Insert scale -- GC, poly-T, repeats, stop codons
    Measured over a region already extracted from its transcript. The unit is
    the **whole insert**, because that is what gets cloned: these numbers
    describe the fragment going into the construct, not the individual spacers
    that will later be chosen from inside it.
Guide scale -- GC band, homopolymer runs
    Measured over one candidate guide, asked of the guide itself. This is
    deliberately *not* a sliding window over an insert: a window is a guess at
    where a guide might go, whereas a guide is a sequence someone has committed
    to, and only the second is worth passing or failing. `check_guide` is the
    entry point.

The two scales do not share thresholds. A run of 4 identical bases is
unremarkable in 2 kb of insert and disqualifying in a 30 nt spacer, and the GC
band that an insert must merely sit inside is much wider than the band a guide
is held to -- hence `DEFAULT_GC_MIN`/`DEFAULT_GC_MAX` for the insert and the
separate `DEFAULT_GUIDE_GC_*` for a guide. They do share *scanners*:
`check_guide` calls the same `gc_fraction` and `homopolymer_runs` the insert
checks use, so the two scales cannot disagree about what a run or a GC fraction
is, only about which values are acceptable.

The insert-scale properties, in detail:

GC content
    One fraction over the whole region -- never a scan at spacer scale. A
    guide's own GC is a property of that guide and is judged by `check_guide`
    against its own band, not read off a window of the insert.
poly-T stretches
    Runs of T at or above a threshold. Four is the default because a run of
    four or more T is the classic premature-termination signal in a U6-driven
    construct.
repeats and homopolymers
    Homopolymers of A, C and G (T is covered by the poly-T scan, so it is not
    reported twice), short tandem repeats of a 2-4 nt motif, and the region's
    Shannon entropy, which catches a composition so skewed that a motif scan
    finds nothing to point at.
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

DEFAULT_GC_MIN = 0.30
DEFAULT_GC_MAX = 0.70
DEFAULT_POLYT_MIN = 4
DEFAULT_HOMOPOLYMER_MIN = 5
DEFAULT_MIN_ENTROPY = 1.5      # bits, out of a maximum of 2 over 4 bases
DEFAULT_REPEAT_MOTIFS = (2, 3, 4)
DEFAULT_REPEAT_COPIES = 3

# Guide scale. Separate from the insert thresholds above on purpose: see the
# module docstring. A 30 nt spacer is held to a tighter GC band and a shorter
# homopolymer than 2 kb of insert, and all four bases count for a guide --
# there is no poly-T scan to hand T off to, because a guide's U run is a
# property of the spacer, not a U6 termination signal in the construct.
DEFAULT_GUIDE_GC_MIN = 0.40
DEFAULT_GUIDE_GC_MAX = 0.60
DEFAULT_GUIDE_HOMOPOLYMER_MIN = 4

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
    """Base-composition entropy in bits, 0 for a homopolymer and 2 for uniform.

    Over a whole insert this is a coarse measure by construction -- a short
    low-complexity stretch barely moves it. That stretch is what
    `homopolymer_runs` and `tandem_repeats` are for; this number only catches a
    composition skewed across the whole fragment.
    """
    seq = sequence.upper().replace("U", "T")
    counts = [seq.count(base) for base in "ACGT"]
    total = sum(counts)
    if not total:
        return 0.0
    return -sum(
        (n / total) * math.log2(n / total) for n in counts if n
    )


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


# ---------------------------------------------------------------------------
# Guide scale
#
# Everything below takes one candidate guide as its unit. See the module
# docstring for why that is kept distinct from the insert-scale functions
# above, and why the thresholds differ.
# ---------------------------------------------------------------------------


class GuideQC(NamedTuple):
    """One guide's QC record, carrying the thresholds that judged it.

    The thresholds travel with the measurements so that `flags` is derivable
    from the record alone -- the same reason the insert summary echoes its
    thresholds into every row. A flag can then be audited without knowing what
    the CLI was invoked with.
    """

    guide_id: str
    sequence: str
    gc: float
    n_ambiguous: int
    runs: tuple[Run, ...]
    gc_min: float
    gc_max: float
    homopolymer_min: int

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def gc_low(self) -> bool:
        return self.gc < self.gc_min

    @property
    def gc_high(self) -> bool:
        return self.gc > self.gc_max

    @property
    def gc_in_band(self) -> bool:
        """Band edges are **inclusive**: exactly 40% GC passes a 0.40 minimum.

        Matching the insert-scale check, which also compares with ``<=``. The
        edge is reachable in practice -- 12 G/C in a 30 nt guide is exactly
        0.40 -- so which side of it passes is a real decision, not a rounding
        detail.
        """
        return self.gc_min <= self.gc <= self.gc_max

    @property
    def longest_run(self) -> Run | None:
        """The longest homopolymer run, earliest one first on a tie."""
        return min(
            self.runs, key=lambda r: (-r.length, r.start), default=None
        )

    @property
    def flags(self) -> tuple[str, ...]:
        """Why this guide is flagged; empty for a guide that passes.

        ``gc_low`` and ``gc_high`` are separate rather than one
        ``gc_out_of_band``, because the two have opposite remedies and a count
        of each is what says whether a band is set wrong or a transcript is
        simply GC-skewed.
        """
        flags = []
        if self.gc_low:
            flags.append("gc_low")
        if self.gc_high:
            flags.append("gc_high")
        if self.runs:
            flags.append("homopolymer")
        return tuple(flags)

    @property
    def passes(self) -> bool:
        return not self.flags


def guide_homopolymer_runs(
    sequence: str, min_length: int = DEFAULT_GUIDE_HOMOPOLYMER_MIN
) -> list[Run]:
    """Runs of >= ``min_length`` identical bases, named in the input alphabet.

    A thin wrapper over `homopolymer_runs` over **all four** bases -- there is
    no poly-T scan at guide scale to report T separately, so T is not excluded
    the way it is for an insert.

    The wrapper exists for the alphabet. `homopolymer_runs` normalises U to T
    before scanning, which is right for finding the run but wrong for naming
    it: a guide RNA's run of four U would come back as ``base='T'``, and
    ``'T' * 4`` is not a sequence that appears anywhere in the guide. Here the
    reported base is read back out of the input, so a UUUU run reports ``U``
    and a GGGG run reports ``G``.

    A run that is only identical *after* normalisation -- ``UUTT``, from a
    sequence mixing the two alphabets -- is reported under the normalised base
    ``T``, since no single input character describes it. Such sequences are a
    known rough edge of this pipeline rather than a supported input; the span
    is still correct, and the span is the actionable part.
    """
    original = sequence.upper()
    runs = []
    for run in homopolymer_runs(original, min_length, bases="ACGT"):
        spelt = set(original[run.start - 1:run.end])
        base = spelt.pop() if len(spelt) == 1 else run.base
        runs.append(Run(base, run.start, run.end))
    return runs


def check_guide(
    guide_id: str,
    sequence: str,
    gc_min: float = DEFAULT_GUIDE_GC_MIN,
    gc_max: float = DEFAULT_GUIDE_GC_MAX,
    homopolymer_min: int = DEFAULT_GUIDE_HOMOPOLYMER_MIN,
) -> GuideQC:
    """QC one guide: GC against a band, and homopolymer runs against a length.

    Accepts either alphabet, since `gc_fraction` counts U with T. The sequence
    is kept verbatim on the record apart from upper-casing, following the
    ``target_seq`` convention in src/tile.py.

    Bases outside ACGTU are excluded from the GC denominator and counted in
    ``n_ambiguous`` rather than scored as not-GC, so an N-containing guide's GC
    is the GC of its real bases. They are reported, not flagged: the count is
    on the record for a caller that wants to reject on it, but an ambiguous
    base is a property of the input FASTA rather than of the guide's chemistry.
    """
    gc, ambiguous = gc_fraction(sequence)
    return GuideQC(
        guide_id=guide_id,
        sequence=sequence.upper(),
        gc=gc,
        n_ambiguous=ambiguous,
        runs=tuple(guide_homopolymer_runs(sequence, homopolymer_min)),
        gc_min=gc_min,
        gc_max=gc_max,
        homopolymer_min=homopolymer_min,
    )
