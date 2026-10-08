"""Classify tiles by which members of a gene family they target.

Given the tiles from stage 1 mapped against a reference of the family's
transcripts, label each tile by the *set* of genes it hits:

``specific``
    One gene only.
``dual``
    Exactly the designed pair, :data:`DUAL_PAIR` (``PIK3CD``/``PIK3CG``).
``pan``
    Every gene in the reference -- conserved across the family.
``other``
    Any other combination: an unintended cross-reaction, not a design goal.
``unmapped``
    No hits at all. Should never happen for a tile taken from the reference
    (see :func:`check_self_hits`), so it means something is wrong upstream.

Scope
-----
The reference is a handful of *paralogous genes*, one RefSeq transcript variant
each, so this measures **within-family cross-reactivity only**. ``specific``
means specific among the sequences in the reference -- not genome-wide, and not
isoform-level within a gene unless the other transcript variants of that gene
are added to the reference. Genome-wide off-target screening is a later stage.

Strand
------
The reference is mRNA sense sequence and a tile is ``target_seq``, verbatim from
the input, so a real binding site is a **forward-strand** hit. A reverse-strand
hit means the tile's reverse complement occurs somewhere, which is not a site a
guide against the sense strand would bind; :func:`hit_sets` drops those. The
caller should still count them, because a nonzero count usually means the input
FASTA holds guides rather than targets.

Validation scope
----------------
Hit sets at ``NM=0`` can be checked exhaustively: :func:`exact_hit_sets` finds
every perfect match by substring search, which for a few transcripts and 30-mers
is both fast and provably complete. It is written without reference to the
aligner so it fails, rather than agreeing, when the aligner is wrong.
:func:`exact_hit_loci` does the same for the *positions* of those matches, so the
``NM=0`` locus columns are validated to the same standard as the hit sets.

**There is no equivalent check above NM=0.** Enumerating every <=n-mismatch match
of a 30-mer is far too large, so mismatch-tolerant hit sets rest on the aligner's
heuristic alone and are advisory: treat a missing hit as "not found", never as
"does not exist".

Separation
----------
Nothing in this module runs an aligner. The classification logic consumes an
iterable of :class:`Hit`, so it is testable without bwa; :func:`hits_from_sam` is
a thin pysam-backed reader that produces them.
"""

from __future__ import annotations

import re
from collections import defaultdict
from typing import Iterable, Iterator, Mapping, NamedTuple

try:
    import pysam
except ImportError as exc:  # pragma: no cover - environment problem, not logic
    raise ImportError(
        "pysam is unavailable. Create and activate the project env:\n"
        "    mamba env create -f environment.yml\n"
        "    conda activate pik3-guide\n"
        f"(underlying error: {exc})"
    ) from exc

#: The one multi-gene combination that is a design goal rather than a surprise.
DUAL_PAIR = frozenset({"PIK3CD", "PIK3CG"})

#: RefSeq headers carry the gene symbol parenthesised, e.g. ``... subunit gamma
#: (PIK3CG), transcript variant 1, mRNA``.
GENE_SYMBOL_PATTERN = re.compile(r"\(([A-Za-z0-9][A-Za-z0-9_.-]*)\)")

#: ``XA`` holds bwa's alternative hits as ``ref,±pos,CIGAR,NM;`` repeated.
_XA_FIELDS = 4

#: CIGAR operations that consume reference bases, for measuring an alignment's span.
_CIGAR_CONSUMES_REFERENCE = frozenset("MDN=X")

_CIGAR_PATTERN = re.compile(r"(\d+)([MIDNSHP=X])")

__all__ = [
    "DUAL_PAIR",
    "GENE_SYMBOL_PATTERN",
    "Hit",
    "check_self_hits",
    "classify",
    "exact_hit_loci",
    "exact_hit_sets",
    "format_loci",
    "gene_symbols",
    "hit_loci",
    "hit_sets",
    "hits_from_sam",
    "reference_span",
    "source_transcript",
]


class Hit(NamedTuple):
    """One alignment of a tile to a reference transcript.

    ``position`` is the 1-based start of the alignment on the reference, matching
    the coordinate convention of ``tile_id`` (see ``src/tile.py``). It defaults to
    ``0`` -- not a valid coordinate -- for hits built in tests of the pure
    set logic, which does not read it; :func:`hits_from_sam` always fills it in.

    ``overhangs`` marks an alignment that runs off the end of its reference, which
    is an artifact of bwa's concatenated index rather than a real match -- see
    :func:`hits_from_sam`. :func:`hit_sets` and :func:`hit_loci` drop these.
    """

    query_id: str
    reference: str
    nm: int
    is_reverse: bool = False
    position: int = 0
    overhangs: bool = False


def source_transcript(tile_id: str) -> str:
    """The transcript a ``tile_id`` came from, e.g. ``NM_006218.4:1-30`` -> ``NM_006218.4``.

    Splits on the last ``:`` so accessions containing one are safe.
    """
    transcript_id, _, coordinates = tile_id.rpartition(":")
    if not transcript_id or not coordinates:
        raise ValueError(f"not a tile_id of the form transcript:start-end: {tile_id!r}")
    return transcript_id


def gene_symbols(records: Iterable) -> dict[str, str]:
    """Map ``transcript_id`` -> gene symbol, parsed from each record's description.

    Takes :class:`tile.Record`-alikes (anything with ``transcript_id`` and
    ``description``). Raises ``ValueError`` unless exactly one parenthesised
    symbol is present: guessing, or falling back to the accession, would quietly
    corrupt every label downstream, so this fails loudly instead.
    """
    symbols: dict[str, str] = {}
    for record in records:
        found = GENE_SYMBOL_PATTERN.findall(record.description or "")
        if len(found) != 1:
            raise ValueError(
                f"{record.transcript_id}: expected exactly one parenthesised gene "
                f"symbol in the FASTA description, found {len(found)} "
                f"({found or 'none'}). Description: {record.description!r}"
            )
        symbols[record.transcript_id] = found[0]
    return symbols


def reference_span(cigar: str) -> int:
    """How many reference bases a CIGAR string covers, for bounds-checking a hit.

    Only ``M/D/N/=/X`` consume reference; insertions, clips and padding do not.
    Used on ``XA`` entries, whose CIGAR is a string rather than a parsed alignment.
    """
    operations = _CIGAR_PATTERN.findall(cigar)
    if not operations or "".join(f"{count}{op}" for count, op in operations) != cigar:
        raise ValueError(f"malformed CIGAR string: {cigar!r}")
    return sum(
        int(count) for count, op in operations if op in _CIGAR_CONSUMES_REFERENCE
    )


def hits_from_sam(path: str) -> Iterator[Hit]:
    """Read every hit from a bwa ``samse`` SAM file, primary alignments and ``XA`` alike.

    bwa reports one record per read; the additional places it mapped are packed
    into the ``XA`` tag, so the complete hit set is the primary alignment plus
    every ``XA`` entry. Each ``XA`` entry carries its **own** ``NM``, which is
    parsed from the tag rather than inherited from the primary -- an alternative
    hit is usually a worse match than the primary, and treating them alike would
    silently promote mismatched hits into the exact set.

    Records are gated on whether an alignment was *reported* (``RNAME`` set), not
    on the unmapped flag. bwa sets flag ``0x4`` on reads with many equally-best
    hits while still filling in ``RNAME``, ``POS``, ``NM`` and a full ``XA`` list
    -- low-complexity tiles (a poly-A tail, a tandem repeat) hit this routinely.
    Trusting ``is_unmapped`` there discards perfect alignments, including a tile's
    own source transcript, which trips :func:`check_self_hits`. ``XA`` is parsed
    whatever the flag says, for the same reason. The ``NM=0`` oracle is what makes
    accepting these records safe.

    Hits whose alignment runs past the end of the reference are flagged
    ``overhangs`` rather than dropped here, so the dump can show them. bwa indexes
    the references as one concatenated sequence, and an ``XA`` entry can report an
    alignment that crosses the junction into the *next* transcript while still
    being attributed to the previous one, with a full-length CIGAR and ``NM:i:0``::

        NM_000002.1:32-61 ... XA:Z:NM_000002.1,+212,30M,0;   <- reference is 240 nt

    A 30M at 212 ends at 241, so one base came from another transcript: within this
    reference the match does not exist, and the substring oracle rightly disagrees.
    Such a hit can inject a reference the tile never matched into its gene set, and
    it is invisible to a set-level check whenever that reference is present anyway.
    :func:`hit_sets` and :func:`hit_loci` therefore skip them.

    Reads with no alignment at all yield nothing. Deduplication is
    :func:`hit_sets`' job.
    """
    with pysam.AlignmentFile(path, "r") as samfile:
        lengths = dict(zip(samfile.references, samfile.lengths))

        def end_of(reference: str) -> int:
            """The reference's length, or an error -- never a default.

            A missing length would make every hit on that reference look like it
            ran off the end, silently emptying its hit sets.
            """
            if reference not in lengths:
                raise ValueError(
                    f"{path}: reference {reference!r} is not in the SAM header, so "
                    "its length is unknown and hits on it cannot be bounds-checked"
                )
            return lengths[reference]

        for read in samfile:
            if read.reference_name is not None:
                yield Hit(
                    query_id=read.query_name,
                    reference=read.reference_name,
                    nm=int(read.get_tag("NM")) if read.has_tag("NM") else 0,
                    is_reverse=bool(read.is_reverse),
                    # pysam's reference_start is 0-based; tile coordinates are
                    # 1-based inclusive, and XA positions are already 1-based.
                    position=read.reference_start + 1,
                    overhangs=(read.reference_end or 0) > end_of(read.reference_name),
                )
            if not read.has_tag("XA"):
                continue
            for entry in str(read.get_tag("XA")).split(";"):
                if not entry:
                    continue
                fields = entry.split(",")
                if len(fields) != _XA_FIELDS:
                    raise ValueError(
                        f"{read.query_name}: malformed XA entry {entry!r} "
                        f"(expected {_XA_FIELDS} comma-separated fields)"
                    )
                reference, position, cigar, nm = fields
                start = abs(int(position))
                yield Hit(
                    query_id=read.query_name,
                    reference=reference,
                    nm=int(nm),
                    is_reverse=position.startswith("-"),
                    position=start,
                    overhangs=start + reference_span(cigar) - 1 > end_of(reference),
                )


def hit_sets(
    hits: Iterable[Hit], max_nm: int = 0
) -> tuple[dict[str, frozenset[str]], dict[str, int], int]:
    """Collapse hits into one reference set per tile, keeping only good enough ones.

    Returns ``(sets, worst_nm, reverse_dropped)``: the forward-strand references
    each tile hit with ``nm <= max_nm``, the largest such ``nm`` per tile, and how
    many reverse-strand hits were discarded.

    Reverse-strand hits are dropped because the reference is mRNA sense sequence
    (see the module docstring), and hits running off the end of their reference
    because they are index artifacts rather than matches (see
    :func:`hits_from_sam`); only the former are counted, being the diagnostic one. A large ``reverse_dropped`` usually means guides
    were mapped where targets were intended.
    """
    if max_nm < 0:
        raise ValueError(f"max_nm must be >= 0, got {max_nm}")

    sets: dict[str, set[str]] = defaultdict(set)
    worst_nm: dict[str, int] = {}
    reverse_dropped = 0
    for hit in hits:
        if hit.is_reverse:
            reverse_dropped += 1
            continue
        if hit.overhangs or hit.nm > max_nm:
            continue
        sets[hit.query_id].add(hit.reference)
        worst_nm[hit.query_id] = max(worst_nm.get(hit.query_id, 0), hit.nm)
    return {query: frozenset(refs) for query, refs in sets.items()}, worst_nm, reverse_dropped


def hit_loci(
    hits: Iterable[Hit], max_nm: int = 0
) -> dict[str, dict[str, list[tuple[int, int]]]]:
    """Where each tile hit, per reference: ``{tile_id: {reference: [(position, nm)]}}``.

    The companion to :func:`hit_sets`, which deduplicates by reference name and so
    collapses a tile matching the *same* reference at several positions into one
    entry. This keeps every position, which is what distinguishes a self-repeat
    (two loci on one reference) from genuine cross-family reactivity (loci on
    different references) -- a distinction ``max_nm`` alone conflates.

    Applies the same filters as :func:`hit_sets` -- reverse-strand and
    reference-overhanging hits dropped, ``nm <= max_nm`` kept -- so the references appearing here are exactly those in
    the corresponding hit set. Identical ``(position, nm)`` pairs are deduplicated,
    because bwa can report the primary alignment again as an ``XA`` entry. Loci are
    sorted by position.
    """
    if max_nm < 0:
        raise ValueError(f"max_nm must be >= 0, got {max_nm}")

    loci: dict[str, dict[str, set[tuple[int, int]]]] = defaultdict(lambda: defaultdict(set))
    for hit in hits:
        if hit.is_reverse or hit.overhangs or hit.nm > max_nm:
            continue
        loci[hit.query_id][hit.reference].add((hit.position, hit.nm))
    return {
        query: {reference: sorted(pairs) for reference, pairs in sorted(by_reference.items())}
        for query, by_reference in loci.items()
    }


def format_loci(by_reference: Mapping[str, Iterable[tuple[int, int]]], length: int) -> str:
    """Render one tile's loci as ``accession:start-end(nm)``, for a TSV cell.

    ``;`` separates references, ``,`` separates several loci on the same reference.
    Coordinates are 1-based inclusive, the same convention as ``tile_id``, so a
    locus string can be read as a coordinate on the reference directly::

        NM_006218.4:2724-2753(0);NM_006219.3:2765-2794(3)

    ``length`` is the tile's length, used to derive the end coordinate; alignments
    are ungapped (``bwa aln -o 0``), so the end is always ``start + length - 1``.
    """
    return ";".join(
        ",".join(
            f"{reference}:{position}-{position + length - 1}({nm})" for position, nm in loci
        )
        for reference, loci in sorted(by_reference.items())
    )


def classify(hit_genes: frozenset[str], all_genes: frozenset[str]) -> str:
    """Label a tile from the set of genes it hits. See the module docstring.

    ``pan`` is tested before ``dual`` so that a reference containing *only*
    ``PIK3CD`` and ``PIK3CG`` labels a both-genes tile ``pan`` -- "conserved
    across the whole reference" is the stronger and more honest claim. With the
    full four-gene family the two cannot collide.
    """
    unknown = hit_genes - all_genes
    if unknown:
        raise ValueError(f"hit genes absent from the reference: {sorted(unknown)}")
    if not hit_genes:
        return "unmapped"
    if hit_genes == all_genes:
        return "pan"
    if len(hit_genes) == 1:
        return "specific"
    if hit_genes == DUAL_PAIR:
        return "dual"
    return "other"


def exact_hit_sets(
    tiles: Mapping[str, str], references: Mapping[str, str]
) -> dict[str, frozenset[str]]:
    """Every perfect forward-strand match, by substring search.

    The oracle the ``NM=0`` columns are validated against: complete by
    construction, and written without reference to the aligner so it fails rather
    than agrees when the aligner is wrong. Both arguments map name -> sequence.

    Only sound at ``NM=0``; there is deliberately no mismatch-tolerant
    counterpart (see the module docstring).
    """
    upper_references = {name: sequence.upper() for name, sequence in references.items()}
    return {
        tile_id: frozenset(
            name for name, sequence in upper_references.items() if tile.upper() in sequence
        )
        for tile_id, tile in tiles.items()
    }


def exact_hit_loci(
    tiles: Mapping[str, str], references: Mapping[str, str]
) -> dict[str, dict[str, list[int]]]:
    """Every perfect forward-strand match **with its position**, by substring search.

    The oracle for the ``NM=0`` locus columns, the position-aware counterpart of
    :func:`exact_hit_sets`: ``{tile_id: {reference: [1-based positions]}}``, with
    overlapping occurrences included, so a tandem repeat yields every copy. As
    complete by construction as the hit-set oracle, and likewise written without
    reference to the aligner.

    Only sound at ``NM=0``; enumerating mismatched positions is exactly the
    intractable problem described in the module docstring.
    """
    upper_references = {name: sequence.upper() for name, sequence in references.items()}
    found: dict[str, dict[str, list[int]]] = {}
    for tile_id, tile in tiles.items():
        sequence = tile.upper()
        per_reference: dict[str, list[int]] = {}
        for name, reference in upper_references.items():
            positions = []
            index = reference.find(sequence)
            while index != -1:
                positions.append(index + 1)
                index = reference.find(sequence, index + 1)
            if positions:
                per_reference[name] = positions
        found[tile_id] = per_reference
    return found


def check_self_hits(hit_sets_exact: Mapping[str, frozenset[str]]) -> list[str]:
    """Tile ids that failed to hit their own source transcript at ``NM=0``.

    Every tile was cut out of a reference transcript, so it must match that
    transcript perfectly. Anything listed here is an alignment the mapper lost,
    which invalidates the hit sets rather than merely trimming them -- the caller
    should fail rather than publish confident labels over incomplete data.
    """
    return sorted(
        tile_id
        for tile_id, references in hit_sets_exact.items()
        if source_transcript(tile_id) not in references
    )
