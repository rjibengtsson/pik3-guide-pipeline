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

__all__ = [
    "DUAL_PAIR",
    "GENE_SYMBOL_PATTERN",
    "Hit",
    "check_self_hits",
    "classify",
    "exact_hit_sets",
    "gene_symbols",
    "hit_sets",
    "hits_from_sam",
    "source_transcript",
]


class Hit(NamedTuple):
    """One alignment of a tile to a reference transcript."""

    query_id: str
    reference: str
    nm: int
    is_reverse: bool = False


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

    Reads with no alignment at all yield nothing. Deduplication is
    :func:`hit_sets`' job.
    """
    with pysam.AlignmentFile(path, "r") as samfile:
        for read in samfile:
            if read.reference_name is not None:
                yield Hit(
                    query_id=read.query_name,
                    reference=read.reference_name,
                    nm=int(read.get_tag("NM")) if read.has_tag("NM") else 0,
                    is_reverse=bool(read.is_reverse),
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
                reference, position, _cigar, nm = fields
                yield Hit(
                    query_id=read.query_name,
                    reference=reference,
                    nm=int(nm),
                    is_reverse=position.startswith("-"),
                )


def hit_sets(
    hits: Iterable[Hit], max_nm: int = 0
) -> tuple[dict[str, frozenset[str]], dict[str, int], int]:
    """Collapse hits into one reference set per tile, keeping only good enough ones.

    Returns ``(sets, worst_nm, reverse_dropped)``: the forward-strand references
    each tile hit with ``nm <= max_nm``, the largest such ``nm`` per tile, and how
    many reverse-strand hits were discarded.

    Reverse-strand hits are dropped because the reference is mRNA sense sequence
    (see the module docstring). A large ``reverse_dropped`` usually means guides
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
        if hit.nm > max_nm:
            continue
        sets[hit.query_id].add(hit.reference)
        worst_nm[hit.query_id] = max(worst_nm.get(hit.query_id, 0), hit.nm)
    return {query: frozenset(refs) for query, refs in sets.items()}, worst_nm, reverse_dropped


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
