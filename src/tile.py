"""Tiling k-mer extraction from mRNA transcript sequences.

Coordinate convention
---------------------
All positions are **1-based, inclusive** on the forward strand of the
transcript as it appears in the input FASTA. The first tile of a sequence is
therefore ``start=1, end=k``, and the relationship to Python slicing is::

    sequence[start - 1:end] == tile

Direction
---------
Every sequence in and out of this module is written **5' -> 3'**, the
conventional direction for reporting and ordering an oligo. Input FASTA
sequences are 5' -> 3', tiles inherit that, and :func:`guide_rna` returns the
spacer 5' -> 3' as well -- so the guide's 5' end pairs with the target's 3' end
(they are antiparallel). The test suite's ``pairs_antiparallel`` helper asserts
that relationship independently of how this module computes it.

mRNA is single-stranded, so tiles are only ever taken from the ``+`` strand.

Alphabets
---------
Nothing here infers whether a sequence is DNA or RNA. Tiles are reported
verbatim from the input, in whatever alphabet the FASTA used; the caller
declares the molecule type. Only :func:`guide_rna` converts, because a guide has
to base-pair as RNA.

Biopython
---------
Sequence manipulation is Biopython's: ``SeqIO`` parses FASTA, and ``transcribe``
(T -> U), ``back_transcribe`` (U -> T), ``reverse_complement`` (DNA alphabet) and
``reverse_complement_rna`` do the alphabet and strand work; they are re-exported
here so callers need one import for the whole stage. Local to this module are
the tiling itself, the provenance tracking and the guide convention.

This requires Biopython *and* numpy -- ``Bio/SeqIO/__init__.py`` imports TwoBitIO
unconditionally, so a missing numpy breaks plain FASTA parsing too. Use the
``pik3-guide`` env from ``environment.yml``; the base conda env has no numpy and
cannot import ``Bio.SeqIO`` at all.
"""

from __future__ import annotations

import os
from typing import Iterable, Iterator, NamedTuple

try:
    from Bio import SeqIO
    from Bio.Data import IUPACData
    from Bio.Seq import (
        back_transcribe,
        reverse_complement,
        reverse_complement_rna,
        transcribe,
    )
except ImportError as exc:  # pragma: no cover - environment problem, not logic
    # Biopython's own message here talks about TwoBit files, which is misleading:
    # Bio/SeqIO/__init__.py imports TwoBitIO unconditionally, so a missing numpy
    # breaks plain FASTA parsing too.
    raise ImportError(
        "Bio.SeqIO is unavailable. It needs Biopython *and* numpy. Create and "
        "activate the project env:\n"
        "    mamba env create -f environment.yml\n"
        "    conda activate pik3-guide\n"
        f"(underlying error: {exc})"
    ) from exc

FASTA_EXTENSIONS = (".fasta", ".fa", ".fna", ".txt")

#: Unambiguous bases, DNA and RNA ({A, C, G, T, U}). Anything else in a tile --
#: N or any other IUPAC ambiguity code -- sets ``has_ambiguous``.
UNAMBIGUOUS_BASES = frozenset(
    IUPACData.unambiguous_dna_letters + IUPACData.unambiguous_rna_letters
)

#: ``transcribe`` (T -> U), ``back_transcribe`` (U -> T) and
#: ``reverse_complement`` (DNA alphabet, A pairs with T) are re-exported from
#: Biopython so callers need one import for the whole stage.
__all__ = [
    "FASTA_EXTENSIONS",
    "UNAMBIGUOUS_BASES",
    "Record",
    "Tile",
    "back_transcribe",
    "find_fasta_files",
    "guide_rna",
    "parse_fasta",
    "read_fasta",
    "reverse_complement",
    "tile_record",
    "tile_sequence",
    "transcribe",
]


class Record(NamedTuple):
    """One FASTA record."""

    transcript_id: str
    description: str
    sequence: str
    source_file: str


class Tile(NamedTuple):
    """One tiled k-mer, with its provenance."""

    transcript_id: str
    source_file: str
    start: int
    end: int
    sequence: str

    @property
    def tile_id(self) -> str:
        """Stable join key for downstream results, e.g. ``NM_006218.4:1-30``."""
        return f"{self.transcript_id}:{self.start}-{self.end}"

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def has_ambiguous(self) -> bool:
        return not set(self.sequence.upper()) <= UNAMBIGUOUS_BASES


def guide_rna(target: str) -> str:
    """The Cas13 spacer for an mRNA target, written 5' -> 3'.

    ``target`` is a tile written 5' -> 3'; the returned guide is its reverse
    complement in the RNA alphabet, also written **5' -> 3'**, which is the
    conventional direction for ordering and reporting an oligo.

    Guide and target are antiparallel, so the guide's 5' end pairs with the
    target's 3' end::

        target  5'-AAAAACCCCCGGGGGUUUUUAAGGCCUUAC-3'
        guide   3'-UUUUUGGGGGCCCCCAAAAUUCCGGAAUG-5'   (same molecule, flipped)
        guide   5'-GUAAGGCCUUAAAAACCCCCGGGGGUUUUU-3'   <- what is returned

    ``A`` in the target therefore becomes ``U`` in the guide, not ``T`` as
    ``reverse_complement`` alone would give. A thin wrapper over Biopython's
    ``reverse_complement_rna``: it exists to name the biology (this is a guide,
    not just a transformed string) and to carry the direction convention above.

    ``target`` may be in either alphabet -- a DNA-alphabet target, as RefSeq
    writes mRNA, gives the same guide as its RNA equivalent.
    """
    return reverse_complement_rna(target)


def parse_fasta(handle: Iterable[str], source_file: str = "") -> Iterator[Record]:
    """Parse FASTA from a handle (or any iterable of lines) via ``Bio.SeqIO``.

    Sequences are uppercased; Biopython handles line wrapping, blank lines,
    CRLF and internal whitespace. The transcript id is ``SeqRecord.id``, the
    first whitespace-delimited token of the header, which for RefSeq FASTA is
    the accession (``NM_006218.4``). ``description`` is the remainder of the
    header, with the id stripped off the front -- Biopython repeats the id in
    ``SeqRecord.description``, which this deliberately does not.

    Raises ``ValueError`` on malformed input, including a record whose header is
    just ``>``: Biopython returns an empty id there, which would produce useless
    ``tile_id`` values like ``:1-30``.
    """
    where = source_file or "<input>"
    for seq_record in SeqIO.parse(handle, "fasta"):
        transcript_id = seq_record.id or ""
        if not transcript_id:
            raise ValueError(f"{where}: FASTA record with empty header")
        description = seq_record.description or ""
        if description.startswith(transcript_id):
            description = description[len(transcript_id) :].strip()
        yield Record(
            transcript_id=transcript_id,
            description=description,
            sequence=str(seq_record.seq).upper(),
            source_file=source_file,
        )


def read_fasta(path: str) -> Iterator[Record]:
    """Parse every record in the FASTA file at ``path``."""
    with open(path, "r", encoding="utf-8") as handle:
        yield from parse_fasta(handle, source_file=os.path.basename(path))


def find_fasta_files(input_dir: str) -> list[str]:
    """Return FASTA paths in ``input_dir``, sorted for deterministic output."""
    if not os.path.isdir(input_dir):
        raise NotADirectoryError(f"not a directory: {input_dir}")
    return sorted(
        os.path.join(input_dir, name)
        for name in os.listdir(input_dir)
        if name.lower().endswith(FASTA_EXTENSIONS)
        and os.path.isfile(os.path.join(input_dir, name))
    )


def tile_sequence(
    sequence: str,
    k: int = 30,
    step: int = 1,
    transcript_id: str = "",
    source_file: str = "",
) -> Iterator[Tile]:
    """Yield every ``k``-mer of ``sequence``, advancing ``step`` bases at a time.

    With ``step=1`` the number of tiles is ``len(sequence) - k + 1``. A sequence
    shorter than ``k`` yields nothing.
    """
    if k < 1:
        raise ValueError(f"k must be >= 1, got {k}")
    if step < 1:
        raise ValueError(f"step must be >= 1, got {step}")

    for start0 in range(0, len(sequence) - k + 1, step):
        yield Tile(
            transcript_id=transcript_id,
            source_file=source_file,
            start=start0 + 1,
            end=start0 + k,
            sequence=sequence[start0 : start0 + k],
        )


def tile_record(record: Record, k: int = 30, step: int = 1) -> Iterator[Tile]:
    """Tile one :class:`Record`, carrying its provenance onto each tile."""
    yield from tile_sequence(
        record.sequence,
        k=k,
        step=step,
        transcript_id=record.transcript_id,
        source_file=record.source_file,
    )
