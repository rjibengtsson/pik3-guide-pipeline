"""Read exon and CDS annotation off a RefSeq mRNA GenBank record.

Coordinates follow src/tile.py: 1-based inclusive, so a feature spans
``end - start + 1`` bases and ``sequence[start - 1:end]`` is its sequence.
Biopython reports a feature location as 0-based with an exclusive end, so
``start`` is shifted by one here and ``end`` is taken as-is.

**On an mRNA record the exons are contiguous.** The introns have been spliced
out, so the features tile the transcript from base 1 with no gaps --
NM_006218.4 runs 1..247, 248..675, 676..885 and so on. The useful signal is
therefore not "exon versus intron", which does not exist on this coordinate
system, but where one exon *ends* and the next begins. `exons_from_genbank`
returns them in transcript order with a 1-based `index` so a caller can make
those boundaries visible.

That contiguity is asserted, not assumed: a record whose exons overlap, leave
a gap or arrive out of order is a genomic record or a mis-parse rather than an
mRNA, and anything drawn from it would be quietly misplaced.

The exon table may still stop short of the sequence it annotates --
NM_002649.3 is 7218 nt but its last exon ends at 7210 -- so `genbank_length`
is offered separately from the exon span, and the two disagreeing is a gap in
NCBI's annotation rather than an error here.

This is library code: it parses no command line and prints nothing.
"""

from __future__ import annotations

import glob
import os
from typing import NamedTuple

from Bio import SeqIO


class Exon(NamedTuple):
    index: int                  # 1-based position in transcript order
    start: int                  # 1-based inclusive
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start + 1


class Cds(NamedTuple):
    start: int                  # 1-based inclusive
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start + 1


def _only_record(path: str):
    """One record per file, matching the one-transcript-per-file convention."""
    records = list(SeqIO.parse(path, "genbank"))
    if len(records) != 1:
        raise ValueError(
            f"{path}: expected exactly one GenBank record, found "
            f"{len(records)}. One transcript per file is the convention "
            "across this pipeline."
        )
    return records[0]


def find_genbank(accession: str, directories: list[str]) -> str | None:
    """Locate ``{accession}.gbk``, trying the directories in order.

    Falls back to a version-insensitive match, so ``NM_002649`` finds
    ``NM_002649.3.gbk``: the accessions get written down without their version
    suffix routinely, and silently drawing no annotation because of a missing
    ``.3`` is worse than matching on the stem. An exact match always wins, and
    an ambiguous stem raises rather than picking one.
    """
    for directory in directories:
        exact = os.path.join(directory, f"{accession}.gbk")
        if os.path.isfile(exact):
            return exact

    stem = accession.split(".")[0]
    matches = sorted(
        {
            path
            for directory in directories
            for path in glob.glob(os.path.join(directory, f"{stem}.*.gbk"))
        }
    )
    if len(matches) > 1:
        raise ValueError(
            f"{accession}: {len(matches)} GenBank files match the unversioned "
            f"accession ({', '.join(os.path.basename(m) for m in matches)}). "
            "Name the version explicitly."
        )
    return matches[0] if matches else None


def exons_from_genbank(path: str) -> list[Exon]:
    """Every ``exon`` feature, in transcript order, 1-based inclusive.

    Raises if the record holds no exon features: a FASTA renamed to ``.gbk``,
    or a record stripped of its feature table, would otherwise yield an empty
    annotation that reads as "this transcript has one exon".
    """
    record = _only_record(path)

    exons = [
        (int(f.location.start) + 1, int(f.location.end))
        for f in record.features
        if f.type == "exon"
    ]
    if not exons:
        raise ValueError(
            f"{path}: no 'exon' features in the record. A RefSeq mRNA .gbk "
            "carries one per exon; a FASTA renamed to .gbk carries none."
        )

    # Sorted rather than trusted in file order: NM_006218.4 lists exons 1 and 2,
    # then the CDS and its misc_features, then exons 3 onward, so file order is
    # not transcript order in general even though it happens to be here.
    exons.sort()

    for (prev_start, prev_end), (start, _) in zip(exons, exons[1:]):
        if start != prev_end + 1:
            gap = start - prev_end - 1
            raise ValueError(
                f"{path}: exons {prev_start}..{prev_end} and {start}.. are "
                f"{'separated by %d nt' % gap if gap > 0 else 'overlapping'}, "
                "so they do not tile the transcript. This module expects a "
                "spliced mRNA record, not a genomic one."
            )

    return [Exon(i, start, end) for i, (start, end) in enumerate(exons, start=1)]


def cds_from_genbank(path: str) -> Cds | None:
    """The single ``CDS`` feature, or None if the record has none."""
    record = _only_record(path)

    cds = [f for f in record.features if f.type == "CDS"]
    if not cds:
        return None
    if len(cds) > 1:
        raise ValueError(
            f"{path}: {len(cds)} CDS features. One coding sequence per mRNA "
            "is assumed; pick the intended one explicitly."
        )

    location = cds[0].location
    return Cds(int(location.start) + 1, int(location.end))


def gene_symbol_from_genbank(path: str) -> str | None:
    """The ``/gene`` qualifier, or None if the record carries none.

    Taken from the ``gene`` feature, falling back to the ``CDS``: the two agree
    on every record in this repo, but a record can carry one without the other.
    Read from the qualifier rather than parsed out of the DESCRIPTION line,
    which is prose and varies between records.

    Returns None rather than guessing. A caller naming a file after this must
    cope with its absence, since a symbol invented from the accession would
    label the file with something that is not the gene.
    """
    record = _only_record(path)

    for feature_type in ("gene", "CDS"):
        for feature in record.features:
            if feature.type != feature_type:
                continue
            symbols = feature.qualifiers.get("gene") or []
            if symbols:
                return str(symbols[0])
    return None


def genbank_length(path: str) -> int:
    """Length of the record's own sequence, i.e. its LOCUS ``bp``.

    Kept separate from the exon span because the two can disagree *within one
    record*: NM_002649.3 is 7218 bp but its last exon ends at 7210, so eight
    3' bases carry no exon feature. That is a gap in NCBI's annotation, not a
    coordinate error -- telling it apart from a genuine version mismatch
    against another file needs this number.
    """
    return len(_only_record(path).seq)
