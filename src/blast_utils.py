"""Reusable blastn helpers
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, asdict
from typing import Iterable, Sequence

# Explicit -outfmt 6 field list, in order. Keep in sync with BlastHit.from_row.
OUTFMT_FIELDS = (
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qlen",
    "qcovhsp",
)

# For long sequences use megablast.
DEFAULT_TASK = "megablast"
DEFAULT_WORD_SIZE = 28
DEFAULT_EVALUE = 1e-5
DEFAULT_MAX_TARGET_SEQS = 100

class BlastError(RuntimeError):
    """Raised when a BLAST+ executable is missing or exits non-zero."""


@dataclass(frozen=True)
class BlastHit:
    """One HSP, with coordinates normalised to always run low -> high."""

    query_id: str
    subject_id: str
    pident: float
    aln_length: int
    mismatch: int
    gapopen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    qlen: int
    qcovhsp: float


    @classmethod
    def from_row(cls, row: str) -> "BlastHit":
        f = row.rstrip("\n").split("\t")
        if len(f) != len(OUTFMT_FIELDS):
            raise BlastError(
                f"expected {len(OUTFMT_FIELDS)} blastn columns, got {len(f)}: {row!r}"
            )
        return cls(
            query_id=f[0],
            subject_id=f[1],
            pident=float(f[2]),
            aln_length=int(f[3]),
            mismatch=int(f[4]),
            gapopen=int(f[5]),
            qstart=int(f[6]),
            qend=int(f[7]),
            sstart=int(f[8]),
            send=int(f[9]),
            evalue=float(f[10]),
            bitscore=float(f[11]),
            qlen=int(f[12]),
            qcovhsp=float(f[13]),
        )

    @property
    def strand(self) -> str:
        """'+' or '-', from the orientation of the subject coordinates."""
        return "-" if self.sstart > self.send else "+"

    @property
    def ref_start(self) -> int:
        """1-based start on the reference, always <= ref_end."""
        return min(self.sstart, self.send)

    @property
    def ref_end(self) -> int:
        return max(self.sstart, self.send)

    def as_record(self) -> dict:
        rec = asdict(self)
        rec.update(strand=self.strand, ref_start=self.ref_start, ref_end=self.ref_end)
        return rec


def _require(exe: str) -> str:
    path = shutil.which(exe)
    if path is None:
        raise BlastError(f"{exe} not found on PATH; install BLAST+ (e.g. conda install -c bioconda blast)")
    return path


def make_blast_db(fasta: str, db_dir: str, force: bool = False) -> str:
    """Build a nucleotide BLAST database for *fasta* inside *db_dir*.

    Returns the database prefix. Rebuilds only when the database is missing,
    older than the FASTA, or *force* is set.
    """
    if not os.path.isfile(fasta):
        raise FileNotFoundError(f"reference FASTA not found: {fasta}")
    if any(c.isspace() for c in db_dir):
        raise BlastError(
            f"BLAST cannot read a database path containing whitespace: {db_dir!r}; "
            "pass a whitespace-free --db-dir"
        )
    os.makedirs(db_dir, exist_ok=True)
    db_prefix = os.path.join(db_dir, "_".join(os.path.basename(fasta).split()))

    index = db_prefix + ".nin"
    if not force and os.path.exists(index) and os.path.getmtime(index) >= os.path.getmtime(fasta):
        return db_prefix

    # BLAST+ parses -in as a space-separated list of files, so a path containing
    # whitespace has to reach it through a whitespace-free symlink.
    fasta_in = os.path.abspath(fasta)
    if any(c.isspace() for c in fasta_in):
        fasta_in = db_prefix + ".source.fasta"
        if os.path.islink(fasta_in) or os.path.exists(fasta_in):
            os.unlink(fasta_in)
        os.symlink(os.path.abspath(fasta), fasta_in)

    proc = subprocess.run(
        [_require("makeblastdb"), "-in", fasta_in, "-dbtype", "nucl", "-out", db_prefix],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise BlastError(f"makeblastdb failed:\n{proc.stderr.strip()}")
    return db_prefix



def write_fasta(records: Sequence[tuple[str, str]], handle) -> int:
    written = 0
    for seq_id, seq in records:
        if not seq:
            continue
        handle.write(f">{seq_id}\n{seq}\n")
        written += 1
    handle.flush()
    return written



def run_blastn(
    records: Iterable[tuple[str, str]],
    db: str,
    task: str = DEFAULT_TASK,
    word_size: int = DEFAULT_WORD_SIZE,
    evalue: float = DEFAULT_EVALUE,
    max_target_seqs: int = DEFAULT_MAX_TARGET_SEQS,
    threads: int = 1,
    dust: str = "no",
    extra_args: Sequence[str] | None = None,
) -> list[BlastHit]:
    """Blast every ``(query_id, sequence)`` in *records* against *db*.

    All sequences go into one blastn process, so query IDs must be unique and
    free of whitespace. Returns every HSP reported, in blastn's output order.
    """
    records = [(str(i), s) for i, s in records]
    if not records:
        return []

    ids = [i for i, _ in records]
    if len(set(ids)) != len(ids):
        raise ValueError("query IDs passed to run_blastn must be unique")

    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as fh:
        query_path = fh.name
        n = write_fasta(records, fh)
    try:
        if n == 0:
            return []
        cmd = [
            _require("blastn"),
            "-query", query_path,
            "-db", db,
            "-task", task,
            "-word_size", str(word_size),
            "-evalue", str(evalue),
            "-max_target_seqs", str(max_target_seqs),
            "-num_threads", str(threads),
            "-dust", dust,
            "-outfmt", "6 " + " ".join(OUTFMT_FIELDS),
        ]
        if extra_args:
            cmd.extend(extra_args)

        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise BlastError(f"blastn failed:\n{proc.stderr.strip()}")
        return [BlastHit.from_row(line) for line in proc.stdout.splitlines() if line.strip()]
    finally:
        os.unlink(query_path)


def blast_one(seq_id: str, sequence: str, db: str, **kwargs) -> list[BlastHit]:
    """Blast a single sequence. Convenience wrapper over :func:`run_blastn`."""
    return run_blastn([(seq_id, sequence)], db, **kwargs)


def best_hit_per_query(hits: Iterable[BlastHit]) -> dict[str, BlastHit]:
    """Pick one hit per query: highest bitscore, then lowest e-value, then longest."""
    best: dict[str, BlastHit] = {}
    for hit in hits:
        current = best.get(hit.query_id)
        if current is None or (
            (hit.bitscore, -hit.evalue, hit.aln_length)
            > (current.bitscore, -current.evalue, current.aln_length)
        ):
            best[hit.query_id] = hit
    return best


def count_hits_per_query(hits: Iterable[BlastHit]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for hit in hits:
        counts[hit.query_id] = counts.get(hit.query_id, 0) + 1
    return counts