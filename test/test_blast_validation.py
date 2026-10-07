"""Independent check of the classification with blastn instead of bwa.

    python -m unittest discover -s test -p 'test_blast_validation.py' -v

A third opinion on the ``NM=0`` hit sets. ``classify.exact_hit_sets`` already
validates them by substring search, but that oracle and the classifier share this
repo's assumptions; blastn shares none of them. If the two agree, a hit set is
wrong only if bwa, a substring search and blastn are wrong in the same way.

Scope: ``NM=0`` only, matching the rest of the pipeline. Per-tile labels and the
permissive ``NM<=n`` columns are not checked here (see ``src/classify.py``).

Short-query blast has three traps, all of which silently *lose* true hits:

``-task blastn-short``
    Default blastn uses a word size of 11 and will not reliably find a 30 bp
    match at all. ``blastn-short`` drops it to 7.
``-dust no -soft_masking false``
    DUST masks low-complexity query sequence by default. The PIK3CG tandem repeat
    around nt 5950-6100 is exactly that, so leaving masking on drops real hits
    and manufactures a disagreement that looks like a classifier bug.
``length == 30 and mismatch == 0 and gapopen == 0``
    blastn reports *any* seeded HSP, including 8-9 bp fragments -- the full run
    emits ~2.8M HSPs for 28k tiles, of which ~28k are full-length. Filtering on
    ``pident`` alone is not enough, because a 9 bp alignment is also 100%
    identical. The length filter is what makes this equivalent to NM=0.

Plus-strand only, to match the pipeline: the reference is mRNA sense sequence, so
a minus-strand hit is not a site a guide against the sense strand would bind.

By default a deterministic sample of tiles is checked, because the full set takes
about 85 s. To check every tile::

    PIK3_BLAST_SAMPLE=all python -m unittest discover -s test \
        -p 'test_blast_validation.py' -v
"""

from __future__ import annotations

import csv
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "src"))

import classify as classify_mod  # noqa: E402
import tile as tile_mod  # noqa: E402

DATA_DIR = os.path.join(_ROOT, "data")
TILES_FASTA = os.path.join(_ROOT, "results", "tiles.fasta")
CLASSIFICATION = os.path.join(_ROOT, "results", "classification.tsv")

#: Tiles to check unless PIK3_BLAST_SAMPLE says otherwise. 400 keeps the test
#: under a few seconds while still covering every transcript.
DEFAULT_SAMPLE = 400

#: Fields blastn is asked for, in order.
OUTFMT = "6 qseqid sseqid pident length mismatch gapopen sstrand"

K = 30


def sample_size() -> int | None:
    """``None`` means every tile. Set ``PIK3_BLAST_SAMPLE=all`` for the full run."""
    requested = os.environ.get("PIK3_BLAST_SAMPLE", "").strip().lower()
    if not requested:
        return DEFAULT_SAMPLE
    if requested == "all":
        return None
    return int(requested)


@unittest.skipUnless(
    shutil.which("blastn") and shutil.which("makeblastdb"), "blastn/makeblastdb not on PATH"
)
@unittest.skipUnless(
    os.path.isfile(TILES_FASTA) and os.path.isfile(CLASSIFICATION),
    "run tile_transcripts.py and classify_guides.py first",
)
class TestBlastAgreesWithClassification(unittest.TestCase):
    """blastn must reproduce the NM=0 hit sets in results/classification.tsv."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()

        records = [
            record
            for path in tile_mod.find_fasta_files(DATA_DIR)
            for record in tile_mod.read_fasta(path)
        ]
        cls.genes = classify_mod.gene_symbols(records)
        cls.sequences = {record.transcript_id: record.sequence for record in records}

        with open(TILES_FASTA, encoding="utf-8") as handle:
            tiles = {r.transcript_id: r.sequence for r in tile_mod.parse_fasta(handle)}
        with open(CLASSIFICATION, newline="", encoding="utf-8") as handle:
            rows = {r["tile_id"]: r for r in csv.DictReader(handle, delimiter="\t")}

        missing = sorted(set(rows) - set(tiles))
        assert not missing, f"classification has tiles absent from tiles.fasta: {missing[:3]}"

        size = sample_size()
        chosen = sorted(rows)
        if size is not None and size < len(chosen):
            # Seeded, so a failure is reproducible rather than a different sample
            # on every run.
            chosen = sorted(random.Random(0).sample(chosen, size))
        cls.rows = {tile_id: rows[tile_id] for tile_id in chosen}
        cls.tiles = {tile_id: tiles[tile_id] for tile_id in chosen}

        cls.blast_hits = cls.run_blast(cls.tiles)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    @classmethod
    def run_blast(cls, tiles: dict[str, str]) -> dict[str, set[str]]:
        """blastn every tile against the transcripts; return gene sets of exact hits."""
        work = cls.tmp.name
        reference = os.path.join(work, "ref.fasta")
        with open(reference, "w", encoding="utf-8") as handle:
            for transcript_id, sequence in cls.sequences.items():
                handle.write(f">{transcript_id}\n{sequence}\n")

        query = os.path.join(work, "query.fasta")
        with open(query, "w", encoding="utf-8") as handle:
            for tile_id, sequence in tiles.items():
                handle.write(f">{tile_id}\n{sequence}\n")

        database = os.path.join(work, "refdb")
        cls._run_command(["makeblastdb", "-in", reference, "-dbtype", "nucl", "-out", database])
        completed = cls._run_command(
            [
                "blastn",
                "-task", "blastn-short",   # word size 7; default 11 misses 30-mers
                "-query", query,
                "-db", database,
                "-dust", "no",             # never mask the query: see module docstring
                "-soft_masking", "false",
                "-perc_identity", "100",
                "-evalue", "1000",         # short exact hits still score low
                "-max_target_seqs", "100",
                "-outfmt", OUTFMT,
            ],
            capture=True,
        )

        hits: dict[str, set[str]] = {}
        for line in completed.stdout.decode().splitlines():
            if not line.strip():
                continue
            query_id, subject, pident, length, mismatch, gapopen, strand = line.split("\t")
            # Full-length, ungapped, no substitutions: the blastn equivalent of NM=0.
            # pident alone would accept a 9 bp fragment.
            if int(length) != K or int(mismatch) != 0 or int(gapopen) != 0:
                continue
            if float(pident) != 100.0 or strand != "plus":
                continue
            hits.setdefault(query_id, set()).add(cls.genes[subject])
        return hits

    @staticmethod
    def _run_command(command: list[str], capture: bool = False):
        # Not `run`: that would shadow TestCase.run and unittest would call it
        # with a TestResult.
        completed = subprocess.run(
            command,
            stdout=subprocess.PIPE if capture else subprocess.DEVNULL,
            stderr=subprocess.PIPE,
        )
        if completed.returncode != 0:
            raise AssertionError(
                f"{command[0]} failed (exit {completed.returncode}): "
                f"{completed.stderr.decode('utf-8', 'replace').strip()}"
            )
        return completed

    def test_every_tile_has_a_full_length_exact_hit(self):
        """A tile was cut from a transcript, so blastn must find it there."""
        unhit = sorted(tile_id for tile_id in self.tiles if not self.blast_hits.get(tile_id))
        self.assertEqual(unhit, [], f"{len(unhit)} tile(s) with no exact blastn hit")

    def test_every_tile_hits_its_own_gene(self):
        for tile_id, row in self.rows.items():
            self.assertIn(
                row["gene"], self.blast_hits.get(tile_id, set()), f"{tile_id} missed its own gene"
            )

    def test_blast_reproduces_the_classification_hit_sets(self):
        """The headline check: identical gene sets, tile by tile."""
        disagreements = []
        for tile_id, row in self.rows.items():
            expected = set(filter(None, row["hit_genes_exact"].split(";")))
            found = self.blast_hits.get(tile_id, set())
            if expected != found:
                disagreements.append((tile_id, sorted(expected), sorted(found)))
        self.assertEqual(
            disagreements,
            [],
            f"{len(disagreements)} tile(s) where blastn and the classifier disagree, "
            f"e.g. {disagreements[:3]} (tile_id, classifier, blastn)",
        )

    def test_blast_confirms_every_hit_is_specific(self):
        """What the real data actually shows: one gene per tile, nothing shared.

        Asserted against blastn rather than our own output, so it is an independent
        statement about the sequences. If a future reference gains a conserved
        region this fails loudly rather than drifting -- see the finding note in
        CLAUDE.md.
        """
        shared = {
            tile_id: sorted(genes)
            for tile_id, genes in self.blast_hits.items()
            if len(genes) != 1
        }
        self.assertEqual(
            shared, {}, f"{len(shared)} tile(s) hit more than one gene, e.g. {list(shared)[:3]}"
        )
        self.assertEqual(
            {row["class_exact"] for row in self.rows.values()},
            {"specific"},
            "classification reports a class other than 'specific'",
        )

    def test_blast_agrees_with_the_substring_oracle(self):
        """Close the triangle: blastn vs the oracle directly, not via the TSV."""
        oracle = classify_mod.exact_hit_sets(self.tiles, self.sequences)
        for tile_id in self.tiles:
            expected = {self.genes[ref] for ref in oracle[tile_id]}
            self.assertEqual(self.blast_hits.get(tile_id, set()), expected, tile_id)


if __name__ == "__main__":
    unittest.main()
