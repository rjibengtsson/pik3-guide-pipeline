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

#: Inputs and outputs are per dataset: ``data/{DATASET}`` is tiled into
#: ``results/{DATASET}``. This check is specific to the PI3K family -- it asserts
#: one gene per tile over four paralogs -- so it names that dataset rather than
#: following whatever happens to be in ``results/``. Pointing it at another
#: dataset is these three lines plus the expectations below.
DATA_DIR = os.path.join(_ROOT, "data", "pik3-iso")
PPIB_DIR = os.path.join(_ROOT, "data", "ppib")
RESULTS_DIR = os.path.join(_ROOT, "results", "pik3")
TILES_FASTA = os.path.join(RESULTS_DIR, "tiles.fasta")
CLASSIFICATION = os.path.join(RESULTS_DIR, "classification.tsv")

#: The reference this check is pinned to: the four PI3K paralogs that
#: ``results/pik3/`` was generated from, one RefSeq variant each.
#:
#: Named accession by accession rather than discovered with
#: ``find_fasta_files(DATA_DIR)``, because that directory holds more than these
#: four and what else it holds has changed over time. It also carries other
#: transcript variants of the same four genes (``PIK3CA.fasta`` is
#: XM_054346848.1, not NM_006218.4) and gapped alignment output
#: (``PIK3CA_aligned.fasta`` is 27,176 nt of which 17,490 are ``-``). Globbing
#: either of those into the reference breaks this check in a way that reads as
#: a classifier bug: extra variants of a gene add hits the committed
#: classification never saw, and a gapped record is not sequence at all.
#:
#: The failure is also silent in the worst case -- ``PIK3CA.fasta`` and
#: ``PIK3CA_aligned.fasta`` share an accession, so a dict keyed on
#: ``transcript_id`` keeps whichever sorts last and the gapped one wins.
#:
#: PPIB rides along as a **negative control**, not as a fifth family member. It
#: is the housekeeping transcript the knockdown work measures against, so a PI3K
#: guide that also cut PPIB would invalidate its own control -- worth asserting
#: rather than assuming. It is deliberately *not* part of the tiled dataset: no
#: tile originates from it, so the self-hit invariant says nothing about it, and
#: its only job here is to stay absent from every hit set.
#:
#: Adding it cannot perturb the agreement with ``classification.tsv``, which was
#: produced over the four paralogs alone: no PI3K tile occurs in PPIB on either
#: strand, and the longest exact k-mer any of the 28,032 tiles shares with it is
#: 16 nt, half a spacer. Should that ever stop being true, the hit-set test is
#: what fails, and it fails loudly.
#:
#: Each accession carries its own directory because PPIB belongs to a different
#: dataset: ``data/ppib`` is tiled into ``results/ppib``, independently of this.
#:
#: Pointing this check at another dataset means these accessions, the
#: directories above, and the expectations below.
REFERENCE_FASTAS = {
    "NM_002649.3": DATA_DIR,   # PIK3CG
    "NM_005026.5": DATA_DIR,   # PIK3CD
    "NM_006218.4": DATA_DIR,   # PIK3CA
    "NM_006219.3": DATA_DIR,   # PIK3CB
    "NM_000942.5": PPIB_DIR,   # PPIB -- control, no tiles of its own
}

#: The genes a tile may legitimately hit. PPIB is in the reference but not in
#: here, so a hit on it is a failure rather than a new label.
FAMILY_GENES = frozenset({"PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG"})
CONTROL_GENES = frozenset({"PPIB"})

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
    """blastn must reproduce the NM=0 hit sets in results/pik3/classification.tsv."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()

        # Read by name, and fail rather than skip on a missing one: a silently
        # smaller reference would drop real cross-family hits and the
        # disagreement would land on the classifier.
        paths = {
            accession: os.path.join(directory, f"{accession}.fasta")
            for accession, directory in REFERENCE_FASTAS.items()
        }
        absent = sorted(p for p in paths.values() if not os.path.isfile(p))
        assert not absent, (
            "reference FASTA(s) missing: "
            + ", ".join(os.path.relpath(p, _ROOT) for p in absent)
        )

        records = [
            record for path in paths.values() for record in tile_mod.read_fasta(path)
        ]

        # The filename is not the accession: a file named NM_006218.4.fasta
        # holding some other record would give a reference that silently is
        # not what this check claims to validate.
        found = tuple(sorted(record.transcript_id for record in records))
        assert found == tuple(sorted(paths)), (
            f"reference holds {found}, expected {tuple(sorted(paths))} -- one "
            "of the reference FASTAs does not contain the accession its name "
            "claims"
        )

        cls.genes = classify_mod.gene_symbols(records)

        # The control has to actually be in the reference for its absence from
        # the hit sets to mean anything. Asserted rather than trusted, because
        # a typo in REFERENCE_FASTAS would otherwise turn the control into no
        # test at all -- it would simply never be hit.
        symbols = set(cls.genes.values())
        assert symbols == FAMILY_GENES | CONTROL_GENES, (
            f"reference spans {sorted(symbols)}, expected "
            f"{sorted(FAMILY_GENES | CONTROL_GENES)}"
        )
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

    def test_no_tile_hits_the_control_gene(self):
        """PPIB is in the reference so that this can be asserted, not assumed.

        PPIB is the housekeeping transcript the knockdown work normalises
        against, so a PI3K guide with an exact match in it would be measuring
        its own control away. It carries no tiles of its own: it is here purely
        to be missed.

        This is the one test in the file whose expected answer is an empty set,
        so `setUpClass` separately asserts PPIB is actually *in* the reference.
        Without that, a typo would make this pass by testing nothing.
        """
        hits = {
            tile_id: sorted(genes & CONTROL_GENES)
            for tile_id, genes in self.blast_hits.items()
            if genes & CONTROL_GENES
        }
        self.assertEqual(
            hits,
            {},
            f"{len(hits)} tile(s) exactly match a control gene, e.g. "
            f"{list(hits.items())[:3]}",
        )

    def test_blast_agrees_with_the_substring_oracle(self):
        """Close the triangle: blastn vs the oracle directly, not via the TSV."""
        oracle = classify_mod.exact_hit_sets(self.tiles, self.sequences)
        for tile_id in self.tiles:
            expected = {self.genes[ref] for ref in oracle[tile_id]}
            self.assertEqual(self.blast_hits.get(tile_id, set()), expected, tile_id)


if __name__ == "__main__":
    unittest.main()
