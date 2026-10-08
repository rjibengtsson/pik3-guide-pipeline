"""Tests for src/classify.py and scripts/classify_guides.py.

    python -m unittest discover -s test -v

The CLI tests shell out to bwa and are skipped when it is not on PATH.
"""

from __future__ import annotations

import csv
import io
import os
import random
import shutil
import sys
import tempfile
import unittest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "src"))
sys.path.insert(0, os.path.join(_ROOT, "scripts"))

import classify as classify_mod  # noqa: E402
import classify_guides  # noqa: E402
import tile as tile_mod  # noqa: E402
import tile_transcripts  # noqa: E402

ALL_GENES = frozenset({"PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG"})

# Regions planted so the right answer is known by construction. Each is 40 nt, so
# with k=30 there are 11 tiles lying wholly inside one and the flanking filler
# cannot leak into them.
PAN = "ACGTTGCAGGTTCCAAGTACGTAAGGCCTTACGTTGCAGT"
DUAL_DG = "TTGGCCAATTCCGGAATTCCGGTTAACCGGTTAACCGGAA"
PAIR_AB = "GGAATTCCAAGGTTCCTTAAGGCCAATTGGCCTTAAGGCC"
TRIPLE_ABD = "CCTTAAGGCCAATTGGAACCTTGGAACCTTGGCCAATTGG"

def _filler(seed: int, length: int = 60) -> str:
    """Deterministic gene-unique filler.

    Pseudo-random rather than a tandem repeat: a repeat would give bwa many
    equally-best hits and so exercise its multi-hit handling (see
    ``TestCliLowComplexity``) instead of the classification logic under test here.
    Two independent 30-mers colliding by chance is a 4**-30 event.
    """
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


#: Distinct per-gene filler: no 30-mer of one appears in another.
FILLER = {gene: _filler(index) for index, gene in enumerate(sorted(ALL_GENES))}

TRANSCRIPTS = {
    "NM_000001.1": ("PIK3CA", FILLER["PIK3CA"] + PAN + PAIR_AB + TRIPLE_ABD + FILLER["PIK3CA"]),
    "NM_000002.1": ("PIK3CB", FILLER["PIK3CB"] + PAN + PAIR_AB + TRIPLE_ABD + FILLER["PIK3CB"]),
    "NM_000003.1": ("PIK3CD", FILLER["PIK3CD"] + PAN + DUAL_DG + TRIPLE_ABD + FILLER["PIK3CD"]),
    "NM_000004.1": ("PIK3CG", FILLER["PIK3CG"] + PAN + DUAL_DG + FILLER["PIK3CG"]),
}


class TestSourceTranscript(unittest.TestCase):
    def test_splits_off_the_coordinates(self):
        self.assertEqual(classify_mod.source_transcript("NM_006218.4:1-30"), "NM_006218.4")

    def test_splits_on_the_last_colon(self):
        self.assertEqual(classify_mod.source_transcript("weird:id:71-100"), "weird:id")

    def test_rejects_non_tile_ids(self):
        for bad in ("NM_006218.4", "", ":1-30", "NM_006218.4:"):
            with self.subTest(tile_id=bad), self.assertRaises(ValueError):
                classify_mod.source_transcript(bad)


class TestGeneSymbols(unittest.TestCase):
    def record(self, transcript_id: str, description: str):
        return tile_mod.Record(transcript_id, description, "ACGT", "f.fasta")

    def test_parses_the_parenthesised_symbol(self):
        records = [
            self.record("NM_006218.4", "Homo sapiens ... subunit alpha (PIK3CA), mRNA"),
            self.record("NM_002649.3", "Homo sapiens ... gamma (PIK3CG), transcript variant 1, mRNA"),
        ]
        self.assertEqual(
            classify_mod.gene_symbols(records),
            {"NM_006218.4": "PIK3CA", "NM_002649.3": "PIK3CG"},
        )

    def test_no_symbol_is_an_error_rather_than_a_fallback(self):
        with self.assertRaises(ValueError):
            classify_mod.gene_symbols([self.record("NM_1.1", "no symbol here, mRNA")])

    def test_ambiguous_multiple_symbols_is_an_error(self):
        with self.assertRaises(ValueError):
            classify_mod.gene_symbols([self.record("NM_1.1", "a (PIK3CA) b (PIK3CB)")])

    def test_real_data_headers_parse(self):
        """The actual PI3K FASTA headers must yield one symbol each.

        Inputs are per dataset -- ``data/pik3-iso`` is the four-paralog worked
        example -- so this names that folder rather than scanning ``data/``, which
        also holds unrelated sets whose symbols are not ALL_GENES.
        """
        data_dir = os.path.join(_ROOT, "data", "pik3-iso")
        records = [
            record
            for path in tile_mod.find_fasta_files(data_dir)
            for record in tile_mod.read_fasta(path)
        ]
        self.assertEqual(set(classify_mod.gene_symbols(records).values()), ALL_GENES)


class TestClassify(unittest.TestCase):
    def label(self, *genes: str) -> str:
        return classify_mod.classify(frozenset(genes), ALL_GENES)

    def test_one_gene_is_specific(self):
        for gene in sorted(ALL_GENES):
            with self.subTest(gene=gene):
                self.assertEqual(self.label(gene), "specific")

    def test_the_designed_pair_is_dual(self):
        self.assertEqual(self.label("PIK3CD", "PIK3CG"), "dual")

    def test_every_gene_is_pan(self):
        self.assertEqual(self.label(*ALL_GENES), "pan")

    def test_any_other_pair_is_other_not_dual(self):
        for pair in (
            ("PIK3CA", "PIK3CB"),
            ("PIK3CA", "PIK3CD"),
            ("PIK3CA", "PIK3CG"),
            ("PIK3CB", "PIK3CD"),
            ("PIK3CB", "PIK3CG"),
        ):
            with self.subTest(pair=pair):
                self.assertEqual(self.label(*pair), "other")

    def test_three_genes_is_other(self):
        self.assertEqual(self.label("PIK3CA", "PIK3CB", "PIK3CD"), "other")

    def test_no_hits_is_unmapped(self):
        self.assertEqual(self.label(), "unmapped")

    def test_hits_outside_the_reference_are_an_error(self):
        with self.assertRaises(ValueError):
            classify_mod.classify(frozenset({"KRAS"}), ALL_GENES)

    def test_pan_wins_over_dual_on_a_two_gene_reference(self):
        """Documented precedence: the stronger claim about the whole reference."""
        self.assertEqual(
            classify_mod.classify(classify_mod.DUAL_PAIR, classify_mod.DUAL_PAIR), "pan"
        )


class TestHitSets(unittest.TestCase):
    def test_collapses_to_one_set_per_tile(self):
        hits = [
            classify_mod.Hit("t1", "NM_1.1", 0),
            classify_mod.Hit("t1", "NM_2.1", 0),
            classify_mod.Hit("t1", "NM_1.1", 0),  # duplicate reference
            classify_mod.Hit("t2", "NM_1.1", 0),
        ]
        sets, _, _ = classify_mod.hit_sets(hits)
        self.assertEqual(sets["t1"], frozenset({"NM_1.1", "NM_2.1"}))
        self.assertEqual(sets["t2"], frozenset({"NM_1.1"}))

    def test_mismatch_cutoff_is_applied(self):
        hits = [
            classify_mod.Hit("t1", "NM_1.1", 0),
            classify_mod.Hit("t1", "NM_2.1", 2),
            classify_mod.Hit("t1", "NM_3.1", 4),
        ]
        exact, _, _ = classify_mod.hit_sets(hits, max_nm=0)
        permissive, worst, _ = classify_mod.hit_sets(hits, max_nm=3)
        self.assertEqual(exact["t1"], frozenset({"NM_1.1"}))
        self.assertEqual(permissive["t1"], frozenset({"NM_1.1", "NM_2.1"}))
        self.assertEqual(worst["t1"], 2)

    def test_reverse_strand_hits_are_dropped_and_counted(self):
        hits = [
            classify_mod.Hit("t1", "NM_1.1", 0),
            classify_mod.Hit("t1", "NM_2.1", 0, is_reverse=True),
        ]
        sets, _, dropped = classify_mod.hit_sets(hits)
        self.assertEqual(sets["t1"], frozenset({"NM_1.1"}))
        self.assertEqual(dropped, 1)

    def test_negative_cutoff_rejected(self):
        with self.assertRaises(ValueError):
            classify_mod.hit_sets([], max_nm=-1)


class TestHitsFromSam(unittest.TestCase):
    """Reading hits out of SAM, including bwa's awkward multi-hit records."""

    HEADER = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:NM_1.1\tLN:240\n@SQ\tSN:NM_2.1\tLN:240\n"
    SEQ = "A" * 30

    def hits(self, body: str) -> list[classify_mod.Hit]:
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "x.sam")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(self.HEADER + body)
            return list(classify_mod.hits_from_sam(path))

    def row(self, name: str, flag: int, ref: str, pos: int, tags: str) -> str:
        return f"{name}\t{flag}\t{ref}\t{pos}\t0\t30M\t*\t0\t0\t{self.SEQ}\t*\t{tags}\n"

    def test_reads_primary_and_every_xa_entry(self):
        hits = self.hits(
            self.row("t1", 0, "NM_1.1", 1, "NM:i:0\tXA:Z:NM_2.1,+61,30M,2;NM_1.1,+99,30M,0;")
        )
        self.assertEqual(
            sorted((h.reference, h.nm) for h in hits),
            [("NM_1.1", 0), ("NM_1.1", 0), ("NM_2.1", 2)],
        )

    def test_xa_entries_keep_their_own_nm_not_the_primarys(self):
        """An XA hit is usually worse than the primary; inheriting NM would promote it."""
        hits = self.hits(self.row("t1", 0, "NM_1.1", 1, "NM:i:0\tXA:Z:NM_2.1,+61,30M,3;"))
        by_reference = {h.reference: h.nm for h in hits}
        self.assertEqual(by_reference, {"NM_1.1": 0, "NM_2.1": 3})

    def test_unmapped_flag_with_a_reported_alignment_is_still_a_hit(self):
        """bwa sets flag 0x4 on reads with many equally-best hits, yet fills in
        RNAME/POS/NM/XA. Trusting is_unmapped there silently drops perfect
        alignments -- including a tile's own source transcript."""
        hits = self.hits(
            self.row("t1", 4, "NM_1.1", 211, "NM:i:0\tXA:Z:NM_1.1,+1,30M,0;NM_2.1,+5,30M,0;")
        )
        self.assertEqual({h.reference for h in hits}, {"NM_1.1", "NM_2.1"})
        self.assertTrue(all(h.nm == 0 for h in hits))

    def test_truly_unaligned_read_yields_nothing(self):
        self.assertEqual(self.hits(f"t1\t4\t*\t0\t0\t*\t*\t0\t0\t{self.SEQ}\t*\n"), [])

    def test_reverse_strand_is_taken_from_the_xa_sign(self):
        hits = self.hits(self.row("t1", 0, "NM_1.1", 1, "NM:i:0\tXA:Z:NM_2.1,-61,30M,0;"))
        reverse = {h.reference: h.is_reverse for h in hits}
        self.assertEqual(reverse, {"NM_1.1": False, "NM_2.1": True})

    def test_malformed_xa_entry_is_an_error(self):
        with self.assertRaises(ValueError):
            self.hits(self.row("t1", 0, "NM_1.1", 1, "NM:i:0\tXA:Z:NM_2.1,+61,30M;"))

    def test_positions_are_1_based_for_both_primary_and_xa(self):
        """POS is 1-based in SAM, pysam's reference_start is not, and XA is. All
        three must come out on the tile_id convention: 1-based inclusive."""
        hits = self.hits(self.row("t1", 0, "NM_1.1", 7, "NM:i:0\tXA:Z:NM_2.1,+61,30M,2;"))
        self.assertEqual(
            sorted((h.reference, h.position) for h in hits),
            [("NM_1.1", 7), ("NM_2.1", 61)],
        )

    def test_reverse_xa_position_keeps_its_magnitude(self):
        hits = self.hits(self.row("t1", 0, "NM_1.1", 1, "NM:i:0\tXA:Z:NM_2.1,-61,30M,0;"))
        reverse = [h for h in hits if h.is_reverse]
        self.assertEqual([h.position for h in reverse], [61])

    def test_alignment_running_past_the_reference_end_is_flagged(self):
        """bwa indexes the references concatenated, so an XA entry can report an
        alignment crossing into the next transcript while still naming the previous
        one. Within that reference the match does not exist."""
        hits = self.hits(
            self.row("t1", 0, "NM_1.1", 1, "NM:i:0\tXA:Z:NM_1.1,+212,30M,0;NM_2.1,+211,30M,0;")
        )
        flagged = {(h.reference, h.position): h.overhangs for h in hits}
        # LN:240, so 212 + 30 - 1 = 241 is past the end; 211 ends exactly on it.
        self.assertTrue(flagged[("NM_1.1", 212)])
        self.assertFalse(flagged[("NM_2.1", 211)])
        self.assertFalse(flagged[("NM_1.1", 1)])

    def test_overhanging_primary_alignment_is_flagged_too(self):
        hits = self.hits(self.row("t1", 0, "NM_1.1", 212, "NM:i:0"))
        self.assertEqual([h.overhangs for h in hits], [True])


class TestReferenceSpan(unittest.TestCase):
    """Only reference-consuming CIGAR operations count towards a hit's end."""

    def test_match_only(self):
        self.assertEqual(classify_mod.reference_span("30M"), 30)

    def test_deletions_and_skips_consume_reference(self):
        """30 reference bases from a 25 base query, so the span is not the read length."""
        self.assertEqual(classify_mod.reference_span("10M5D15M"), 30)
        self.assertEqual(classify_mod.reference_span("10M5N15M"), 30)

    def test_insertions_and_clips_do_not(self):
        self.assertEqual(classify_mod.reference_span("10M5I15M"), 25)
        self.assertEqual(classify_mod.reference_span("5S25M"), 25)
        self.assertEqual(classify_mod.reference_span("5H25M"), 25)

    def test_malformed_cigar_is_an_error(self):
        for bad in ("", "30", "M", "30Z", "30M?"):
            with self.subTest(cigar=bad), self.assertRaises(ValueError):
                classify_mod.reference_span(bad)


class TestHitLoci(unittest.TestCase):
    """Positions survive per reference, which is what hit_sets throws away."""

    def test_keeps_every_position_on_the_same_reference(self):
        hits = [
            classify_mod.Hit("t1", "a", 0, position=5),
            classify_mod.Hit("t1", "a", 0, position=100),
            classify_mod.Hit("t1", "b", 0, position=7),
        ]
        self.assertEqual(
            classify_mod.hit_loci(hits),
            {"t1": {"a": [(5, 0), (100, 0)], "b": [(7, 0)]}},
        )

    def test_identical_hits_are_deduplicated(self):
        """bwa can report the primary alignment again as an XA entry."""
        hits = [classify_mod.Hit("t1", "a", 0, position=5)] * 3
        self.assertEqual(classify_mod.hit_loci(hits), {"t1": {"a": [(5, 0)]}})

    def test_mismatch_cutoff_matches_hit_sets(self):
        hits = [
            classify_mod.Hit("t1", "a", 0, position=5),
            classify_mod.Hit("t1", "b", 3, position=9),
        ]
        self.assertEqual(classify_mod.hit_loci(hits, max_nm=0), {"t1": {"a": [(5, 0)]}})
        loci = classify_mod.hit_loci(hits, max_nm=3)
        self.assertEqual(set(loci["t1"]), set(classify_mod.hit_sets(hits, max_nm=3)[0]["t1"]))

    def test_reverse_and_overhanging_hits_are_dropped(self):
        hits = [
            classify_mod.Hit("t1", "a", 0, position=5),
            classify_mod.Hit("t1", "b", 0, is_reverse=True, position=9),
            classify_mod.Hit("t1", "c", 0, position=11, overhangs=True),
        ]
        self.assertEqual(classify_mod.hit_loci(hits), {"t1": {"a": [(5, 0)]}})

    def test_negative_cutoff_rejected(self):
        with self.assertRaises(ValueError):
            classify_mod.hit_loci([], max_nm=-1)


class TestFormatLoci(unittest.TestCase):
    def test_renders_accession_coordinates_and_nm(self):
        self.assertEqual(
            classify_mod.format_loci({"NM_1.1": [(2724, 0)]}, 30),
            "NM_1.1:2724-2753(0)",
        )

    def test_semicolon_between_references_comma_within_one(self):
        rendered = classify_mod.format_loci(
            {"NM_2.1": [(5, 3)], "NM_1.1": [(1, 0), (61, 2)]}, 30
        )
        self.assertEqual(rendered, "NM_1.1:1-30(0),NM_1.1:61-90(2);NM_2.1:5-34(3)")

    def test_no_hits_is_empty(self):
        self.assertEqual(classify_mod.format_loci({}, 30), "")

    def test_coordinates_are_1_based_inclusive_like_tile_ids(self):
        """A locus string must be readable as a coordinate on the reference, so a
        tile taken from position 1 of its own transcript renders as its own id."""
        self.assertEqual(
            classify_mod.format_loci({"NM_1.1": [(1, 0)]}, 30), "NM_1.1:1-30(0)"
        )


class TestExactLocusOracle(unittest.TestCase):
    """The position-aware oracle, independent of the aligner like its sibling."""

    def test_finds_every_occurrence_including_overlapping_ones(self):
        found = classify_mod.exact_hit_loci({"t": "AA"}, {"a": "AAAA"})
        self.assertEqual(found["t"], {"a": [1, 2, 3]})

    def test_positions_are_1_based(self):
        found = classify_mod.exact_hit_loci({"t": "CCCC"}, {"a": "AAAACCCCGGGG"})
        self.assertEqual(found["t"], {"a": [5]})

    def test_reports_each_reference_separately(self):
        found = classify_mod.exact_hit_loci(
            {"t": "CCCC"}, {"a": "AAAACCCC", "b": "CCCCTT", "c": "GGGG"}
        )
        self.assertEqual(found["t"], {"a": [5], "b": [1]})

    def test_is_case_insensitive(self):
        found = classify_mod.exact_hit_loci({"t": "acgt"}, {"a": "TTACGTTT"})
        self.assertEqual(found["t"], {"a": [3]})

    def test_agrees_with_the_set_oracle_on_which_references_hit(self):
        references = {"a": "AAAACCCCGGGG", "b": "TTTTCCCCTTTT", "c": "GGGGGGGG"}
        tiles = {"t_cccc": "CCCC", "t_gggg": "GGGG", "t_none": "ACGTACGT"}
        loci = classify_mod.exact_hit_loci(tiles, references)
        sets = classify_mod.exact_hit_sets(tiles, references)
        for tile_id in tiles:
            self.assertEqual(frozenset(loci[tile_id]), sets[tile_id])


class TestExactOracle(unittest.TestCase):
    def test_finds_every_perfect_match(self):
        references = {"a": "AAAACCCCGGGG", "b": "TTTTCCCCTTTT", "c": "GGGGGGGG"}
        tiles = {"t_cccc": "CCCC", "t_aaaa": "AAAA", "t_gggg": "GGGG", "t_none": "ACGTACGT"}
        found = classify_mod.exact_hit_sets(tiles, references)
        self.assertEqual(found["t_cccc"], frozenset({"a", "b"}))
        self.assertEqual(found["t_aaaa"], frozenset({"a"}))
        self.assertEqual(found["t_gggg"], frozenset({"a", "c"}))
        self.assertEqual(found["t_none"], frozenset())

    def test_is_case_insensitive(self):
        found = classify_mod.exact_hit_sets({"t": "acgt"}, {"a": "TTACGTTT"})
        self.assertEqual(found["t"], frozenset({"a"}))

    def test_does_not_find_reverse_complements(self):
        """Forward strand only -- the reference is mRNA sense sequence."""
        found = classify_mod.exact_hit_sets({"t": "AAAA"}, {"a": "GGTTTTGG"})
        self.assertEqual(found["t"], frozenset())


class TestCheckSelfHits(unittest.TestCase):
    def test_passes_when_every_tile_hits_its_own_transcript(self):
        self.assertEqual(
            classify_mod.check_self_hits(
                {
                    "NM_1.1:1-30": frozenset({"NM_1.1"}),
                    "NM_2.1:1-30": frozenset({"NM_1.1", "NM_2.1"}),
                }
            ),
            [],
        )

    def test_reports_tiles_missing_their_own_transcript(self):
        self.assertEqual(
            classify_mod.check_self_hits(
                {
                    "NM_1.1:1-30": frozenset({"NM_2.1"}),
                    "NM_2.1:1-30": frozenset(),
                    "NM_3.1:1-30": frozenset({"NM_3.1"}),
                }
            ),
            ["NM_1.1:1-30", "NM_2.1:1-30"],
        )


@unittest.skipUnless(shutil.which("bwa"), "bwa not on PATH")
class TestCli(unittest.TestCase):
    """End to end: tile the synthetic family, then classify it.

    Runs the real stage-1 script so the tile ids, FASTA layout and self-hit
    invariant are exercised exactly as in production rather than hand-faked.
    """

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.input_dir = os.path.join(cls.tmp.name, "data")
        os.makedirs(cls.input_dir)
        for transcript_id, (gene, sequence) in TRANSCRIPTS.items():
            path = os.path.join(cls.input_dir, f"{transcript_id}.fasta")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(
                    f">{transcript_id} Homo sapiens synthetic subunit ({gene}), mRNA\n"
                    f"{sequence}\n"
                )

        cls.tiles_fasta = os.path.join(cls.tmp.name, "results", "tiles.fasta")
        cls.guides_fasta = os.path.join(cls.tmp.name, "results", "guides.fasta")
        cls.tiles_tsv = os.path.join(cls.tmp.name, "results", "tiles.tsv")
        cls.run_quietly(
            tile_transcripts.main,
            [
                "--input-dir", cls.input_dir,
                "--output", cls.tiles_tsv,
                "--input-type", "mrna",
                "-k", "30",
                "--fasta-out", cls.tiles_fasta,
                "--fasta-field", "target_seq",
            ],
        )
        cls.run_quietly(
            tile_transcripts.main,
            [
                "--input-dir", cls.input_dir,
                "--output", cls.tiles_tsv,
                "--input-type", "mrna",
                "-k", "30",
                "--fasta-out", cls.guides_fasta,
            ],
        )

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    @staticmethod
    def run_quietly(main, argv: list[str]) -> str:
        """Call a main(), capturing its stderr summary."""
        stderr = io.StringIO()
        old = sys.stderr
        sys.stderr = stderr
        try:
            assert main(argv) == 0
        finally:
            sys.stderr = old
        return stderr.getvalue()

    def classify(self, *extra: str, tiles: str | None = None) -> list[dict[str, str]]:
        output = os.path.join(self.tmp.name, "results", "classification.tsv")
        self.summary = self.run_quietly(
            classify_guides.main,
            [
                "--tiles", tiles or self.tiles_fasta,
                "--input-dir", self.input_dir,
                "--output", output,
                *extra,
            ],
        )
        with open(output, newline="", encoding="utf-8") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_every_tile_gets_exactly_one_row(self):
        rows = self.classify("--max-mismatches", "0")
        with open(self.tiles_fasta, encoding="utf-8") as handle:
            tile_ids = [r.transcript_id for r in tile_mod.parse_fasta(handle)]
        self.assertEqual([r["tile_id"] for r in rows], tile_ids)
        self.assertEqual(len(set(r["tile_id"] for r in rows)), len(tile_ids))

    def test_oracle_validation_passes_on_real_bwa_output(self):
        """The headline check: bwa's NM=0 hit sets equal the substring oracle's."""
        self.classify("--max-mismatches", "0")
        self.assertIn("oracle check: passed", self.summary)

    def test_planted_regions_get_the_expected_classes(self):
        rows = self.classify("--max-mismatches", "0")
        by_class: dict[str, set[str]] = {}
        for row in rows:
            by_class.setdefault(row["class_exact"], set()).add(row["hit_genes_exact"])

        self.assertIn(";".join(sorted(ALL_GENES)), by_class["pan"])
        self.assertEqual(by_class["dual"], {"PIK3CD;PIK3CG"})
        self.assertEqual(by_class["specific"], {gene for gene in ALL_GENES})
        self.assertIn("PIK3CA;PIK3CB", by_class["other"])
        self.assertIn("PIK3CA;PIK3CB;PIK3CD", by_class["other"])
        self.assertNotIn("unmapped", by_class)

    def test_pan_tile_count_matches_the_planted_region(self):
        """A 40 nt region shared by all four yields 40 - 30 + 1 tiles lying wholly
        inside it, *per transcript* -- same sequence, four distinct tile_ids each."""
        rows = self.classify("--max-mismatches", "0")
        pan = [r for r in rows if r["class_exact"] == "pan"]
        inside = len(PAN) - 30 + 1
        self.assertEqual(len(pan), inside * len(TRANSCRIPTS))
        self.assertEqual({r["n_hits_exact"] for r in pan}, {"4"})
        per_transcript = {transcript_id: 0 for transcript_id in TRANSCRIPTS}
        for row in pan:
            per_transcript[row["transcript_id"]] += 1
        self.assertEqual(set(per_transcript.values()), {inside})
        # Identical sequence throughout, so only `inside` distinct tiles exist.
        self.assertEqual(len({r["tile_id"].split(":")[1] for r in pan}), inside)

    def test_classes_agree_with_the_oracle_independently(self):
        """Recompute every label from the oracle alone and compare to the output."""
        rows = self.classify("--max-mismatches", "0")
        with open(self.tiles_fasta, encoding="utf-8") as handle:
            tiles = {r.transcript_id: r.sequence for r in tile_mod.parse_fasta(handle)}
        sequences = {tid: sequence for tid, (_, sequence) in TRANSCRIPTS.items()}
        genes = {tid: gene for tid, (gene, _) in TRANSCRIPTS.items()}
        oracle = classify_mod.exact_hit_sets(tiles, sequences)
        for row in rows:
            expected = classify_mod.classify(
                frozenset(genes[ref] for ref in oracle[row["tile_id"]]), ALL_GENES
            )
            self.assertEqual(row["class_exact"], expected, row["tile_id"])

    def test_every_tile_self_hits(self):
        rows = self.classify("--max-mismatches", "0")
        for row in rows:
            self.assertIn(row["gene"], row["hit_genes_exact"].split(";"), row["tile_id"])

    def test_permissive_columns_are_a_superset_of_exact(self):
        rows = self.classify("--max-mismatches", "3")
        for row in rows:
            exact = set(filter(None, row["hit_genes_exact"].split(";")))
            permissive = set(filter(None, row["hit_genes_nm"].split(";")))
            self.assertTrue(exact <= permissive, row["tile_id"])

    def test_locus_columns_name_the_same_references_as_the_gene_columns(self):
        """The two views of a hit set must not diverge: every accession in
        hit_locs_* must map back to a symbol in hit_genes_*, and vice versa."""
        rows = self.classify("--max-mismatches", "3")
        genes = {
            transcript_id: gene for transcript_id, (gene, _seq) in TRANSCRIPTS.items()
        }
        for row in rows:
            for locus_column, gene_column in (
                ("hit_locs_exact", "hit_genes_exact"),
                ("hit_locs_nm", "hit_genes_nm"),
            ):
                accessions = {
                    locus.rsplit(":", 1)[0]
                    for group in filter(None, row[locus_column].split(";"))
                    for locus in group.split(",")
                }
                self.assertEqual(
                    {genes[accession] for accession in accessions},
                    set(filter(None, row[gene_column].split(";"))),
                    f"{row['tile_id']} {locus_column}",
                )

    def test_every_tile_locus_includes_where_it_was_cut_from(self):
        """A tile's own coordinates must appear among its NM=0 loci -- the locus-level
        form of the self-hit invariant, and what pins the 1-based convention end to
        end: the tile_id itself is a reference coordinate."""
        for row in self.classify():
            self.assertIn(
                f"{row['tile_id']}(0)",
                row["hit_locs_exact"],
                row["tile_id"],
            )

    def test_locus_coordinates_round_trip_to_the_tile_sequence(self):
        """Slice the reference at each reported NM=0 locus and get the tile back."""
        tiles = {}
        with open(self.tiles_fasta, "r", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith(">"):
                    current = line[1:].strip()
                else:
                    tiles[current] = line.strip()
        sequences = {
            transcript_id: sequence for transcript_id, (_gene, sequence) in TRANSCRIPTS.items()
        }
        for row in self.classify():
            for group in filter(None, row["hit_locs_exact"].split(";")):
                for locus in group.split(","):
                    accession, span = locus.rsplit(":", 1)
                    coordinates, _, _nm = span.partition("(")
                    start, end = (int(value) for value in coordinates.split("-"))
                    self.assertEqual(
                        sequences[accession][start - 1 : end],
                        tiles[row["tile_id"]],
                        locus,
                    )

    def test_summary_flags_the_permissive_columns_as_unvalidated(self):
        self.classify("--max-mismatches", "3")
        self.assertIn("NOT validated", self.summary)

    def test_guides_fasta_is_rejected(self):
        """Guides are RNA and reverse-complemented; mapping them would be wrong."""
        with self.assertRaises(SystemExit) as caught:
            self.classify(tiles=self.guides_fasta)
        self.assertIn("looks like guides.fasta", str(caught.exception))

    def kept_work_dir(self, summary: str) -> str:
        marker = "kept intermediates in "
        line = next((l for l in summary.splitlines() if l.startswith(marker)), None)
        self.assertIsNotNone(line, f"summary did not report a kept work dir:\n{summary}")
        return line[len(marker):].strip()

    def test_default_run_discards_intermediates(self):
        self.classify("--max-mismatches", "0")
        self.assertIn("intermediates discarded", self.summary)
        self.assertNotIn("kept intermediates in", self.summary)

    def test_keep_temp_leaves_the_intermediates_to_inspect(self):
        self.classify("--max-mismatches", "0", "--keep-temp")
        work_dir = self.kept_work_dir(self.summary)
        self.addCleanup(shutil.rmtree, work_dir, True)
        for name in (
            "transcriptome.fasta",
            "transcriptome.fasta.bwt",
            "transcriptome.fasta.sa",
            "tiles.sai",
            "tiles.sam",
            "hits.tsv",
        ):
            self.assertTrue(os.path.isfile(os.path.join(work_dir, name)), name)

    def test_explicit_work_dir_is_kept_without_keep_temp(self):
        work_dir = os.path.join(self.tmp.name, "work")
        self.classify("--max-mismatches", "0", "--work-dir", work_dir)
        self.assertEqual(self.kept_work_dir(self.summary), work_dir)
        self.assertTrue(os.path.isfile(os.path.join(work_dir, "tiles.sam")))

    def test_dumped_hits_explain_the_hit_sets(self):
        """hits.tsv is pre-dedup, so it must account for every counted hit."""
        rows = self.classify("--max-mismatches", "3", "--keep-temp")
        work_dir = self.kept_work_dir(self.summary)
        self.addCleanup(shutil.rmtree, work_dir, True)
        with open(os.path.join(work_dir, "hits.tsv"), newline="", encoding="utf-8") as handle:
            hits = list(csv.DictReader(handle, delimiter="\t"))

        self.assertTrue(hits)
        # Reverse-strand hits must be visible but never counted.
        for hit in hits:
            if hit["strand"] == "-":
                self.assertEqual((hit["counted_exact"], hit["counted_nm"]), ("0", "0"))
            if int(hit["nm"]) > 0:
                self.assertEqual(hit["counted_exact"], "0")

        # Collapsing the dump by reference must reproduce the exact hit sets.
        genes = {tid: gene for tid, (gene, _) in TRANSCRIPTS.items()}
        collapsed: dict[str, set[str]] = {}
        for hit in hits:
            if hit["counted_exact"] == "1":
                collapsed.setdefault(hit["tile_id"], set()).add(genes[hit["reference"]])
        for row in rows:
            expected = set(filter(None, row["hit_genes_exact"].split(";")))
            self.assertEqual(collapsed.get(row["tile_id"], set()), expected, row["tile_id"])

    def test_mismatched_reference_is_rejected(self):
        other_dir = os.path.join(self.tmp.name, "other")
        os.makedirs(other_dir, exist_ok=True)
        with open(os.path.join(other_dir, "x.fasta"), "w", encoding="utf-8") as handle:
            handle.write(">NM_999.9 Homo sapiens unrelated (KRAS), mRNA\n" + "ACGT" * 20 + "\n")
        with self.assertRaises(SystemExit) as caught:
            self.classify("--input-dir", other_dir)
        self.assertIn("absent from", str(caught.exception))


@unittest.skipUnless(shutil.which("bwa"), "bwa not on PATH")
class TestCliLowComplexity(unittest.TestCase):
    """Regression: low-complexity tiles must not be lost.

    A poly-A tail and a tandem repeat give bwa many equally-best hits, and bwa
    then flags those records unmapped while still reporting valid alignments.
    Real mRNA has both, so this is not a synthetic-only concern.
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.input_dir = os.path.join(self.tmp.name, "data")
        os.makedirs(self.input_dir)
        shared = _filler(101, 40)
        transcripts = {
            "NM_100001.1": ("PIK3CA", _filler(11) + shared + "CAG" * 20 + "A" * 50),
            "NM_100002.1": ("PIK3CB", _filler(12) + shared + "CAG" * 20 + "A" * 50),
        }
        for transcript_id, (gene, sequence) in transcripts.items():
            path = os.path.join(self.input_dir, f"{transcript_id}.fasta")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(f">{transcript_id} synthetic ({gene}), mRNA\n{sequence}\n")

        self.tiles_fasta = os.path.join(self.tmp.name, "tiles.fasta")
        TestCli.run_quietly(
            tile_transcripts.main,
            [
                "--input-dir", self.input_dir,
                "--output", os.path.join(self.tmp.name, "tiles.tsv"),
                "--input-type", "mrna",
                "-k", "30",
                "--fasta-out", self.tiles_fasta,
                "--fasta-field", "target_seq",
            ],
        )

    def test_poly_a_and_repeat_tiles_still_self_hit_and_validate(self):
        output = os.path.join(self.tmp.name, "classification.tsv")
        summary = TestCli.run_quietly(
            classify_guides.main,
            [
                "--tiles", self.tiles_fasta,
                "--input-dir", self.input_dir,
                "--output", output,
                "--max-mismatches", "0",
            ],
        )
        self.assertIn("oracle check: passed", summary)
        with open(output, newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        # The poly-A and CAG tiles are shared, so nothing should be unmapped and
        # every tile must still hit its own transcript.
        self.assertFalse([r for r in rows if r["class_exact"] == "unmapped"])
        for row in rows:
            self.assertIn(row["gene"], row["hit_genes_exact"].split(";"), row["tile_id"])
        # A pure poly-A tile is in both transcripts, so it cannot be specific.
        poly_a = [r for r in rows if r["n_hits_exact"] == "2"]
        self.assertTrue(poly_a)


if __name__ == "__main__":
    unittest.main()
