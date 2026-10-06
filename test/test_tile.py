"""Tests for src/tile.py and scripts/tile_transcripts.py.

    python -m unittest discover -s test -v
"""

from __future__ import annotations

import csv
import io
import os
import sys
import tempfile
import unittest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "src"))
sys.path.insert(0, os.path.join(_ROOT, "scripts"))

import tile as tile_mod  # noqa: E402
import tile_transcripts  # noqa: E402

def pairs_antiparallel(target: str, guide: str) -> bool:
    """True if ``guide`` (5' -> 3') is fully Watson-Crick paired to ``target`` (5' -> 3').

    A test helper, deliberately written without reusing ``tile.reverse_complement``
    so that it fails when that function is wrong instead of agreeing with it.
    Direction is checked explicitly: ``guide[i]`` must pair with the base ``i``
    positions from the *end* of ``target``. Ambiguity codes never count as paired.

    Alphabet-blind by design: both sides are read as RNA first, so a DNA-alphabet
    target (the default) and a ``--guide-alphabet dna`` guide both pass. The
    alphabet itself is covered by separate tests.
    """
    if len(target) != len(guide):
        return False
    pairing = {"A": "U", "U": "A", "G": "C", "C": "G"}
    target_rna = tile_mod.transcribe(target.upper())
    guide_rna = tile_mod.transcribe(guide.upper())
    return all(
        pairing.get(guide_rna[i]) == target_rna[len(target_rna) - 1 - i]
        for i in range(len(guide_rna))
    )


TWO_RECORD_FASTA = """\
>NM_000001.1 Homo sapiens fake gene one, mRNA
ACGTACGTAC
GTACGTACGT

>NM_000002.2 Homo sapiens fake gene two, mRNA
acgtNNacgt
"""


class TestParseFasta(unittest.TestCase):
    def records(self, text: str) -> list[tile_mod.Record]:
        return list(tile_mod.parse_fasta(io.StringIO(text), source_file="fake.fasta"))

    def test_parses_both_records(self):
        records = self.records(TWO_RECORD_FASTA)
        self.assertEqual([r.transcript_id for r in records], ["NM_000001.1", "NM_000002.2"])

    def test_concatenates_wrapped_lines_and_uppercases(self):
        first, second = self.records(TWO_RECORD_FASTA)
        self.assertEqual(first.sequence, "ACGTACGTACGTACGTACGT")
        self.assertEqual(second.sequence, "ACGTNNACGT")

    def test_keeps_description_and_source(self):
        first, _ = self.records(TWO_RECORD_FASTA)
        self.assertEqual(first.description, "Homo sapiens fake gene one, mRNA")
        self.assertEqual(first.source_file, "fake.fasta")

    def test_sequence_before_header_is_an_error(self):
        with self.assertRaises(ValueError):
            self.records("ACGT\n>NM_000001.1 x\nACGT\n")

    def test_empty_header_is_an_error(self):
        with self.assertRaises(ValueError):
            self.records(">\nACGT\n")

    def test_empty_input_yields_nothing(self):
        self.assertEqual(self.records(""), [])


class TestTileSequence(unittest.TestCase):
    def test_tile_count_is_len_minus_k_plus_one(self):
        sequence = "A" * 100
        for k in (1, 2, 30, 99, 100):
            with self.subTest(k=k):
                tiles = list(tile_mod.tile_sequence(sequence, k=k))
                self.assertEqual(len(tiles), len(sequence) - k + 1)

    def test_first_and_last_coordinates(self):
        sequence = "ACGT" * 25  # 100 nt
        tiles = list(tile_mod.tile_sequence(sequence, k=30, step=1))
        self.assertEqual((tiles[0].start, tiles[0].end), (1, 30))
        self.assertEqual((tiles[-1].start, tiles[-1].end), (71, 100))

    def test_slicing_round_trip(self):
        sequence = "ACGTTGCA" * 10
        for tile in tile_mod.tile_sequence(sequence, k=30, step=1):
            self.assertEqual(sequence[tile.start - 1 : tile.end], tile.sequence)
            self.assertEqual(tile.length, 30)

    def test_step_advances_by_step(self):
        tiles = list(tile_mod.tile_sequence("A" * 20, k=5, step=3))
        self.assertEqual([t.start for t in tiles], [1, 4, 7, 10, 13, 16])
        self.assertTrue(all(t.end == t.start + 4 for t in tiles))

    def test_sequence_shorter_than_k_yields_nothing(self):
        self.assertEqual(list(tile_mod.tile_sequence("ACGT", k=30)), [])

    def test_sequence_exactly_k_yields_one_tile(self):
        tiles = list(tile_mod.tile_sequence("A" * 30, k=30))
        self.assertEqual(len(tiles), 1)
        self.assertEqual((tiles[0].start, tiles[0].end), (1, 30))

    def test_invalid_k_and_step_rejected(self):
        with self.assertRaises(ValueError):
            list(tile_mod.tile_sequence("ACGT", k=0))
        with self.assertRaises(ValueError):
            list(tile_mod.tile_sequence("ACGT", k=2, step=0))


class TestTileMetadata(unittest.TestCase):
    def test_tile_id_format_and_provenance(self):
        record = tile_mod.Record("NM_006218.4", "desc", "ACGT" * 10, "NM_006218.4.fasta")
        tile = next(tile_mod.tile_record(record, k=30))
        self.assertEqual(tile.tile_id, "NM_006218.4:1-30")
        self.assertEqual(tile.transcript_id, "NM_006218.4")
        self.assertEqual(tile.source_file, "NM_006218.4.fasta")

    def test_has_ambiguous_flags_only_tiles_covering_the_n(self):
        # N at 1-based position 11; with k=10 the tiles spanning it are starts 2..11.
        sequence = "A" * 10 + "N" + "A" * 10
        flagged = [t.start for t in tile_mod.tile_sequence(sequence, k=10) if t.has_ambiguous]
        self.assertEqual(flagged, list(range(2, 12)))

    def test_clean_tiles_are_not_flagged(self):
        self.assertFalse(
            any(t.has_ambiguous for t in tile_mod.tile_sequence("ACGTU" * 10, k=30))
        )


class TestSequenceHelpers(unittest.TestCase):
    def test_biopython_helpers_behave_as_this_module_assumes(self):
        """Guards the Biopython behaviour the pipeline relies on across upgrades.

        These are Biopython's functions, not ours -- the point is to fail loudly
        if an upgrade changes the semantics tiling and guide design depend on.
        """
        # transcribe / back_transcribe are pure T <-> U, preserving case and
        # leaving ambiguity codes alone.
        self.assertEqual(tile_mod.transcribe("ACGT"), "ACGU")
        self.assertEqual(tile_mod.transcribe("acgtn"), "acgun")
        self.assertEqual(tile_mod.back_transcribe("ACGU"), "ACGT")
        # reverse_complement is DNA-alphabet (A -> T) and handles IUPAC codes.
        self.assertEqual(tile_mod.reverse_complement("AAAC"), "GTTT")
        self.assertEqual(tile_mod.reverse_complement("ACGTN"), "NACGT")
        sequence = "ACGTTGCANRY"
        self.assertEqual(
            tile_mod.reverse_complement(tile_mod.reverse_complement(sequence)), sequence
        )
        # UNAMBIGUOUS_BASES comes from IUPACData, so pin what it resolves to.
        self.assertEqual(set(tile_mod.UNAMBIGUOUS_BASES), set("ACGTU"))

    def test_guide_rna_pairs_a_with_u_not_t(self):
        self.assertEqual(tile_mod.guide_rna("AAAA"), "UUUU")
        self.assertEqual(tile_mod.guide_rna("ACGU"), "ACGU")
        self.assertNotIn("T", tile_mod.guide_rna("ACGT"))

    def test_guide_is_written_5_to_3(self):
        # Asymmetric target so a reversed guide cannot pass by accident.
        target = "AAAAACCCCCGGGGGUUUUUAAGGCCUUAC"
        guide = tile_mod.guide_rna(target)
        self.assertEqual(guide, "GUAAGGCCUUAAAAACCCCCGGGGGUUUUU")
        # Guide 5' end pairs with target 3' end: antiparallel, not parallel.
        self.assertTrue(pairs_antiparallel(target, guide))
        self.assertFalse(pairs_antiparallel(target, guide[::-1]))

    def test_guide_5_prime_base_pairs_with_target_3_prime_base(self):
        target = "ACGUACGUAC"
        guide = tile_mod.guide_rna(target)
        self.assertEqual(guide[0], "G")  # pairs with target's final C
        self.assertEqual(target[-1], "C")
        self.assertEqual(guide[-1], "U")  # pairs with target's first A
        self.assertEqual(target[0], "A")

    def test_pairs_antiparallel_rejects_mismatches_and_bad_lengths(self):
        self.assertFalse(pairs_antiparallel("ACGU", "ACGG"))
        self.assertFalse(pairs_antiparallel("ACGU", "ACG"))
        self.assertFalse(pairs_antiparallel("ACGN", tile_mod.guide_rna("ACGN")))

    def test_guide_rna_is_an_involution_on_rna(self):
        target = "ACGUUGCA"
        self.assertEqual(tile_mod.guide_rna(tile_mod.guide_rna(target)), target)

    def test_guide_rna_length_is_preserved(self):
        self.assertEqual(len(tile_mod.guide_rna("ACGT" * 7 + "AC")), 30)


class TestCli(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.input_dir = os.path.join(self.tmp.name, "data")
        os.makedirs(self.input_dir)
        self.output = os.path.join(self.tmp.name, "results", "tiles.tsv")

    def write_fasta(self, name: str, text: str) -> None:
        with open(os.path.join(self.input_dir, name), "w", encoding="utf-8") as handle:
            handle.write(text)

    def run_cli(self, *extra: str) -> list[dict[str, str]]:
        stderr = io.StringIO()
        argv = ["--input-dir", self.input_dir, "--output", self.output, *extra]
        if "--input-type" not in extra:
            argv += ["--input-type", "mrna"]
        old = sys.stderr
        sys.stderr = stderr
        try:
            self.assertEqual(tile_transcripts.main(argv), 0)
        finally:
            sys.stderr = old
        with open(self.output, newline="", encoding="utf-8") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_tiles_every_file_in_the_folder(self):
        self.write_fasta("a.fasta", ">NM_1.1 a\n" + "ACGT" * 10 + "\n")
        self.write_fasta("b.fasta", ">NM_2.1 b\n" + "TGCA" * 10 + "\n")
        rows = self.run_cli("-k", "30")
        self.assertEqual(len(rows), 2 * (40 - 30 + 1))
        self.assertEqual({r["transcript_id"] for r in rows}, {"NM_1.1", "NM_2.1"})
        self.assertEqual({r["source_file"] for r in rows}, {"a.fasta", "b.fasta"})
        self.assertEqual({r["strand"] for r in rows}, {"+"})

    def test_output_is_deterministic_and_coordinates_round_trip(self):
        sequence = "ACGTTGCA" * 8
        self.write_fasta("a.fasta", ">NM_1.1 a\n" + sequence + "\n")
        rows = self.run_cli("-k", "30")
        self.assertEqual([int(r["start"]) for r in rows], list(range(1, len(sequence) - 29 + 1)))
        for row in rows:
            start, end = int(row["start"]), int(row["end"])
            # target_seq is verbatim, so coordinates index the input directly.
            self.assertEqual(sequence[start - 1 : end], row["target_seq"])
            self.assertEqual(row["tile_id"], f"NM_1.1:{start}-{end}")

    def test_guide_is_the_rna_reverse_complement_of_the_target(self):
        self.write_fasta("a.fasta", ">NM_1.1 a\n" + "ACGT" * 10 + "\n")
        rows = self.run_cli("-k", "30")
        for row in rows:
            self.assertEqual(row["guide_rna"], tile_mod.guide_rna(row["target_seq"]))
            # Guide pairs with target: complementing back recovers it, in RNA
            # space (the target itself stays in the input's alphabet).
            self.assertEqual(
                tile_mod.guide_rna(row["guide_rna"]), tile_mod.transcribe(row["target_seq"])
            )

    def test_every_output_guide_is_5_to_3_and_antiparallel(self):
        self.write_fasta("a.fasta", ">NM_1.1 a\n" + "AAGGCCTTACGTTGCAGGTTCCAAGTACGT" * 3 + "\n")
        for extra in (("--guide-alphabet", "rna"), ("--guide-alphabet", "dna")):
            with self.subTest(alphabet=extra[1]):
                rows = self.run_cli("-k", "30", *extra)
                self.assertTrue(rows)
                for row in rows:
                    self.assertTrue(
                        pairs_antiparallel(row["target_seq"], row["guide_rna"]),
                        f"guide not antiparallel 5'->3' for {row['tile_id']}",
                    )

    def test_target_keeps_the_input_alphabet_while_guide_is_rna(self):
        self.write_fasta("a.fasta", ">NM_1.1 a\n" + "ACGT" * 10 + "\n")
        rows = self.run_cli("-k", "30")
        for row in rows:
            # DNA in, so the target keeps its T and gains no U.
            self.assertIn("T", row["target_seq"])
            self.assertNotIn("U", row["target_seq"])
            # The guide is the real molecule: RNA.
            self.assertIn("U", row["guide_rna"])
            self.assertNotIn("T", row["guide_rna"])

    def test_target_is_a_verbatim_substring_of_the_input(self):
        sequence = "ACGTTGCAGGTTCCAAGTACGTAAGGCCTTACGT"
        self.write_fasta("a.fasta", ">NM_1.1 a\n" + sequence + "\n")
        for extra in (("--input-type", "dna"), ("--input-type", "mrna")):
            with self.subTest(input_type=extra[1]):
                for row in self.run_cli("-k", "30", *extra):
                    self.assertIn(row["target_seq"], sequence)

    def test_dna_guide_alphabet_for_aligners_leaves_target_alone(self):
        self.write_fasta("a.fasta", ">NM_1.1 a\n" + "ACGT" * 10 + "\n")
        rows = self.run_cli("-k", "30", "--guide-alphabet", "dna")
        for row in rows:
            self.assertNotIn("U", row["guide_rna"])
            self.assertIn("T", row["target_seq"])

    def test_declared_type_does_not_change_the_output(self):
        # Transcription only affects the guide, which is RNA either way.
        self.write_fasta("a.fasta", ">NM_1.1 a\n" + "ACGT" * 10 + "\n")
        self.assertEqual(
            self.run_cli("-k", "30", "--input-type", "dna"),
            self.run_cli("-k", "30", "--input-type", "mrna"),
        )

    def test_rna_input_declared_as_mrna_is_accepted_unchanged(self):
        sequence = "ACGU" * 10
        self.write_fasta("a.fasta", ">NM_1.1 a\n" + sequence + "\n")
        rows = self.run_cli("-k", "30", "--input-type", "mrna")
        self.assertEqual(rows[0]["target_seq"], sequence[:30])

    def test_fasta_field_can_be_the_target(self):
        self.write_fasta("a.fasta", ">NM_1.1 a\n" + "ACGT" * 10 + "\n")
        fasta_out = os.path.join(self.tmp.name, "results", "targets.fasta")
        rows = self.run_cli("-k", "30", "--fasta-out", fasta_out, "--fasta-field", "target_seq")
        with open(fasta_out, encoding="utf-8") as handle:
            emitted = list(tile_mod.parse_fasta(handle))
        self.assertEqual([r.sequence for r in emitted], [r["target_seq"] for r in rows])

    def test_ambiguous_tiles_kept_and_flagged_by_default(self):
        self.write_fasta("a.fasta", ">NM_1.1 a\n" + "A" * 30 + "N" + "A" * 30 + "\n")
        rows = self.run_cli("-k", "30")
        self.assertEqual(len(rows), 61 - 30 + 1)
        self.assertEqual(sum(int(r["has_ambiguous"]) for r in rows), 30)

    def test_drop_ambiguous_omits_them(self):
        self.write_fasta("a.fasta", ">NM_1.1 a\n" + "A" * 30 + "N" + "A" * 30 + "\n")
        rows = self.run_cli("-k", "30", "--drop-ambiguous")
        self.assertEqual(len(rows), 2)
        self.assertTrue(all(int(r["has_ambiguous"]) == 0 for r in rows))

    def test_fasta_out_headers_are_tile_ids(self):
        self.write_fasta("a.fasta", ">NM_1.1 a\n" + "ACGT" * 10 + "\n")
        fasta_out = os.path.join(self.tmp.name, "results", "tiles.fasta")
        rows = self.run_cli("-k", "30", "--fasta-out", fasta_out)
        with open(fasta_out, encoding="utf-8") as handle:
            emitted = list(tile_mod.parse_fasta(handle))
        self.assertEqual([r.transcript_id for r in emitted], [r["tile_id"] for r in rows])
        self.assertEqual([r.sequence for r in emitted], [r["guide_rna"] for r in rows])

    def test_per_transcript_files(self):
        self.write_fasta("a.fasta", ">NM_1.1 a\n" + "ACGT" * 10 + "\n")
        self.write_fasta("b.fasta", ">NM_2.1 b\n" + "TGCA" * 10 + "\n")
        self.run_cli("-k", "30", "--per-transcript")
        per_dir = os.path.join(os.path.dirname(self.output), "tiles")
        self.assertEqual(sorted(os.listdir(per_dir)), ["NM_1.1.tsv", "NM_2.1.tsv"])
        with open(os.path.join(per_dir, "NM_1.1.tsv"), newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 11)
        self.assertEqual({r["transcript_id"] for r in rows}, {"NM_1.1"})

    def test_sequence_shorter_than_k_is_an_error(self):
        self.write_fasta("a.fasta", ">NM_1.1 a\nACGT\n")
        with self.assertRaises(SystemExit):
            self.run_cli("-k", "30")

    def test_duplicate_transcript_id_is_an_error(self):
        self.write_fasta("a.fasta", ">NM_1.1 a\n" + "ACGT" * 10 + "\n")
        self.write_fasta("b.fasta", ">NM_1.1 a again\n" + "ACGT" * 10 + "\n")
        with self.assertRaises(SystemExit):
            self.run_cli("-k", "30")

    def test_empty_input_dir_is_an_error(self):
        with self.assertRaises(SystemExit):
            self.run_cli("-k", "30")

    def test_missing_input_dir_is_an_error(self):
        self.input_dir = os.path.join(self.tmp.name, "nope")
        with self.assertRaises(SystemExit):
            self.run_cli("-k", "30")


if __name__ == "__main__":
    unittest.main()
