"""Tests for src/seqqc.py and scripts/check_target_sites.py.

Every scanner is checked against a hand-built sequence with the feature planted
at a known offset, so a failure names the position rather than a count. The
coordinate convention is pinned separately: offsets out of these functions are
1-based and relative to the region, and the CLI is what shifts them onto the
transcript.
"""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import unittest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "src"))
sys.path.insert(0, os.path.join(_ROOT, "scripts"))

import seqqc  # noqa: E402

CHECK_SCRIPT = os.path.join(_ROOT, "scripts", "check_target_sites.py")


class TestGc(unittest.TestCase):
    def test_gc_fraction_counts_only_unambiguous_bases(self):
        gc, ambiguous = seqqc.gc_fraction("GGCCAATT")
        self.assertAlmostEqual(gc, 0.5)
        self.assertEqual(ambiguous, 0)

    def test_ambiguous_bases_leave_the_denominator(self):
        # 2 of 4 real bases are GC; the Ns must not count as not-GC, which
        # would give 0.25 and understate the GC of the real sequence.
        gc, ambiguous = seqqc.gc_fraction("GCATNNNN")
        self.assertAlmostEqual(gc, 0.5)
        self.assertEqual(ambiguous, 4)

    def test_rna_alphabet_is_not_all_ambiguous(self):
        gc, ambiguous = seqqc.gc_fraction("GGCCAAUU")
        self.assertAlmostEqual(gc, 0.5)
        self.assertEqual(ambiguous, 0)

    def test_empty_sequence_does_not_divide_by_zero(self):
        self.assertEqual(seqqc.gc_fraction(""), (0.0, 0))

    def test_sliding_gc_window_positions_are_one_based_inclusive(self):
        windows = seqqc.sliding_gc("ATATGCGC", window=4)
        self.assertEqual(windows[0].start, 1)
        self.assertEqual(windows[0].end, 4)
        self.assertAlmostEqual(windows[0].value, 0.0)
        self.assertEqual(windows[-1].end, 8)
        self.assertAlmostEqual(windows[-1].value, 1.0)

    def test_step_is_one_so_every_offset_is_measured(self):
        self.assertEqual(len(seqqc.sliding_gc("A" * 100, window=30)), 71)

    def test_short_sequence_yields_one_window_rather_than_none(self):
        windows = seqqc.sliding_gc("GGC", window=30)
        self.assertEqual(len(windows), 1)
        self.assertEqual((windows[0].start, windows[0].end), (1, 3))


class TestHomopolymers(unittest.TestCase):
    def test_run_is_found_at_the_right_offset(self):
        # TTTT at 1-based offsets 5-8.
        runs = seqqc.homopolymer_runs("ACGC" + "TTTT" + "ACGC", 4, bases="T")
        self.assertEqual(len(runs), 1)
        self.assertEqual((runs[0].base, runs[0].start, runs[0].end), ("T", 5, 8))
        self.assertEqual(runs[0].length, 4)

    def test_threshold_is_inclusive(self):
        self.assertEqual(len(seqqc.homopolymer_runs("ACTTTTCA", 4, bases="T")), 1)
        self.assertEqual(len(seqqc.homopolymer_runs("ACTTTCA", 4, bases="T")), 0)

    def test_run_at_the_very_start_and_end_is_found(self):
        runs = seqqc.homopolymer_runs("TTTTACGTAAAA", 4)
        self.assertEqual([(r.base, r.start, r.end) for r in runs],
                         [("T", 1, 4), ("A", 9, 12)])

    def test_bases_filter_excludes_other_homopolymers(self):
        runs = seqqc.homopolymer_runs("AAAAATTTTT", 5, bases="ACG")
        self.assertEqual([r.base for r in runs], ["A"])

    def test_u_is_treated_as_t(self):
        runs = seqqc.homopolymer_runs("ACGCUUUUACGC", 4, bases="T")
        self.assertEqual((runs[0].start, runs[0].end), (5, 8))


class TestTandemRepeats(unittest.TestCase):
    def test_dinucleotide_repeat_is_found_with_its_copy_count(self):
        # (AT) x4 at offsets 4-11.
        repeats = seqqc.tandem_repeats("GGC" + "ATATATAT" + "GGC", min_copies=3)
        self.assertEqual(len(repeats), 1)
        self.assertEqual(
            (repeats[0].motif, repeats[0].start, repeats[0].end, repeats[0].copies),
            ("AT", 4, 11, 4),
        )

    def test_trinucleotide_repeat_is_found(self):
        # The motif is one *rotation* of the repeat -- whichever starts at the
        # earliest qualifying offset -- so this asserts the span and the copy
        # count, and only that the motif is a rotation of CAG. Pinning the
        # rotation would be pinning an artefact of the flanking sequence.
        repeats = seqqc.tandem_repeats("GG" + "CAGCAGCAG" + "GG", min_copies=3)
        self.assertEqual(len(repeats), 1)
        self.assertEqual(repeats[0].copies, 3)
        self.assertEqual(repeats[0].length, 9)
        self.assertIn(repeats[0].motif, {"CAG", "AGC", "GCA"})

    def test_homopolymer_is_not_also_reported_as_a_repeat(self):
        # (AA) x3 is a 6 nt poly-A; homopolymer_runs covers it, and reporting
        # it here too would double-count one feature.
        self.assertEqual(seqqc.tandem_repeats("GC" + "AAAAAA" + "GC"), [])

    def test_below_min_copies_is_not_reported(self):
        self.assertEqual(seqqc.tandem_repeats("GC" + "ATAT" + "GC", min_copies=3), [])

    def test_one_row_per_feature_not_one_per_nested_motif(self):
        # (AT) x6 is also (ATAT) x3 and contains shorter (AT) runs. One feature.
        repeats = seqqc.tandem_repeats("GG" + "AT" * 6 + "GG", min_copies=3)
        self.assertEqual(len(repeats), 1)
        self.assertEqual((repeats[0].start, repeats[0].end), (3, 14))

    def test_two_separate_repeats_are_both_reported(self):
        sequence = "GG" + "ATATATAT" + "CCCC" + "CAGCAGCAG" + "GG"
        motifs = {r.motif for r in seqqc.tandem_repeats(sequence, min_copies=3)}
        self.assertEqual(motifs, {"AT", "CAG"})


class TestEntropy(unittest.TestCase):
    def test_homopolymer_is_zero_bits(self):
        self.assertAlmostEqual(seqqc.shannon_entropy("AAAAAAAA"), 0.0)

    def test_uniform_four_bases_is_two_bits(self):
        self.assertAlmostEqual(seqqc.shannon_entropy("ACGTACGT"), 2.0)

    def test_two_bases_evenly_is_one_bit(self):
        self.assertAlmostEqual(seqqc.shannon_entropy("ATATATAT"), 1.0)

    def test_sliding_entropy_finds_the_low_complexity_window(self):
        sequence = "ACGT" * 10 + "A" * 30 + "ACGT" * 10
        windows = seqqc.sliding_entropy(sequence, window=30)
        lowest = min(windows, key=lambda w: w.value)
        self.assertAlmostEqual(lowest.value, 0.0)
        # The poly-A runs 41-70, so the all-A window starts there.
        self.assertEqual(lowest.start, 41)


class TestStopCodons(unittest.TestCase):
    def test_each_codon_is_found_at_its_offset(self):
        self.assertEqual(
            [(s.codon, s.start) for s in seqqc.stop_codons("GGTAAGG")],
            [("TAA", 3)],
        )

    def test_all_three_codons_are_recognised(self):
        found = {s.codon for s in seqqc.stop_codons("TAAcccTAGcccTGA")}
        self.assertEqual(found, {"TAA", "TAG", "TGA"})

    def test_stop_codons_can_never_overlap(self):
        # Not an assumption: none of TAA/TAG/TGA holds a T except as its first
        # base, so no codon's 2nd or 3rd base can begin another. Proven here by
        # exhaustion rather than asserted, so it fails if the codon set changes.
        for first in seqqc.STOP_CODONS:
            for second in seqqc.STOP_CODONS:
                for shift in (1, 2):
                    self.assertNotEqual(
                        first[shift:], second[:3 - shift],
                        f"{first} and {second} would overlap at shift {shift}",
                    )

    def test_adjacent_stops_are_both_found(self):
        self.assertEqual(
            [(s.codon, s.start) for s in seqqc.stop_codons("TAATGA")],
            [("TAA", 1), ("TGA", 4)],
        )

    def test_frame_is_irrelevant_so_offset_two_counts(self):
        self.assertEqual(len(seqqc.stop_codons("GTAAG")), 1)

    def test_end_coordinate_spans_three_bases(self):
        stop = seqqc.stop_codons("GGTAA")[0]
        self.assertEqual((stop.start, stop.end), (3, 5))

    def test_rna_alphabet_is_scanned(self):
        self.assertEqual(len(seqqc.stop_codons("GGUAAGG")), 1)

    def test_none_in_a_stopless_sequence(self):
        self.assertEqual(seqqc.stop_codons("GGCCGGCCGGCC"), [])


class TestWidestRun(unittest.TestCase):
    def test_counts_the_longest_consecutive_stretch(self):
        self.assertEqual(seqqc.widest_run([True, False, True, True, True, False]), 3)

    def test_scattered_and_consecutive_are_distinguished(self):
        scattered = [True, False, True, False, True, False]
        consecutive = [True, True, True, False, False, False]
        self.assertEqual(sum(scattered), sum(consecutive))
        self.assertNotEqual(
            seqqc.widest_run(scattered), seqqc.widest_run(consecutive)
        )

    def test_empty_and_all_false(self):
        self.assertEqual(seqqc.widest_run([]), 0)
        self.assertEqual(seqqc.widest_run([False, False]), 0)


class TestCli(unittest.TestCase):
    """End to end, including the coordinate shift onto the transcript."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.data = os.path.join(cls.tmp.name, "data")
        cls.out = os.path.join(cls.tmp.name, "out")
        os.makedirs(cls.data)

        # 100 nt of filler, then a region holding one planted feature of each
        # kind, then filler. The region is 1-based 101-160.
        filler = "ACGTGCATGC" * 10
        region = (
            "ACGTACGTAC"      # 101-110, ordinary
            "TTTTT"           # 111-115, poly-T, 5
            "ACGTA"           # 116-120
            "GGGGGG"          # 121-126, homopolymer G, 6
            "ATATATAT"        # 127-134, (AT) x4
            "CGTAACGTAC"      # 135-144, holds TAA at 138
            "GCGCGCGCGCGC"    # 145-156, GC-rich
            "ACGT"            # 157-160
        )
        cls.region_start, cls.region_end = 101, 100 + len(region)
        with open(os.path.join(cls.data, "NM_999999.1.fasta"), "w") as handle:
            handle.write(
                ">NM_999999.1 Homo sapiens test transcript (TESTG), mRNA\n"
                f"{filler}{region}{filler}\n"
            )

        cls.selection = os.path.join(cls.tmp.name, "selection.tsv")
        with open(cls.selection, "w") as handle:
            handle.write("transcript\tgene\tstart\tend\tlength\n")
            handle.write(
                f"NM_999999.1\tTESTG\t{cls.region_start}\t{cls.region_end}"
                f"\t{len(region)}\n"
            )

        cls.completed = subprocess.run(
            [
                sys.executable, CHECK_SCRIPT,
                "--selection", cls.selection,
                "--input-dir", cls.data,
                "--output-dir", cls.out,
            ],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        assert cls.completed.returncode == 0, cls.completed.stderr.decode()
        cls.summary = cls._read(os.path.join(cls.out, "target_site_qc.tsv"))
        cls.features = cls._read(
            os.path.join(cls.out, "target_site_qc_features.tsv")
        )

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    @staticmethod
    def _read(path):
        with open(path, newline="", encoding="utf-8") as handle:
            lines = [line.rstrip("\n").split("\t") for line in handle]
        header = lines[0]
        return [dict(zip(header, row)) for row in lines[1:]]

    def test_writes_all_three_files(self):
        for name in ("target_site_qc.tsv", "target_site_qc_features.tsv",
                     "target_site_qc.md"):
            self.assertTrue(os.path.isfile(os.path.join(self.out, name)), name)

    def test_one_summary_row_per_region(self):
        self.assertEqual(len(self.summary), 1)
        self.assertEqual(self.summary[0]["gene"], "TESTG")

    def test_declared_length_agreeing_is_reported_as_matching(self):
        row = self.summary[0]
        self.assertEqual(row["length_actual"], row["length_declared"])
        self.assertEqual(row["length_matches"], "1")

    def test_feature_coordinates_are_on_the_transcript_not_the_region(self):
        polyt = [f for f in self.features if f["feature_type"] == "polyT"]
        self.assertEqual(len(polyt), 1)
        # Planted at region offset 11, i.e. transcript 111.
        self.assertEqual(polyt[0]["start"], "111")
        self.assertEqual(polyt[0]["offset_in_region"], "11")
        self.assertEqual(polyt[0]["length"], "5")

    def test_offset_and_transcript_coordinate_are_consistent(self):
        for feature in self.features:
            self.assertEqual(
                int(feature["start"]) - int(feature["offset_in_region"]) + 1,
                self.region_start,
                feature,
            )

    def test_homopolymer_and_repeat_are_found_where_planted(self):
        homo = [f for f in self.features if f["feature_type"] == "homopolymer"]
        self.assertEqual([(f["start"], f["length"]) for f in homo], [("121", "6")])
        # Two repeats, not one: the planted (AT) x4, and the GC-rich stretch,
        # which genuinely *is* a dinucleotide tandem repeat and must be
        # reported as one. It starts a base earlier than the GC run, because
        # the preceding base completes a further CG copy.
        repeats = [f for f in self.features if f["feature_type"] == "tandem_repeat"]
        self.assertEqual([(f["start"], f["detail"]) for f in repeats],
                         [("127", "(AT) x4"), ("144", "(CG) x6")])

    def test_stop_codons_are_info_not_warnings(self):
        stops = [f for f in self.features if f["feature_type"] == "stop_codon"]
        self.assertTrue(stops)
        self.assertEqual({f["severity"] for f in stops}, {"info"})

    def test_stop_codons_do_not_raise_a_flag(self):
        # The whole point of reporting them as counts: presence is not a defect.
        self.assertNotIn("stop", self.summary[0]["flags"])

    def test_poly_t_is_not_also_counted_as_a_homopolymer(self):
        homo = [f for f in self.features if f["feature_type"] == "homopolymer"]
        self.assertEqual([f["detail"][0] for f in homo], ["G"])

    def test_stop_counts_add_up(self):
        row = self.summary[0]
        self.assertEqual(
            int(row["n_stops_total"]),
            int(row["n_stop_TAA"]) + int(row["n_stop_TAG"]) + int(row["n_stop_TGA"]),
        )

    def test_thresholds_are_echoed_into_the_row(self):
        row = self.summary[0]
        self.assertEqual(row["gc_band"], "0.3-0.7")
        self.assertEqual(row["polyT_min"], "4")
        self.assertEqual(row["homopolymer_min"], "5")


class TestCliRejectsBadInput(unittest.TestCase):
    """The failure modes this table has actually shown."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.data = os.path.join(self.tmp.name, "data")
        os.makedirs(self.data)
        with open(os.path.join(self.data, "NM_999999.1.fasta"), "w") as handle:
            handle.write(">NM_999999.1 test (TESTG), mRNA\n" + "ACGT" * 100 + "\n")

    def tearDown(self):
        self.tmp.cleanup()

    def _run(self, text):
        path = os.path.join(self.tmp.name, "selection.tsv")
        with open(path, "w") as handle:
            handle.write(text)
        return subprocess.run(
            [
                sys.executable, CHECK_SCRIPT,
                "--selection", path,
                "--input-dir", self.data,
                "--output-dir", os.path.join(self.tmp.name, "out"),
            ],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )

    def test_header_with_spaces_instead_of_a_tab_is_rejected(self):
        # The real failure: "start    end" as one field, four header fields
        # over five-field rows, which DictReader reads silently and wrongly.
        completed = self._run(
            "transcript\tgene\tstart    end\tlength\n"
            "NM_999999.1\tTESTG\t1\t100\t100\n"
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn(b"tabs", completed.stderr)

    def test_length_mismatch_warns_but_still_reports(self):
        # A table disagreeing with itself is worth a loud warning, not a
        # refusal to describe the region: the span is what the coordinates say.
        completed = self._run(
            "transcript\tgene\tstart\tend\tlength\n"
            "NM_999999.1\tTESTG\t1\t100\t99\n"
        )
        self.assertEqual(completed.returncode, 0, completed.stderr.decode())
        self.assertIn(b"disagrees with itself", completed.stderr)
        rows = TestCli._read(
            os.path.join(self.tmp.name, "out", "target_site_qc.tsv")
        )
        self.assertEqual(rows[0]["length_actual"], "100")
        self.assertEqual(rows[0]["length_declared"], "99")
        self.assertEqual(rows[0]["length_matches"], "0")

    def test_region_running_past_the_transcript_is_rejected(self):
        completed = self._run(
            "transcript\tgene\tstart\tend\tlength\n"
            "NM_999999.1\tTESTG\t1\t5000\t5000\n"
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn(b"only", completed.stderr)

    def test_missing_fasta_is_rejected_naming_the_gene(self):
        completed = self._run(
            "transcript\tgene\tstart\tend\tlength\n"
            "NM_000000.0\tOTHER\t1\t100\t100\n"
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn(b"OTHER", completed.stderr)


if __name__ == "__main__":
    unittest.main()
