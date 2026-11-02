"""Tests for the QC filter in src/visualise_scatterlines.py.

Scope is `drop_qc_flagged` and its interaction with `transcript_lengths`. The
plotly layout is not covered here -- this file is about the one piece of that
module that *removes data*, where being silently wrong costs guides off a plot
with nothing to show it happened.

The fixtures spell out each outcome as its own guide, so a failure names the
guide rather than a count.
"""

from __future__ import annotations

import os
import sys
import unittest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "src"))
sys.path.insert(0, os.path.join(_ROOT, "scripts"))

try:
    import pandas as pd
except ImportError:  # pragma: no cover - environment.yml installs pandas
    pd = None

import seqqc  # noqa: E402

if pd is not None:
    from visualise_scatterlines import drop_qc_flagged, transcript_lengths


# One guide per outcome. Lengths are 30, the spacer width stage 1 tiles at.
GUIDES = [
    # (tile_id, transcript_id, start, spacer_seq, expected flags)
    ("T1:1-30",   "T1",  1, "GCAUGCAUGCAUGCAUGCAUGCAUGCAUGA", ()),
    ("T1:2-31",   "T1",  2, "GCGCGCGCGCGCAUAUAUAUAUAUAUAUAU", ()),       # 40% exactly
    ("T1:3-32",   "T1",  3, "GCGCGCGCGCGCGCGCGCAUAUAUAUAUAU", ()),       # 60% exactly
    ("T1:4-33",   "T1",  4, "GCGCGCGCGCGAUAUAUAUAUAUAUAUAUA", ("gc_low",)),
    ("T1:5-34",   "T1",  5, "GCGCGCGCGCGCGCGCGCGAUAUAUAUAUA", ("gc_high",)),
    ("T1:6-35",   "T1",  6, "UUUUGCAUGCAUGCAUGCAUGCAUGCAUAU",
     ("homopolymer_t",)),
    ("T1:7-36",   "T1",  7, "GGGGAUAUGCAUGCAUAUAUGCAUGCAUAU",
     ("homopolymer_v",)),
    ("T1:8-37",   "T1",  8, "UUUUGGGGAUGCAUGCAUGCAUGCAUGCAU",
     ("homopolymer_t", "homopolymer_v")),
    ("T2:1-30",   "T2",  1, "AAAAAUAUAUAUAUAUAUAUAUAUAUAUGC",
     ("gc_low", "homopolymer_v")),
    ("T2:100-129", "T2", 100, "GCAUGCAUGCAUGCAUGCAUGCAUGCAUGA", ()),
]


def frame(guides=GUIDES):
    return pd.DataFrame(
        {
            "tile_id": [g[0] for g in guides],
            "transcript_id": [g[1] for g in guides],
            "start": [g[2] for g in guides],
            "end": [g[2] + 29 for g in guides],
            "spacer_seq": [g[3] for g in guides],
            "predicted_efficacy": [0.5] * len(guides),
        }
    )


@unittest.skipIf(pd is None, "pandas not installed")
class TestDropQcFlagged(unittest.TestCase):
    def setUp(self):
        self.df = frame()
        self.qc = drop_qc_flagged(self.df)

    def test_the_fixture_matches_what_seqqc_says_about_it(self):
        # Pins the fixture itself. If these expectations drift from
        # seqqc.check_guide, every assertion below is testing the wrong thing.
        for tile_id, _, _, spacer, expected in GUIDES:
            with self.subTest(tile_id=tile_id):
                self.assertEqual(
                    seqqc.check_guide(tile_id, spacer).flags, expected
                )

    def test_clean_guides_are_kept(self):
        kept = set(self.qc.kept["tile_id"])
        for tile_id, _, _, _, expected in GUIDES:
            if not expected:
                with self.subTest(tile_id=tile_id):
                    self.assertIn(tile_id, kept)

    def test_flagged_guides_are_dropped(self):
        kept = set(self.qc.kept["tile_id"])
        for tile_id, _, _, _, expected in GUIDES:
            if expected:
                with self.subTest(tile_id=tile_id):
                    self.assertNotIn(tile_id, kept)

    def test_both_band_edges_survive(self):
        # Inclusive edges, as seqqc defines them: a guide at exactly 40% or 60%
        # GC must not be dropped off the plot.
        kept = set(self.qc.kept["tile_id"])
        self.assertIn("T1:2-31", kept)
        self.assertIn("T1:3-32", kept)

    def test_counts_add_up(self):
        self.assertEqual(self.qc.n_input, len(GUIDES))
        self.assertEqual(self.qc.n_kept + self.qc.n_dropped, self.qc.n_input)
        self.assertEqual(self.qc.n_kept, 4)
        self.assertEqual(self.qc.n_dropped, 6)

    def test_per_flag_counts_overlap_rather_than_summing(self):
        # T2:1-30 raises two flags, so the per-flag counts sum to more than
        # the number of guides dropped. A summary that added them would
        # overstate the damage.
        self.assertEqual(
            self.qc.by_flag,
            {"gc_low": 2, "gc_high": 1, "homopolymer_t": 2, "homopolymer_v": 3},
        )
        self.assertGreater(sum(self.qc.by_flag.values()), self.qc.n_dropped)

    def test_t_and_v_runs_are_counted_under_separate_flags(self):
        # A filter that merged them would report one count of 5 here.
        self.assertEqual(self.qc.by_flag["homopolymer_t"], 2)
        self.assertEqual(self.qc.by_flag["homopolymer_v"], 3)

    def test_a_flag_nobody_raised_is_absent_not_zero(self):
        qc = drop_qc_flagged(frame(GUIDES[:1]))
        self.assertEqual(qc.by_flag, {})
        # So the CLI must read it with .get, which is what it does.
        self.assertEqual(qc.by_flag.get("gc_low", 0), 0)

    def test_dropped_per_transcript(self):
        self.assertEqual(self.qc.dropped_per_transcript, {"T1": 5, "T2": 1})

    def test_a_transcript_losing_nothing_is_absent_from_the_mapping(self):
        qc = drop_qc_flagged(frame([g for g in GUIDES if not g[4]]))
        self.assertEqual(qc.dropped_per_transcript, {})

    def test_row_order_and_index_are_preserved(self):
        kept = self.qc.kept
        self.assertEqual(list(kept["tile_id"]),
                         [g[0] for g in GUIDES if not g[4]])
        self.assertEqual(list(kept.index), sorted(kept.index))

    def test_columns_are_untouched(self):
        self.assertEqual(list(self.qc.kept.columns), list(self.df.columns))

    def test_the_input_frame_is_not_mutated(self):
        drop_qc_flagged(self.df)
        self.assertEqual(len(self.df), len(GUIDES))

    def test_thresholds_are_recorded(self):
        self.assertEqual(self.qc.gc_min, seqqc.DEFAULT_GUIDE_GC_MIN)
        self.assertEqual(self.qc.gc_max, seqqc.DEFAULT_GUIDE_GC_MAX)
        self.assertEqual(self.qc.homopolymer_t_min,
                         seqqc.DEFAULT_GUIDE_HOMOPOLYMER_T_MIN)
        self.assertEqual(self.qc.homopolymer_v_min,
                         seqqc.DEFAULT_GUIDE_HOMOPOLYMER_V_MIN)

    def test_thresholds_are_honoured_not_just_recorded(self):
        loose = drop_qc_flagged(self.df, gc_min=0.0, gc_max=1.0,
                                homopolymer_t_min=99, homopolymer_v_min=99)
        self.assertEqual(loose.n_dropped, 0)
        self.assertEqual(loose.by_flag, {})

        strict = drop_qc_flagged(self.df, gc_min=0.99, gc_max=1.0)
        self.assertEqual(strict.n_kept, 0)

    def test_the_two_run_thresholds_are_independent(self):
        # Relaxing T must leave the V drops standing, and vice versa.
        t_loose = drop_qc_flagged(self.df, homopolymer_t_min=99)
        self.assertNotIn("homopolymer_t", t_loose.by_flag)
        self.assertEqual(t_loose.by_flag["homopolymer_v"], 3)

        v_loose = drop_qc_flagged(self.df, homopolymer_v_min=99)
        self.assertNotIn("homopolymer_v", v_loose.by_flag)
        self.assertEqual(v_loose.by_flag["homopolymer_t"], 2)

    def test_an_everything_dropped_filter_returns_an_empty_frame(self):
        # Not an exception: the CLI is what decides that plotting nothing is
        # an error, because that is a presentation decision.
        qc = drop_qc_flagged(self.df, gc_min=0.99, gc_max=1.0)
        self.assertTrue(qc.kept.empty)
        self.assertEqual(list(qc.kept.columns), list(self.df.columns))

    def test_the_filter_agrees_with_check_guide_row_for_row(self):
        # The contract: this filter must never disagree with the module that
        # produces guide_qc.tsv, or a plot and its QC table describe different
        # guide sets.
        kept = set(self.qc.kept["tile_id"])
        for tile_id, _, _, spacer, _ in GUIDES:
            with self.subTest(tile_id=tile_id):
                self.assertEqual(
                    tile_id in kept, seqqc.check_guide(tile_id, spacer).passes
                )

    def test_rna_and_dna_spacers_give_the_same_verdict(self):
        rna = frame()
        dna = frame()
        dna["spacer_seq"] = dna["spacer_seq"].str.replace("U", "T")
        self.assertEqual(
            list(drop_qc_flagged(rna).kept["tile_id"]),
            list(drop_qc_flagged(dna).kept["tile_id"]),
        )


@unittest.skipIf(pd is None, "pandas not installed")
class TestFilteringDoesNotShrinkTheAxis(unittest.TestCase):
    """The trap the CLI has to avoid, pinned here.

    ``transcript_lengths`` reads the largest tile end. Taking it after the QC
    filter would stop each axis at the last surviving guide, so a QC-poor 3'
    end would read as *absent* rather than as filtered out -- the opposite of
    what drawing the full transcript is for.
    """

    def test_length_from_the_unfiltered_frame_is_the_full_transcript(self):
        df = frame()
        self.assertEqual(transcript_lengths(df), {"T1": 37, "T2": 129})

    def test_the_last_guide_on_a_transcript_being_dropped_would_shrink_it(self):
        # T2's surviving guide starts at 100; its dropped one at 1. Reverse the
        # situation so the *high* coordinate is the one dropped, and show the
        # two lengths differ -- which is why order matters in the CLI.
        guides = [
            ("T3:1-30", "T3", 1, "GCAUGCAUGCAUGCAUGCAUGCAUGCAUGA", ()),
            ("T3:500-529", "T3", 500, "AAAAAUAUAUAUAUAUAUAUAUAUAUAUGC",
             ("gc_low", "homopolymer_v")),
        ]
        df = frame(guides)
        unfiltered = transcript_lengths(df)
        filtered = transcript_lengths(drop_qc_flagged(df).kept)
        self.assertEqual(unfiltered["T3"], 529)
        self.assertEqual(filtered["T3"], 30)
        self.assertNotEqual(unfiltered, filtered)


if __name__ == "__main__":
    unittest.main()
