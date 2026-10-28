"""Tests for src/hotspots.py and scripts/find_hotspot_regions.py.

    python -m unittest discover -s test -v

Everything here is synthetic and hand-checkable: the fixtures plant hotspots at
known coordinates so the right answer is known by construction rather than by
re-running the code under test. Nothing needs an aligner or the committed
results.
"""

from __future__ import annotations

import io
import os
import sys
import tempfile
import unittest
from contextlib import redirect_stderr

import pandas as pd

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "src"))
sys.path.insert(0, os.path.join(_ROOT, "scripts"))

import find_hotspot_regions  # noqa: E402
import hotspots as hs  # noqa: E402

K = 30


def tiles(transcript: str, length: int, scores) -> pd.DataFrame:
    """Step-1 tiling of `length` bases, with `scores(start) -> efficacy`.

    Mirrors stage 1: one 30-mer starting at every position that fits, 1-based
    inclusive coordinates.
    """
    rows = []
    for start in range(1, length - K + 2):
        end = start + K - 1
        rows.append(
            {
                "tile_id": f"{transcript}:{start}-{end}",
                "transcript_id": transcript,
                "start": start,
                "end": end,
                "predicted_efficacy": float(scores(start)),
            }
        )
    return pd.DataFrame(rows)


def plateau(spans, high=0.95, low=0.10):
    """Scorer giving `high` to tiles starting inside any span, else `low`."""

    def score(start: int) -> float:
        return high if any(a <= start <= b for a, b in spans) else low

    return score


class TestMergeIntervals(unittest.TestCase):
    def test_overlapping_spans_merge(self):
        self.assertEqual(
            hs.merge_intervals([(1, 30), (15, 44), (100, 129)]),
            [(1, 44), (100, 129)],
        )

    def test_abutting_spans_merge_at_gap_zero(self):
        self.assertEqual(hs.merge_intervals([(1, 30), (31, 60)]), [(1, 60)])

    def test_a_one_base_gap_does_not_merge_at_gap_zero(self):
        self.assertEqual(
            hs.merge_intervals([(1, 30), (32, 61)]), [(1, 30), (32, 61)]
        )

    def test_gap_counts_the_bases_strictly_between(self):
        # 31 is the only base between them, so gap=1 merges and gap=0 does not.
        self.assertEqual(hs.merge_intervals([(1, 30), (32, 61)], gap=1), [(1, 61)])

    def test_negative_gap_requires_genuine_overlap(self):
        self.assertEqual(
            hs.merge_intervals([(1, 30), (31, 60)], gap=-1), [(1, 30), (31, 60)]
        )

    def test_unsorted_input_is_handled(self):
        self.assertEqual(
            hs.merge_intervals([(100, 129), (1, 30), (15, 44)]),
            [(1, 44), (100, 129)],
        )

    def test_a_span_contained_in_another_does_not_shrink_it(self):
        self.assertEqual(hs.merge_intervals([(1, 100), (10, 20)]), [(1, 100)])

    def test_no_intervals(self):
        self.assertEqual(hs.merge_intervals([]), [])


class TestFindHotspots(unittest.TestCase):
    """One planted plateau must give exactly one hotspot, not one per tile."""

    def test_a_single_plateau_is_one_hotspot_not_many_guides(self):
        # Tiles starting in 101..150 score high: 50 overlapping guides, one site.
        df = tiles("NM_1", 400, plateau([(101, 150)]))
        found = hs.find_hotspots(df, min_efficacy=0.9)
        self.assertEqual(len(found), 1)
        spot = found[0]
        self.assertEqual(spot.n_guides, 50)
        self.assertEqual((spot.start, spot.end), (101, 150 + K - 1))
        self.assertEqual(spot.length, 79)

    def test_two_plateaus_far_apart_are_two_hotspots(self):
        df = tiles("NM_1", 600, plateau([(101, 130), (401, 430)]))
        found = hs.find_hotspots(df, min_efficacy=0.9)
        self.assertEqual([(h.start, h.end) for h in found], [(101, 159), (401, 459)])

    def test_plateaus_whose_tiles_overlap_merge_into_one(self):
        # Tiles at 101 and 120 overlap (both 30 nt), so this is one site.
        df = tiles("NM_1", 400, plateau([(101, 105), (120, 125)]))
        self.assertEqual(len(hs.find_hotspots(df, min_efficacy=0.9)), 1)

    def test_merge_gap_controls_how_many_sites_are_reported(self):
        # Tiles 101-110 end at 139; the next run starts at 171. 31 bases
        # between, so gap=0 keeps them apart and gap=31 fuses them.
        df = tiles("NM_1", 400, plateau([(101, 110), (171, 180)]))
        self.assertEqual(len(hs.find_hotspots(df, min_efficacy=0.9, merge_gap=0)), 2)
        self.assertEqual(len(hs.find_hotspots(df, min_efficacy=0.9, merge_gap=31)), 1)

    def test_threshold_is_inclusive(self):
        df = tiles("NM_1", 100, lambda s: 0.9)
        self.assertEqual(len(hs.find_hotspots(df, min_efficacy=0.9)), 1)
        self.assertEqual(hs.find_hotspots(df, min_efficacy=0.9000001), [])

    def test_below_threshold_guides_do_not_extend_a_hotspot(self):
        df = tiles("NM_1", 400, plateau([(101, 150)]))
        spot = hs.find_hotspots(df, min_efficacy=0.9)[0]
        self.assertEqual(spot.n_guides, 50)
        self.assertEqual(spot.mean_efficacy, 0.95)

    def test_transcripts_are_kept_separate(self):
        df = pd.concat(
            [
                tiles("NM_1", 400, plateau([(101, 130)])),
                tiles("NM_2", 400, plateau([(101, 130)])),
            ],
            ignore_index=True,
        )
        found = hs.find_hotspots(df, min_efficacy=0.9)
        self.assertEqual([h.transcript for h in found], ["NM_1", "NM_2"])

    def test_best_guide_names_the_highest_scoring_tile(self):
        def score(start: int) -> float:
            if start == 110:
                return 0.99
            return 0.95 if 101 <= start <= 130 else 0.1

        spot = hs.find_hotspots(tiles("NM_1", 400, score), min_efficacy=0.9)[0]
        self.assertEqual(spot.best_guide, "NM_1:110-139")
        self.assertAlmostEqual(spot.max_efficacy, 0.99)

    def test_output_is_sorted_by_transcript_then_start(self):
        df = pd.concat(
            [
                tiles("NM_2", 600, plateau([(401, 430), (101, 130)])),
                tiles("NM_1", 600, plateau([(201, 230)])),
            ],
            ignore_index=True,
        )
        found = hs.find_hotspots(df, min_efficacy=0.9)
        self.assertEqual(
            [(h.transcript, h.start) for h in found],
            [("NM_1", 201), ("NM_2", 101), ("NM_2", 401)],
        )


class TestWindowBounds(unittest.TestCase):
    def test_windows_stay_inside_the_region(self):
        region = hs.Region("NM_1", 1000, 4000)
        for start, end, _ in hs.window_bounds(region):
            self.assertGreaterEqual(start, region.start)
            self.assertLessEqual(end, region.end)

    def test_lengths_are_swept_and_the_largest_fitting_one_is_included(self):
        # 1500..2000 step 250 gives 1500, 1750, 2000; the region fits all three.
        region = hs.Region("NM_1", 1, 3000)
        lengths = {e - s + 1 for s, e, _ in hs.window_bounds(region)}
        self.assertEqual(lengths, {1500, 1750, 2000})

    def test_the_largest_fitting_length_is_included_even_off_the_step(self):
        # 1900 is not 1500 + 250k, but it is the widest window that fits.
        region = hs.Region("NM_1", 1, 1900)
        lengths = {e - s + 1 for s, e, _ in hs.window_bounds(region)}
        self.assertIn(1900, lengths)
        self.assertEqual(max(lengths), 1900)

    def test_the_region_tail_is_always_scanned(self):
        # 2600 - 1500 + 1 = 1101 start positions; step 50 lands on 1051, not
        # 1101, so the flush-right window must be added explicitly.
        region = hs.Region("NM_1", 1, 2600)
        bounds = [(s, e) for s, e, _ in hs.window_bounds(region, max_window=1500)]
        self.assertIn((1101, 2600), bounds)

    def test_a_short_region_yields_one_undersized_window_rather_than_none(self):
        region = hs.Region("NM_1", 100, 900)  # 801 bp, below --min-window
        self.assertEqual(
            list(hs.window_bounds(region)), [(100, 900, True)]
        )

    def test_a_region_exactly_min_window_long_is_not_undersized(self):
        region = hs.Region("NM_1", 1, 1500)
        self.assertEqual(list(hs.window_bounds(region)), [(1, 1500, False)])


class TestScoreWindow(unittest.TestCase):
    def setUp(self):
        self.region = hs.Region("NM_1", 1, 5000)
        self.spots = [
            hs.Hotspot("NM_1", 1000, 1099, 10, 0.95, 0.93, "a"),
            hs.Hotspot("NM_1", 2000, 2099, 20, 0.99, 0.95, "b"),
            hs.Hotspot("NM_1", 4000, 4099, 30, 0.91, 0.90, "c"),
        ]

    def test_only_wholly_contained_hotspots_count(self):
        # 1050 clips the first hotspot in half, so it must not count.
        w = hs.score_window("NM_1", self.region, 1050, 3000, self.spots)
        self.assertEqual(w.n_hotspots, 1)
        self.assertEqual(w.n_guides, 20)

    def test_a_hotspot_flush_with_the_window_edge_counts(self):
        w = hs.score_window("NM_1", self.region, 1000, 2099, self.spots)
        self.assertEqual(w.n_hotspots, 2)

    def test_a_hotspot_one_base_past_the_edge_does_not(self):
        w = hs.score_window("NM_1", self.region, 1000, 2098, self.spots)
        self.assertEqual(w.n_hotspots, 1)

    def test_coverage_and_largest_gap(self):
        w = hs.score_window("NM_1", self.region, 1000, 2099, self.spots)
        self.assertEqual(w.covered_bp, 200)
        self.assertEqual(w.length, 1100)
        # 1100..1999 is the empty stretch between the two hotspots, and both
        # ends are flush, so it is the largest gap.
        self.assertEqual(w.largest_gap, 900)

    def test_largest_gap_includes_the_window_edges(self):
        # One hotspot at 2000-2099 inside 1-3000: the tail 2100..3000 is 901.
        w = hs.score_window("NM_1", self.region, 1, 3000, self.spots[1:2])
        self.assertEqual(w.largest_gap, 1999)  # the head, 1..1999

    def test_an_empty_window_scores_zero_rather_than_failing(self):
        w = hs.score_window("NM_1", self.region, 3000, 3500, self.spots)
        self.assertEqual(w.n_hotspots, 0)
        self.assertEqual(w.max_efficacy, 0.0)
        self.assertEqual(w.mean_hotspot_efficacy, 0.0)
        self.assertEqual(w.largest_gap, 501)
        self.assertEqual(w.hotspot_loci, "")

    def test_loci_string_follows_the_classification_tsv_convention(self):
        w = hs.score_window("NM_1", self.region, 1000, 2099, self.spots)
        self.assertEqual(w.hotspot_loci, "1000-1099(10);2000-2099(20)")


class TestRanking(unittest.TestCase):
    def _w(self, start, end, n, gap=0, transcript="NM_1"):
        length = end - start + 1
        return hs.Window(
            transcript, 1, 10000, start, end, length, False, n, n * 10,
            n * 1000 / length, 0.9, 0.9, 0, 0.0, gap, "",
        )

    def test_most_hotspots_wins(self):
        ranked = hs.rank_windows([self._w(1, 2000, 3), self._w(3000, 5000, 7)])
        self.assertEqual(ranked[0].n_hotspots, 7)

    def test_a_tie_goes_to_the_shorter_window(self):
        """Otherwise the length sweep always crowns the longest window."""
        ranked = hs.rank_windows([self._w(1, 2000, 5), self._w(1, 1500, 5)])
        self.assertEqual(ranked[0].length, 1500)

    def test_a_tie_on_length_goes_to_the_more_evenly_covered_window(self):
        ranked = hs.rank_windows(
            [self._w(1, 1500, 5, gap=900), self._w(3000, 4499, 5, gap=200)]
        )
        self.assertEqual(ranked[0].largest_gap, 200)

    def test_ranking_is_deterministic_on_a_full_tie(self):
        a, b = self._w(3000, 4499, 5), self._w(1, 1500, 5)
        self.assertEqual(hs.rank_windows([a, b]), hs.rank_windows([b, a]))

    def test_transcripts_are_grouped(self):
        ranked = hs.rank_windows(
            [self._w(1, 1500, 1, transcript="NM_2"), self._w(1, 1500, 9)]
        )
        self.assertEqual([w.transcript for w in ranked], ["NM_1", "NM_2"])


class TestNonOverlappingSelection(unittest.TestCase):
    def _w(self, start, end, n, transcript="NM_1"):
        length = end - start + 1
        return hs.Window(
            transcript, 1, 10000, start, end, length, False, n, n * 10,
            n * 1000 / length, 0.9, 0.9, 0, 0.0, 0, "",
        )

    def test_shifted_copies_of_one_region_collapse_to_one(self):
        windows = [self._w(s, s + 1499, 5) for s in range(1, 500, 50)]
        self.assertEqual(len(hs.select_non_overlapping(windows)), 1)

    def test_a_disjoint_second_region_survives(self):
        windows = [self._w(1, 1500, 5), self._w(2000, 3499, 4)]
        self.assertEqual(len(hs.select_non_overlapping(windows)), 2)

    def test_abutting_windows_do_not_overlap(self):
        windows = [self._w(1, 1500, 5), self._w(1501, 3000, 4)]
        self.assertEqual(len(hs.select_non_overlapping(windows)), 2)

    def test_one_base_of_overlap_is_overlap(self):
        windows = [self._w(1, 1500, 5), self._w(1500, 2999, 4)]
        self.assertEqual(len(hs.select_non_overlapping(windows)), 1)

    def test_top_is_applied_per_transcript(self):
        windows = [
            self._w(1, 1500, 5), self._w(2000, 3499, 4), self._w(4000, 5499, 3),
            self._w(1, 1500, 5, transcript="NM_2"),
            self._w(2000, 3499, 4, transcript="NM_2"),
        ]
        kept = hs.select_non_overlapping(windows, top=2)
        self.assertEqual(sum(1 for w in kept if w.transcript == "NM_1"), 2)
        self.assertEqual(sum(1 for w in kept if w.transcript == "NM_2"), 2)

    def test_the_best_window_is_always_kept(self):
        windows = [self._w(1, 1500, 2), self._w(100, 1599, 9)]
        self.assertEqual(hs.select_non_overlapping(windows, top=1)[0].n_hotspots, 9)


class TestFilterByCore(unittest.TestCase):
    """Even coverage as a hard requirement, not a tie-break."""

    def _w(self, n, gap, transcript="NM_1"):
        return hs.Window(
            transcript, 1, 10000, 1, 2000, 2000, False, n, n * 10,
            1.0, 0.9, 0.9, 0, 0.0, gap, "",
        )

    def test_none_keeps_everything(self):
        windows = [self._w(5, 50), self._w(5, 5000)]
        self.assertEqual(len(hs.filter_by_core(windows, None)), 2)

    def test_a_window_with_too_long_a_dead_zone_is_dropped(self):
        windows = [self._w(5, 50), self._w(9, 500)]
        kept = hs.filter_by_core(windows, 100)
        self.assertEqual([w.largest_gap for w in kept], [50])

    def test_the_limit_is_inclusive(self):
        self.assertEqual(len(hs.filter_by_core([self._w(5, 100)], 100)), 1)
        self.assertEqual(len(hs.filter_by_core([self._w(5, 101)], 100)), 0)

    def test_it_can_override_the_count_ranking(self):
        """A gappy window winning on count is exactly what this exists to stop."""
        gappy_winner, even_runner_up = self._w(23, 500), self._w(21, 50)
        both = [gappy_winner, even_runner_up]
        self.assertEqual(hs.rank_windows(both)[0].n_hotspots, 23)
        kept = hs.filter_by_core(both, 100)
        self.assertEqual(hs.rank_windows(kept)[0].n_hotspots, 21)


class TestRegions(unittest.TestCase):
    def test_minus_strand_rows_are_ordered(self):
        frame = pd.DataFrame(
            [{"reference": "NM_1", "ref_start": 900, "ref_end": 100}]
        )
        self.assertEqual(hs.regions_from_table(frame), [hs.Region("NM_1", 100, 900)])

    def test_unplaced_rows_are_skipped(self):
        frame = pd.DataFrame(
            [
                {"reference": "NM_1", "ref_start": 100, "ref_end": 900},
                {"reference": "NM_2", "ref_start": None, "ref_end": None},
            ]
        )
        self.assertEqual(hs.regions_from_table(frame), [hs.Region("NM_1", 100, 900)])

    def test_duplicate_rows_collapse(self):
        frame = pd.DataFrame(
            [{"reference": "NM_1", "ref_start": 1, "ref_end": 900}] * 3
        )
        self.assertEqual(len(hs.regions_from_table(frame)), 1)

    def test_whole_transcript_fallback_spans_the_tiling(self):
        df = tiles("NM_1", 500, lambda s: 0.5)
        self.assertEqual(
            hs.whole_transcript_regions(df), [hs.Region("NM_1", 1, 500)]
        )


class TestUnplaceableHotspots(unittest.TestCase):
    def setUp(self):
        self.regions = [hs.Region("NM_1", 1000, 3000)]

    def test_a_hotspot_outside_every_region_is_reported(self):
        spot = hs.Hotspot("NM_1", 5000, 5099, 1, 0.95, 0.95, "x")
        self.assertEqual(hs.unplaceable_hotspots([spot], self.regions), [spot])

    def test_a_hotspot_straddling_a_boundary_is_reported(self):
        spot = hs.Hotspot("NM_1", 2950, 3050, 1, 0.95, 0.95, "x")
        self.assertEqual(hs.unplaceable_hotspots([spot], self.regions), [spot])

    def test_a_contained_hotspot_is_not_reported(self):
        spot = hs.Hotspot("NM_1", 1500, 1599, 1, 0.95, 0.95, "x")
        self.assertEqual(hs.unplaceable_hotspots([spot], self.regions), [])

    def test_a_hotspot_on_another_transcript_is_reported(self):
        spot = hs.Hotspot("NM_2", 1500, 1599, 1, 0.95, 0.95, "x")
        self.assertEqual(hs.unplaceable_hotspots([spot], self.regions), [spot])

    def test_oversized_hotspots_are_reported(self):
        wide = hs.Hotspot("NM_1", 1, 2500, 1, 0.95, 0.95, "x")
        narrow = hs.Hotspot("NM_1", 1, 100, 1, 0.95, 0.95, "y")
        self.assertEqual(hs.oversized_hotspots([wide, narrow], 2000), [wide])


class TestScanRegions(unittest.TestCase):
    """The end-to-end library path, on a transcript with a planted answer."""

    def setUp(self):
        # 6000 nt. Hotspots every 500 nt from 1000 to 2500 (4 of them, dense),
        # then every 1500 nt from 4000 (2 of them, sparse). The densest 1500 bp
        # window must land on the first block.
        dense = [(s, s + 9) for s in (1000, 1500, 2000, 2500)]
        sparse = [(s, s + 9) for s in (4000, 5500)]
        self.df = tiles("NM_1", 6000, plateau(dense + sparse))
        self.hotspots = hs.find_hotspots(self.df, min_efficacy=0.9)

    def test_the_planted_hotspots_are_all_found(self):
        self.assertEqual(len(self.hotspots), 6)

    def test_the_best_window_lands_on_the_dense_block(self):
        regions = [hs.Region("NM_1", 1, 6000)]
        best = hs.select_non_overlapping(
            hs.scan_regions(self.hotspots, regions), top=1
        )[0]
        self.assertEqual(best.n_hotspots, 4)
        self.assertLessEqual(best.start, 1000)
        self.assertGreaterEqual(best.end, 2539)

    def test_a_region_confines_the_answer_even_when_a_better_window_exists(self):
        """The whole point of --regions: conservation beats raw hotspot count."""
        regions = [hs.Region("NM_1", 3500, 6000)]
        best = hs.select_non_overlapping(
            hs.scan_regions(self.hotspots, regions), top=1
        )[0]
        self.assertGreaterEqual(best.start, 3500)
        self.assertLessEqual(best.end, 6000)
        self.assertEqual(best.n_hotspots, 2)

    def test_the_tie_break_prefers_the_shorter_window(self):
        # The dense block spans 1000..2539, so 1750 and 2000 bp windows both
        # reach 4 hotspots; the 1500 bp one cannot.
        regions = [hs.Region("NM_1", 1, 6000)]
        best = hs.select_non_overlapping(
            hs.scan_regions(self.hotspots, regions), top=1
        )[0]
        self.assertEqual(best.length, 1750)


class TestSweep(unittest.TestCase):
    def test_a_threshold_that_changes_the_answer_shows_up_as_two_spans(self):
        # Block A is merely good (0.92), block B is excellent (0.97). At 0.9
        # both qualify and A's three sites win; at 0.95 only B is left.
        def score(start: int) -> float:
            if any(a <= start <= a + 9 for a in (500, 900, 1300)):
                return 0.92
            if any(a <= start <= a + 9 for a in (4000, 4400)):
                return 0.97
            return 0.1

        df = tiles("NM_1", 6000, score)
        sweep = hs.sweep_thresholds(df, [0.9, 0.95])
        self.assertEqual(len(sweep), 2)
        self.assertEqual(sorted(sweep["n_hotspots"]), [2, 3])
        spans = {(r.window_start, r.window_end) for r in sweep.itertuples()}
        self.assertEqual(len(spans), 2)


class TestSweepHonoursTheCoreFilter(unittest.TestCase):
    def test_max_gap_is_applied_inside_the_sweep(self):
        """Otherwise the sweep reports on a ranking the caller never saw.

        Hotspots every 200 bp across 1..1600, so a 1500 bp window is evenly
        covered and its largest gap is ~170 bp: 300 admits it, 100 does not.
        """
        df = tiles("NM_1", 2000, plateau([(s, s + 9) for s in range(1, 1601, 200)]))
        loose = hs.sweep_thresholds(df, [0.9])
        self.assertEqual(len(loose), 1)
        gap = int(loose.iloc[0]["largest_gap"])
        self.assertTrue(100 < gap <= 300, f"fixture gap {gap} does not straddle")

        self.assertEqual(len(hs.sweep_thresholds(df, [0.9], max_gap=300)), 1)
        # Rejected everywhere, so the transcript simply has no row.
        self.assertTrue(hs.sweep_thresholds(df, [0.9], max_gap=100).empty)


class TestFrames(unittest.TestCase):
    def test_rank_restarts_per_transcript(self):
        windows = []
        for transcript in ("NM_1", "NM_2"):
            for i, n in enumerate((5, 3, 1)):
                start = 1 + i * 2000
                windows.append(
                    hs.Window(
                        transcript, 1, 10000, start, start + 1499, 1500, False,
                        n, n * 10, 1.0, 0.9, 0.9, 0, 0.0, 0, "",
                    )
                )
        frame = hs.windows_to_frame(windows)
        self.assertEqual(list(frame[frame.transcript_id == "NM_1"]["rank"]), [1, 2, 3])
        self.assertEqual(list(frame[frame.transcript_id == "NM_2"]["rank"]), [1, 2, 3])
        self.assertEqual(list(frame[frame.transcript_id == "NM_1"]["n_hotspots"]), [5, 3, 1])

    def test_empty_frames_keep_their_columns(self):
        self.assertIn("n_hotspots", hs.windows_to_frame([]).columns)
        self.assertIn("hotspot_start", hs.hotspots_to_frame([]).columns)

    def test_hotspot_frame_carries_length(self):
        frame = hs.hotspots_to_frame(
            [hs.Hotspot("NM_1", 100, 199, 5, 0.95, 0.93, "x")]
        )
        self.assertEqual(frame.loc[0, "length"], 100)


class TestCli(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        dense = [(s, s + 9) for s in (1000, 1500, 2000, 2500)]
        sparse = [(s, s + 9) for s in (4000, 5500)]
        self.predictions = os.path.join(self.tmp, "pred.csv")
        tiles("NM_1", 6000, plateau(dense + sparse)).to_csv(
            self.predictions, index=False
        )
        self.regions = os.path.join(self.tmp, "regions.blastn")
        pd.DataFrame(
            [{"reference": "NM_1", "ref_start": 1, "ref_end": 6000}]
        ).to_csv(self.regions, sep="\t", index=False)
        self.output = os.path.join(self.tmp, "out.tsv")

    def run_cli(self, *extra):
        argv = [
            "--predictions", self.predictions,
            "--regions", self.regions,
            "--output", self.output,
            *extra,
        ]
        err = io.StringIO()
        with redirect_stderr(err):
            code = find_hotspot_regions.main(argv)
        self.assertEqual(code, 0)
        return pd.read_csv(self.output, sep="\t"), err.getvalue()

    def test_it_writes_the_expected_winner(self):
        frame, _ = self.run_cli("--top", "1")
        self.assertEqual(len(frame), 1)
        row = frame.iloc[0]
        self.assertEqual(row["n_hotspots"], 4)
        self.assertEqual(row["rank"], 1)
        self.assertEqual(row["window_length"], 1750)

    def test_loci_round_trip_to_the_hotspot_dump(self):
        dump = os.path.join(self.tmp, "hotspots.tsv")
        frame, _ = self.run_cli("--top", "1", "--hotspots-out", dump)
        spots = pd.read_csv(dump, sep="\t")
        loci = frame.iloc[0]["hotspot_loci"].split(";")
        self.assertEqual(len(loci), 4)
        for locus in loci:
            span = locus.split("(")[0]
            start, end = (int(x) for x in span.split("-"))
            self.assertTrue(
                ((spots.hotspot_start == start) & (spots.hotspot_end == end)).any(),
                f"{locus} is absent from {dump}",
            )

    def test_all_windows_holds_more_rows_than_the_suppressed_output(self):
        every = os.path.join(self.tmp, "all.tsv")
        frame, _ = self.run_cli("--top", "1", "--all-windows", every)
        self.assertGreater(len(pd.read_csv(every, sep="\t")), len(frame))

    def test_running_without_regions_warns_that_conservation_is_not_guaranteed(self):
        err = io.StringIO()
        with redirect_stderr(err):
            find_hotspot_regions.main(
                ["--predictions", self.predictions, "--output", self.output]
            )
        self.assertIn("NOT guaranteed to be conserved", err.getvalue())

    def test_a_region_shorter_than_min_window_is_reported_not_skipped(self):
        short = os.path.join(self.tmp, "short.blastn")
        pd.DataFrame(
            [{"reference": "NM_1", "ref_start": 1000, "ref_end": 1800}]
        ).to_csv(short, sep="\t", index=False)
        err = io.StringIO()
        with redirect_stderr(err):
            find_hotspot_regions.main(
                ["--predictions", self.predictions, "--regions", short,
                 "--output", self.output]
            )
        frame = pd.read_csv(self.output, sep="\t")
        self.assertEqual(len(frame), 1)
        self.assertEqual(frame.iloc[0]["window_length"], 801)
        self.assertTrue(bool(frame.iloc[0]["undersized"]))
        self.assertIn("shorter than --min-window", err.getvalue())

    def test_unplaceable_hotspots_are_counted_in_the_summary(self):
        partial = os.path.join(self.tmp, "partial.blastn")
        pd.DataFrame(
            [{"reference": "NM_1", "ref_start": 1, "ref_end": 3000}]
        ).to_csv(partial, sep="\t", index=False)
        err = io.StringIO()
        with redirect_stderr(err):
            find_hotspot_regions.main(
                ["--predictions", self.predictions, "--regions", partial,
                 "--output", self.output]
            )
        self.assertIn("lie outside every region", err.getvalue())

    def test_a_threshold_above_every_score_fails_loudly(self):
        with self.assertRaises(SystemExit) as caught:
            find_hotspot_regions.main(
                ["--predictions", self.predictions, "--regions", self.regions,
                 "--output", self.output, "--min-efficacy", "0.999"]
            )
        self.assertIn("no guide", str(caught.exception))

    def test_a_missing_column_names_the_file(self):
        bad = os.path.join(self.tmp, "bad.csv")
        pd.DataFrame([{"transcript_id": "NM_1", "start": 1}]).to_csv(bad, index=False)
        with self.assertRaises(SystemExit) as caught:
            find_hotspot_regions.main(
                ["--predictions", bad, "--output", self.output]
            )
        self.assertIn("predicted_efficacy", str(caught.exception))

    def test_a_missing_file_is_reported_before_pandas_sees_it(self):
        with self.assertRaises(SystemExit) as caught:
            find_hotspot_regions.main(
                ["--predictions", os.path.join(self.tmp, "nope.csv"),
                 "--output", self.output]
            )
        self.assertIn("no such file", str(caught.exception))

    def test_min_window_above_max_window_is_rejected(self):
        with self.assertRaises(SystemExit) as caught:
            find_hotspot_regions.main(
                ["--predictions", self.predictions, "--output", self.output,
                 "--min-window", "3000", "--max-window", "2000"]
            )
        self.assertIn("exceeds --max-window", str(caught.exception))

    def test_the_sweep_table_is_written(self):
        sweep = os.path.join(self.tmp, "sweep.tsv")
        self.run_cli("--sweep-thresholds", "0.5,0.9", "--sweep-out", sweep)
        frame = pd.read_csv(sweep, sep="\t")
        self.assertEqual(sorted(frame["min_efficacy"].unique()), [0.5, 0.9])

    def test_sweep_out_defaults_beside_the_output(self):
        self.run_cli("--sweep-thresholds", "0.9")
        self.assertTrue(os.path.isfile(os.path.join(self.tmp, "out.sweep.tsv")))

    def test_an_unchanged_winner_is_called_stable(self):
        _, err = self.run_cli("--sweep-thresholds", "0.5,0.9")
        self.assertIn("a stable choice", err)
        self.assertNotIn("JUMPS", err)

    def test_a_winner_that_moves_to_a_disjoint_span_is_called_out(self):
        # Three merely-good sites early, two excellent ones late: raising the
        # threshold moves the winner to the other end of the transcript.
        def score(start: int) -> float:
            if any(a <= start <= a + 9 for a in (500, 900, 1300)):
                return 0.92
            if any(a <= start <= a + 9 for a in (4000, 4400)):
                return 0.97
            return 0.1

        predictions = os.path.join(self.tmp, "jumpy.csv")
        tiles("NM_1", 6000, score).to_csv(predictions, index=False)
        err = io.StringIO()
        with redirect_stderr(err):
            find_hotspot_regions.main(
                ["--predictions", predictions, "--regions", self.regions,
                 "--output", self.output, "--sweep-thresholds", "0.9,0.95"]
            )
        self.assertIn("JUMPS", err.getvalue())

    def test_a_threshold_with_no_surviving_window_is_not_called_stable(self):
        """One row agrees with itself; that is missing evidence, not stability."""
        _, err = self.run_cli(
            "--max-gap", "500", "--sweep-thresholds", "0.9,0.99"
        )
        self.assertIn("Stability is UNTESTED", err)
        self.assertIn("nothing at 0.99", err)
        self.assertNotIn("a stable choice", err)

    def test_a_malformed_sweep_threshold_list_is_rejected(self):
        with self.assertRaises(SystemExit) as caught:
            find_hotspot_regions.main(
                ["--predictions", self.predictions, "--output", self.output,
                 "--sweep-thresholds", "0.9,high"]
            )
        self.assertIn("comma-separated list", str(caught.exception))

    def test_max_gap_reports_what_it_dropped(self):
        # The planted hotspots sit 461 bp apart, so 500 keeps the dense block
        # and rejects the windows reaching into the sparse tail.
        frame, err = self.run_cli("--max-gap", "500")
        self.assertIn("--max-gap 500:", err)
        self.assertTrue((frame["largest_gap"] <= 500).all())

    def test_an_impossible_max_gap_fails_loudly_rather_than_writing_nothing(self):
        err = io.StringIO()
        with redirect_stderr(err), self.assertRaises(SystemExit) as caught:
            find_hotspot_regions.main(
                ["--predictions", self.predictions, "--regions", self.regions,
                 "--output", self.output, "--max-gap", "100"]
            )
        # The stop happens before anything is written, so the message itself has
        # to say how tight the limit could have been.
        self.assertIn("rejected every one", str(caught.exception))
        self.assertIn("461 bp dead zone", str(caught.exception))
        self.assertEqual(err.getvalue(), "")

    def test_output_is_deterministic(self):
        first, _ = self.run_cli()
        second, _ = self.run_cli()
        pd.testing.assert_frame_equal(first, second)


if __name__ == "__main__":
    unittest.main()
