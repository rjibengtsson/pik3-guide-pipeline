"""Tests for src/visualize.py and scripts/plot_guides.py.

    python -m unittest discover -s test -v

No test here asserts anything about how the plot *looks* -- colours, bar widths
and layout are presentation and will change. What is pinned is the part that can
be silently wrong and would corrupt the reading of the plot: the lane packing,
the 1-based inclusive coordinate convention shared with ``tile_id``, the bounds
check, and the CLI's accession matching and filtering.

Plotly is a hard dependency of the module, so these do not skip; the CLI tests
also need pandas. Both are in ``environment.yml``.
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

import plot_guides  # noqa: E402
import visualize as visualize_mod  # noqa: E402
from visualize import (  # noqa: E402
    Guide,
    check_within_transcript,
    guide_track,
    guides_from_rows,
    pack_lanes,
    write_html,
)

TRANSCRIPT = "NM_000942.5"
LENGTH = 120

#: A spacer and the target it pairs with, kept distinct so a test can tell which
#: of the two reached the hover box.
SPACER = "GCAGCAGCAGGAAGAAGACGGACCCCGCGA"
TARGET = "UCGCGGGGUCCGUCUUCUUCCUGCUGCUGC"

#: Rows shaped like a predictions CSV, with a known right answer: three guides at
#: 30 nt, two of which overlap (1-30 and 15-44) and one that does not (61-90).
ROWS = [
    {
        "tile_id": f"{TRANSCRIPT}:1-30",
        "transcript_id": TRANSCRIPT,
        "start": "1",
        "end": "30",
        "predicted_efficacy": "0.97",
        "class_exact": "specific",
        "class_nm": "specific",
        "spacer_seq": SPACER,
        "target_seq": TARGET,
    },
    {
        "tile_id": f"{TRANSCRIPT}:15-44",
        "transcript_id": TRANSCRIPT,
        "start": "15",
        "end": "44",
        "predicted_efficacy": "0.96",
        "class_exact": "specific",
        "class_nm": "specific",
        "spacer_seq": SPACER,
        "target_seq": TARGET,
    },
    {
        "tile_id": f"{TRANSCRIPT}:61-90",
        "transcript_id": TRANSCRIPT,
        "start": "61",
        "end": "90",
        "predicted_efficacy": "0.41",
        "class_exact": "specific",
        "class_nm": "specific",
        "spacer_seq": SPACER,
        "target_seq": TARGET,
    },
]


def guides(*spans: tuple[int, int]) -> list[Guide]:
    """Guides at the given ``(start, end)`` spans, ids derived from them."""
    return [Guide(f"{TRANSCRIPT}:{start}-{end}", start, end) for start, end in spans]


def lanes_are_non_overlapping(assigned: list[Guide], lanes: list[int], gap: int) -> bool:
    """Whether no two guides sharing a lane are closer than ``gap`` bases.

    Written as a direct pairwise check rather than in terms of :func:`pack_lanes`'
    own logic, for the same reason ``pairs_antiparallel`` avoids
    ``reverse_complement``: it must fail, not agree, when the packer is wrong.
    """
    by_lane: dict[int, list[Guide]] = {}
    for guide, lane in zip(assigned, lanes):
        by_lane.setdefault(lane, []).append(guide)
    for members in by_lane.values():
        members.sort(key=lambda guide: guide.start)
        for earlier, later in zip(members, members[1:]):
            if earlier.end + gap >= later.start:
                return False
    return True


class TestGuide(unittest.TestCase):
    def test_length_is_inclusive_of_both_endpoints(self):
        """The coordinate convention of src/tile.py: sequence[start - 1:end]."""
        self.assertEqual(Guide("t:1-30", 1, 30).length, 30)
        self.assertEqual(Guide("t:5-5", 5, 5).length, 1)


class TestGuidesFromRows(unittest.TestCase):
    def test_coordinates_and_scores_are_coerced_from_strings(self):
        """csv.DictReader yields strings; a DataFrame yields numbers. Both work."""
        built = guides_from_rows(ROWS)
        self.assertEqual([guide.start for guide in built], [1, 15, 61])
        self.assertEqual([guide.end for guide in built], [30, 44, 90])
        self.assertEqual(built[0].score, 0.97)
        self.assertIsInstance(built[0].start, int)
        self.assertIsInstance(built[0].score, float)

    def test_rows_are_sorted_by_start(self):
        built = guides_from_rows(list(reversed(ROWS)))
        self.assertEqual([guide.start for guide in built], [1, 15, 61])

    def test_transcript_id_filter_keeps_only_the_named_transcript(self):
        rows = ROWS + [dict(ROWS[0], transcript_id="NM_999999.9", tile_id="other:1-30")]
        built = guides_from_rows(rows, transcript_id=TRANSCRIPT)
        self.assertEqual(len(built), 3)
        self.assertTrue(all(guide.tile_id.startswith(TRANSCRIPT) for guide in built))

    def test_filtering_without_the_column_is_an_error_not_a_silent_pass(self):
        """Plotting every transcript on one axis would look plausible and be wrong."""
        rows = [{"tile_id": "t:1-30", "start": 1, "end": 30}]
        with self.assertRaises(KeyError):
            guides_from_rows(rows, transcript_id=TRANSCRIPT)

    def test_absent_score_column_is_tolerated_but_empty_value_becomes_none(self):
        """The same call has to work on a predictions file and a plain tiles file."""
        no_column = guides_from_rows([{"tile_id": "t:1-30", "start": 1, "end": 30}])
        self.assertIsNone(no_column[0].score)
        empty = guides_from_rows(
            [{"tile_id": "t:1-30", "start": 1, "end": 30, "predicted_efficacy": ""}]
        )
        self.assertIsNone(empty[0].score)

    def test_label_column_is_read_only_when_asked_for(self):
        self.assertIsNone(guides_from_rows(ROWS)[0].label)
        self.assertEqual(guides_from_rows(ROWS, label_column="class_exact")[0].label, "specific")

    def test_sequence_column_is_read_only_when_asked_for(self):
        """Not defaulted: guessing wrong would put the target in a 'guide' box."""
        self.assertIsNone(guides_from_rows(ROWS)[0].sequence)
        built = guides_from_rows(ROWS, sequence_column="spacer_seq")
        self.assertEqual(built[0].sequence, SPACER)


class TestAnnotations(unittest.TestCase):
    """Extra hover lines from a later stage, e.g. stage 2's classification."""

    def test_named_columns_are_carried_in_order(self):
        built = guides_from_rows(ROWS, annotation_columns=("class_exact", "predicted_efficacy"))
        self.assertEqual(
            built[0].annotations, (("class_exact", "specific"), ("predicted_efficacy", "0.97"))
        )

    def test_no_annotation_columns_means_no_annotations(self):
        self.assertEqual(guides_from_rows(ROWS)[0].annotations, ())

    def test_an_absent_column_is_skipped_not_shown_empty(self):
        """What makes a left join usable: unmatched guides carry fewer lines."""
        built = guides_from_rows(ROWS, annotation_columns=("class_exact", "not_a_column"))
        self.assertEqual(built[0].annotations, (("class_exact", "specific"),))

    def test_an_empty_value_is_skipped(self):
        """A left join fills unmatched rows with blanks; they must not show."""
        rows = [dict(ROWS[0], class_exact="")]
        self.assertEqual(guides_from_rows(rows, annotation_columns=("class_exact",))[0].annotations, ())

    def test_annotations_reach_the_hover_box(self):
        figure = guide_track(
            LENGTH,
            [Guide("t:1-30", 1, 30, 0.9, None, None, (("class_exact", "specific"),))],
            transcript_id=TRANSCRIPT,
        )
        self.assertIn("class_exact", figure.data[0].customdata[0][7])
        self.assertIn("specific", figure.data[0].customdata[0][7])

    def test_a_guide_without_annotations_adds_no_hover_line(self):
        figure = guide_track(LENGTH, guides((1, 30)), transcript_id=TRANSCRIPT)
        self.assertEqual(figure.data[0].customdata[0][7], "")

    def test_a_long_multi_reference_value_is_broken_at_semicolons(self):
        """';' divides references, so one reference per line reads correctly."""
        locus = "NM_005026.5:2754-2783(0);NM_006219.3:2978-3007(3)"
        rendered = visualize_mod._annotation_lines((("hit_locs_nm", locus),))
        self.assertIn(";<br>", rendered)
        self.assertEqual(rendered.count("<br>"), 2)  # the leading one, plus the split

    def test_several_loci_on_one_reference_are_never_broken(self):
        """',' means a self-repeat; splitting it would scatter the distinction.

        This pair -- ';' broken, ',' not -- is the whole point: it is what lets a
        reader tell cross-family reactivity from a tile matching its own
        transcript twice (see the stage 2 notes in CLAUDE.md).
        """
        locus = "NM_002649.3:5952-5981(0),NM_002649.3:6004-6033(3)"
        rendered = visualize_mod._annotation_lines((("hit_locs_nm", locus),))
        self.assertNotIn(",<br>", rendered)
        self.assertIn(locus, rendered)

    def test_a_short_value_is_left_on_one_line(self):
        """A needless break reads as two findings where there is one."""
        rendered = visualize_mod._annotation_lines((("hit_genes_nm", "PIK3CB;PIK3CD"),))
        self.assertNotIn(";<br>", rendered)
        self.assertIn("PIK3CB;PIK3CD", rendered)


class TestPackLanes(unittest.TestCase):
    def test_non_overlapping_guides_share_one_lane(self):
        self.assertEqual(pack_lanes(guides((1, 30), (61, 90))), [0, 0])

    def test_overlapping_guides_are_separated(self):
        self.assertEqual(pack_lanes(guides((1, 30), (15, 44))), [0, 1])

    def test_abutting_guides_are_separated_by_the_default_gap(self):
        """end=30 then start=31 would otherwise draw as one unbroken bar."""
        self.assertEqual(pack_lanes(guides((1, 30), (31, 60)))[1], 1)
        self.assertEqual(pack_lanes(guides((1, 30), (31, 60)), gap=0), [0, 0])

    def test_lane_count_equals_the_depth_of_overlap(self):
        """The claim the plot's y axis makes: lanes at a position == guides there."""
        stepped = guides(*[(start, start + 29) for start in range(1, 31)])
        lanes = pack_lanes(stepped)
        self.assertEqual(max(lanes) + 1, 30)

    def test_no_two_guides_in_a_lane_overlap_on_a_dense_tiling(self):
        """Checked independently of how pack_lanes decides (see module docstring)."""
        dense = guides(*[(start, start + 29) for start in range(1, 200, 3)])
        for gap in (0, 1, 5):
            lanes = pack_lanes(dense, gap=gap)
            self.assertTrue(
                lanes_are_non_overlapping(dense, lanes, gap),
                f"overlap within a lane at gap={gap}",
            )

    def test_lanes_are_returned_in_input_order(self):
        """Lanes are positionally matched to guides, which may be unsorted."""
        unsorted = guides((61, 90), (1, 30), (15, 44))
        self.assertEqual(len(pack_lanes(unsorted)), 3)
        self.assertEqual(pack_lanes(unsorted)[0], 0)

    def test_empty_input_uses_no_lanes(self):
        self.assertEqual(pack_lanes([]), [])

    def test_negative_gap_is_rejected(self):
        with self.assertRaises(ValueError):
            pack_lanes(guides((1, 30)), gap=-1)


class TestCheckWithinTranscript(unittest.TestCase):
    """A guide off the end means the guides and the length are different sequences."""

    def test_guides_inside_the_transcript_pass(self):
        check_within_transcript(guides((1, 30), (91, 120)), LENGTH)

    def test_a_guide_running_past_the_end_is_rejected(self):
        with self.assertRaises(ValueError) as caught:
            check_within_transcript(guides((100, 130)), LENGTH)
        self.assertIn(f"{TRANSCRIPT}:100-130", str(caught.exception))

    def test_a_guide_before_the_first_base_is_rejected(self):
        """Coordinates are 1-based, so 0 is not a position."""
        with self.assertRaises(ValueError):
            check_within_transcript(guides((0, 29)), LENGTH)

    def test_an_inverted_span_is_rejected(self):
        with self.assertRaises(ValueError):
            check_within_transcript(guides((50, 20)), LENGTH)

    def test_a_zero_length_transcript_is_rejected(self):
        with self.assertRaises(ValueError):
            check_within_transcript([], 0)


class TestGuideTrack(unittest.TestCase):
    """The figure. Only the data behind it is asserted, never its appearance."""

    def test_bars_span_each_guide_inclusive_of_both_end_bases(self):
        """base..base+x must be start-0.5..end+0.5, or the plot misreports position."""
        figure = guide_track(LENGTH, guides((1, 30), (61, 90)), transcript_id=TRANSCRIPT)
        bars = figure.data[0]
        self.assertEqual(list(bars.base), [0.5, 60.5])
        self.assertEqual(list(bars.x), [30, 30])
        self.assertEqual([base + width for base, width in zip(bars.base, bars.x)], [30.5, 90.5])

    def test_x_axis_covers_the_whole_transcript_not_just_the_guides(self):
        """Uncovered regions have to read as uncovered -- the point of the plot."""
        figure = guide_track(LENGTH, guides((1, 30)), transcript_id=TRANSCRIPT)
        self.assertEqual(tuple(figure.layout.xaxis.range), (0.5, LENGTH + 0.5))

    def test_hover_carries_the_tile_id_and_span_of_every_guide(self):
        figure = guide_track(LENGTH, guides((1, 30), (61, 90)), transcript_id=TRANSCRIPT)
        ids = [row[0] for row in figure.data[0].customdata]
        self.assertEqual(ids, [f"{TRANSCRIPT}:1-30", f"{TRANSCRIPT}:61-90"])
        # Plotly normalises the tuples it is given to lists.
        self.assertEqual(
            [list(row[1:4]) for row in figure.data[0].customdata], [[1, 30, 30], [61, 90, 30]]
        )

    def test_a_label_is_shown_only_when_present(self):
        """An absent label must add nothing to the hover box, not a blank line."""
        figure = guide_track(
            LENGTH,
            [Guide("t:1-30", 1, 30, 0.9, "specific"), Guide("t:61-90", 61, 90, 0.9, None)],
            transcript_id=TRANSCRIPT,
        )
        self.assertEqual(figure.data[0].customdata[0][5], "<br>specific")
        self.assertEqual(figure.data[0].customdata[1][5], "")

    def test_the_guide_sequence_is_shown_in_the_hover_box(self):
        figure = guide_track(
            LENGTH,
            [Guide("t:1-30", 1, 30, 0.9, None, SPACER)],
            transcript_id=TRANSCRIPT,
        )
        cell = figure.data[0].customdata[0][6]
        self.assertIn(SPACER, cell)
        self.assertIn("5'-", cell)
        self.assertIn("-3'", cell)

    def test_a_guide_without_a_sequence_adds_no_hover_line(self):
        """An absent sequence must add nothing, not a blank line or a stray 5'-."""
        figure = guide_track(LENGTH, guides((1, 30)), transcript_id=TRANSCRIPT)
        self.assertEqual(figure.data[0].customdata[0][6], "")

    def test_the_hovertemplate_references_every_customdata_field(self):
        """A field added to customdata but not the template is invisible."""
        figure = guide_track(
            LENGTH,
            [Guide("t:1-30", 1, 30, 0.9, "specific", SPACER)],
            transcript_id=TRANSCRIPT,
        )
        template = figure.data[0].hovertemplate
        for index in range(len(figure.data[0].customdata[0])):
            self.assertIn(f"customdata[{index}]", template)

    def test_lane_mode_is_the_default_and_puts_overlapping_guides_on_separate_rows(self):
        figure = guide_track(LENGTH, guides((1, 30), (15, 44)), transcript_id=TRANSCRIPT)
        self.assertEqual(list(figure.data[0].y), [0, 1])

    def test_score_mode_puts_each_guide_at_the_height_of_its_score(self):
        scored = [Guide("t:1-30", 1, 30, 0.97), Guide("t:61-90", 61, 90, 0.41)]
        figure = guide_track(LENGTH, scored, transcript_id=TRANSCRIPT, y_mode="score")
        self.assertEqual(list(figure.data[0].y), [0.97, 0.41])

    def test_score_mode_without_scores_is_rejected(self):
        with self.assertRaises(ValueError) as caught:
            guide_track(LENGTH, guides((1, 30)), y_mode="score")
        self.assertIn("score", str(caught.exception))

    def test_an_unknown_y_mode_is_rejected(self):
        with self.assertRaises(ValueError):
            guide_track(LENGTH, guides((1, 30)), y_mode="sideways")

    def test_no_guides_still_plots_the_transcript(self):
        """A threshold nothing passes is a real answer, not an error."""
        figure = guide_track(LENGTH, [], transcript_id=TRANSCRIPT)
        self.assertEqual(len(figure.data), 0)
        self.assertEqual(tuple(figure.layout.xaxis.range), (0.5, LENGTH + 0.5))
        self.assertTrue(any("no guides" in note.text for note in figure.layout.annotations))

    def test_out_of_bounds_guides_are_rejected_before_anything_is_drawn(self):
        with self.assertRaises(ValueError):
            guide_track(LENGTH, guides((100, 130)), transcript_id=TRANSCRIPT)


class TestWriteHtml(unittest.TestCase):
    def test_a_standalone_file_is_written_and_its_path_returned(self):
        figure = guide_track(LENGTH, guides((1, 30)), transcript_id=TRANSCRIPT)
        with tempfile.TemporaryDirectory() as tmp:
            path = write_html(figure, os.path.join(tmp, "guides.html"))
            self.assertTrue(os.path.exists(path))
            with open(path, encoding="utf-8") as handle:
                html = handle.read()
            self.assertIn("<div", html)
            self.assertIn(f"{TRANSCRIPT}:1-30", html)


class TestCli(unittest.TestCase):
    """End to end over a written-out FASTA and predictions CSV."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.fasta = os.path.join(self.tmp.name, f"{TRANSCRIPT}.fasta")
        self.write_fasta({TRANSCRIPT: "ACGT" * (LENGTH // 4)})
        self.predictions = os.path.join(self.tmp.name, "predictions.csv")
        self.write_predictions(ROWS)
        self.output = os.path.join(self.tmp.name, "out", "guides.html")

    def write_fasta(self, sequences: dict[str, str]) -> None:
        with open(self.fasta, "w", encoding="utf-8") as handle:
            for accession, sequence in sequences.items():
                handle.write(f">{accession} Homo sapiens synthetic, mRNA\n{sequence}\n")

    def write_predictions(self, rows: list[dict[str, str]]) -> None:
        with open(self.predictions, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    def run_cli(self, *extra: str) -> str:
        """Call main(), capturing the stderr summary."""
        stderr = io.StringIO()
        old = sys.stderr
        sys.stderr = stderr
        try:
            self.assertEqual(plot_guides.main(list(extra)), 0)
        finally:
            sys.stderr = old
        return stderr.getvalue()

    def read_output(self) -> str:
        with open(self.output, encoding="utf-8") as handle:
            return handle.read()

    def plot(self, *extra: str, efficacy: str = "0.95") -> str:
        return self.run_cli(
            "--fasta", self.fasta,
            "--predictions", self.predictions,
            "--min-efficacy", efficacy,
            "--output", self.output,
            *extra,
        )

    def test_length_comes_from_the_fasta_and_sets_the_axis(self):
        summary = self.plot()
        self.assertIn(f"length: {LENGTH:,} nt", summary)

    def test_the_threshold_filters_the_predictions(self):
        """0.97 and 0.96 pass at 0.95; 0.41 does not."""
        summary = self.plot()
        self.assertIn("kept 2", summary)
        with open(self.output, encoding="utf-8") as handle:
            html = handle.read()
        self.assertIn(f"{TRANSCRIPT}:1-30", html)
        self.assertNotIn(f"{TRANSCRIPT}:61-90", html)

    def test_the_threshold_boundary_is_inclusive(self):
        self.assertIn("kept 1", self.plot(efficacy="0.97"))
        self.assertIn("kept 3", self.plot(efficacy="0.41"))

    def test_only_rows_matching_the_fasta_accession_are_plotted(self):
        other = dict(ROWS[0], transcript_id="NM_999999.9", tile_id="NM_999999.9:1-30")
        self.write_predictions(ROWS + [other])
        summary = self.plot()
        self.assertIn("3 row(s) for this transcript", summary)
        with open(self.output, encoding="utf-8") as handle:
            self.assertNotIn("NM_999999.9", handle.read())

    def test_a_predictions_file_for_another_dataset_is_a_hard_error(self):
        """Nearly always a wrong-file mistake; an empty plot would look plausible."""
        self.write_predictions(
            [dict(row, transcript_id="NM_999999.9") for row in ROWS]
        )
        with self.assertRaises(SystemExit) as caught:
            self.plot()
        self.assertIn("NM_999999.9", str(caught.exception))

    def test_a_threshold_nothing_passes_is_not_an_error(self):
        summary = self.plot(efficacy="0.999")
        self.assertIn("kept 0", summary)
        self.assertIn("0.97", summary)  # the highest available, so the user can retry
        self.assertTrue(os.path.exists(self.output))

    def test_a_multi_record_fasta_is_rejected(self):
        self.write_fasta({TRANSCRIPT: "ACGT" * 30, "NM_999999.9": "TGCA" * 30})
        with self.assertRaises(SystemExit) as caught:
            self.plot()
        self.assertIn("2 records", str(caught.exception))

    def test_a_missing_required_column_is_named(self):
        self.write_predictions([{key: value for key, value in row.items() if key != "end"} for row in ROWS])
        with self.assertRaises(SystemExit) as caught:
            self.plot()
        self.assertIn("end", str(caught.exception))

    def test_a_missing_input_file_is_reported_not_traced(self):
        with self.assertRaises(SystemExit) as caught:
            self.run_cli(
                "--fasta", os.path.join(self.tmp.name, "absent.fasta"),
                "--predictions", self.predictions,
                "--min-efficacy", "0.95",
            )
        self.assertIn("no such FASTA", str(caught.exception))

    def test_an_unknown_label_column_is_rejected(self):
        with self.assertRaises(SystemExit) as caught:
            self.plot("--label-column", "nope")
        self.assertIn("nope", str(caught.exception))

    def test_the_default_output_name_carries_the_accession_and_threshold(self):
        """Two thresholds must not overwrite one another."""
        summary = self.run_cli(
            "--fasta", self.fasta,
            "--predictions", self.predictions,
            "--min-efficacy", "0.96",
        )
        expected = os.path.join(self.tmp.name, f"guides-{TRANSCRIPT}-eff0.96.html")
        self.assertIn(expected, summary)
        self.assertTrue(os.path.exists(expected))

    def test_a_non_numeric_efficacy_is_warned_about_and_excluded(self):
        self.write_predictions([dict(ROWS[0], predicted_efficacy="n/a")] + ROWS[1:])
        summary = self.plot()
        self.assertIn("warning:", summary)
        self.assertIn("kept 1", summary)

    def test_score_mode_is_selectable(self):
        self.assertIn("y axis: score", self.plot("--y-mode", "score"))

    def test_the_spacer_is_auto_detected_and_reaches_the_html(self):
        summary = self.plot()
        self.assertIn("hover sequence: spacer_seq", summary)
        with open(self.output, encoding="utf-8") as handle:
            html = handle.read()
        self.assertIn(SPACER, html)

    def test_the_guide_is_shown_not_the_target(self):
        """target_seq is never auto-detected: it is the target, not the guide.

        The fixture carries both columns, so this fails if the resolver ever
        reaches for target_seq.
        """
        self.plot()
        with open(self.output, encoding="utf-8") as handle:
            html = handle.read()
        self.assertIn(SPACER, html)
        self.assertNotIn(TARGET, html)

    def test_guide_rna_is_used_when_there_is_no_spacer_seq(self):
        """A stage-1 tiles file names the column differently."""
        self.write_predictions(
            [
                {key: value for key, value in dict(row, guide_rna=row["spacer_seq"]).items()
                 if key != "spacer_seq"}
                for row in ROWS
            ]
        )
        self.assertIn("hover sequence: guide_rna", self.plot())

    def test_no_sequence_column_is_noted_rather_than_silently_dropped(self):
        self.write_predictions(
            [
                {key: value for key, value in row.items() if key not in ("spacer_seq",)}
                for row in ROWS
            ]
        )
        summary = self.plot()
        self.assertIn("hover sequence: none", summary)
        self.assertIn("note:", summary)

    def test_the_sequence_can_be_turned_off(self):
        summary = self.plot("--sequence-column", "none")
        self.assertIn("hover sequence: none", summary)
        with open(self.output, encoding="utf-8") as handle:
            self.assertNotIn(SPACER, handle.read())

    def write_classification(self, rows: list[dict[str, str]]) -> str:
        path = os.path.join(self.tmp.name, "classification.tsv")
        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
        return path

    def classification_rows(self, **overrides: str) -> list[dict[str, str]]:
        return [
            dict(
                {
                    "tile_id": row["tile_id"],
                    "class_nm": "specific",
                    "hit_genes_nm": "PIK3CA",
                    "hit_locs_nm": f"{row['tile_id']}(0)",
                },
                **overrides,
            )
            for row in ROWS
        ]

    def test_without_the_flag_nothing_changes(self):
        """The file is optional: the plot must be identical without it."""
        summary = self.plot()
        self.assertNotIn("classification", summary)
        with open(self.output, encoding="utf-8") as handle:
            self.assertNotIn("class_nm", handle.read())

    def test_the_classification_columns_reach_the_hover_box(self):
        path = self.write_classification(self.classification_rows())
        summary = self.plot("--classification", path)
        self.assertIn("2 of 2 plotted guide(s) annotated", summary)
        with open(self.output, encoding="utf-8") as handle:
            html = handle.read()
        for column in plot_guides.CLASSIFICATION_COLUMNS:
            self.assertIn(column, html)

    def test_the_advisory_nature_of_the_nm_columns_is_stated(self):
        """The NM<=3 columns are unvalidated; the summary has to say so."""
        path = self.write_classification(self.classification_rows())
        summary = self.plot("--classification", path)
        self.assertIn("advisory", summary)

    def test_guides_missing_from_the_classification_are_kept_and_warned_about(self):
        """A left join: the plot shows the filtered guides regardless."""
        rows = self.classification_rows()
        path = self.write_classification(rows[:1])
        summary = self.plot("--classification", path)
        self.assertIn("1 of 2 plotted guide(s) annotated", summary)
        self.assertIn("warning:", summary)
        self.assertIn(f"{TRANSCRIPT}:15-44", self.read_output())

    def test_an_unannotated_guide_simply_carries_fewer_lines(self):
        path = self.write_classification(self.classification_rows()[:1])
        self.plot("--classification", path)
        self.assertIn(f"{TRANSCRIPT}:1-30", self.read_output())

    def test_a_classification_from_another_tiling_is_a_hard_error(self):
        """No tile_id overlap means different files; every box would be blank."""
        rows = [dict(row, tile_id="NM_999999.9:1-30") for row in self.classification_rows()]
        path = self.write_classification(rows)
        with self.assertRaises(SystemExit) as caught:
            self.plot("--classification", path)
        self.assertIn("tile_id", str(caught.exception))

    def test_a_missing_classification_column_is_named(self):
        rows = [
            {key: value for key, value in row.items() if key != "hit_locs_nm"}
            for row in self.classification_rows()
        ]
        path = self.write_classification(rows)
        with self.assertRaises(SystemExit) as caught:
            self.plot("--classification", path)
        self.assertIn("hit_locs_nm", str(caught.exception))

    def test_a_missing_classification_file_is_reported_not_traced(self):
        with self.assertRaises(SystemExit) as caught:
            self.plot("--classification", os.path.join(self.tmp.name, "absent.tsv"))
        self.assertIn("no such classification file", str(caught.exception))

    def test_an_empty_guide_set_with_a_classification_is_not_a_mismatch(self):
        """Nothing passing the threshold must not look like a wrong-file error."""
        path = self.write_classification(self.classification_rows())
        summary = self.plot("--classification", path, efficacy="0.999")
        self.assertIn("kept 0", summary)
        self.assertTrue(os.path.exists(self.output))

    def test_the_classification_file_wins_a_column_name_collision(self):
        """ROWS already carries class_nm, as a predictions file may.

        Left to pandas, the merge suffixes the *incoming* column and the stale
        predictions value keeps the name everything downstream reads -- so the
        hover box and the match count would both report the wrong file.
        """
        path = self.write_classification(self.classification_rows(class_nm="dual"))
        summary = self.plot("--classification", path)
        self.assertIn("2 of 2 plotted guide(s) annotated", summary)
        html = self.read_output()
        self.assertIn("dual", html)
        self.assertNotIn("class_nm_stage2", html)

    def test_a_label_can_name_a_classification_column(self):
        """class_exact is not shown by default; --label-column has to reach it.

        The check runs after the join for this reason, and the join carries the
        named column across even though it is not one of the three annotations.
        """
        # A sentinel, not "pan": short words turn up inside Plotly's own HTML
        # ("span" contains "pan") and the assertion would pass vacuously.
        sentinel = "CLASSIFICATION_FILE_LABEL"
        rows = [dict(row, class_exact=sentinel) for row in self.classification_rows()]
        path = self.write_classification(rows)
        # class_exact dropped from the predictions side, so the label can only
        # come from the classification file.
        self.write_predictions(
            [{key: value for key, value in row.items() if key != "class_exact"} for row in ROWS]
        )
        self.plot("--classification", path, "--label-column", "class_exact")
        self.assertIn(sentinel, self.read_output())

    def test_a_label_in_neither_file_names_both(self):
        path = self.write_classification(self.classification_rows())
        with self.assertRaises(SystemExit) as caught:
            self.plot("--classification", path, "--label-column", "nope")
        message = str(caught.exception)
        self.assertIn("predictions.csv", message)
        self.assertIn("classification.tsv", message)

    def test_duplicate_tile_ids_do_not_multiply_the_guides(self):
        rows = self.classification_rows()
        path = self.write_classification(rows + rows)
        summary = self.plot("--classification", path)
        self.assertIn("kept 2", summary)
        self.assertIn("2 of 2 plotted guide(s) annotated", summary)

    def test_an_unknown_sequence_column_is_rejected(self):
        """Asking for a column and silently getting nothing would be worse."""
        with self.assertRaises(SystemExit) as caught:
            self.plot("--sequence-column", "nope")
        self.assertIn("nope", str(caught.exception))


class TestRealData(unittest.TestCase):
    """The committed ppib dataset, if it is present. Skips rather than fails."""

    FASTA = os.path.join(_ROOT, "data", "ppib", "NM_000942.5.fasta")
    PREDICTIONS = os.path.join(_ROOT, "results", "ppib", "ppib-predictions.csv")

    def setUp(self):
        for path in (self.FASTA, self.PREDICTIONS):
            if not os.path.exists(path):
                self.skipTest(f"{path} not generated")
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def test_every_guide_lies_within_the_transcript(self):
        """The bounds check over real coordinates: tile_id spans must fit the FASTA."""
        accession, length = plot_guides.read_transcript(self.FASTA)
        rows = plot_guides.read_predictions(self.PREDICTIONS, accession)
        built = guides_from_rows(rows.to_dict("records"), transcript_id=accession)
        check_within_transcript(built, length)
        self.assertTrue(built)

    def test_the_cli_runs_over_the_committed_dataset(self):
        output = os.path.join(self.tmp.name, "guides.html")
        stderr = io.StringIO()
        old = sys.stderr
        sys.stderr = stderr
        try:
            self.assertEqual(
                plot_guides.main(
                    [
                        "--fasta", self.FASTA,
                        "--predictions", self.PREDICTIONS,
                        "--min-efficacy", "0.95",
                        "--output", output,
                    ]
                ),
                0,
            )
        finally:
            sys.stderr = old
        self.assertIn("NM_000942.5", stderr.getvalue())
        self.assertTrue(os.path.exists(output))


class TestPlotlyAssumptions(unittest.TestCase):
    """Pin the Plotly semantics this module relies on, so an upgrade fails loudly.

    The same role ``test_biopython_helpers_behave_as_this_module_assumes`` plays
    for Biopython in ``test_tile.py``.
    """

    def test_horizontal_bar_base_offsets_the_start_of_the_bar(self):
        """guide_track draws spans as base=start-0.5 with x=length, not x=end."""
        import plotly.graph_objects as go

        bar = go.Bar(x=[30], base=[0.5], orientation="h", y=[0])
        self.assertEqual(list(bar.base), [0.5])
        self.assertEqual(list(bar.x), [30])

    def test_the_module_exposes_graph_objects_under_the_expected_name(self):
        self.assertTrue(hasattr(visualize_mod.go, "Bar"))
        self.assertTrue(hasattr(visualize_mod.go, "Figure"))


if __name__ == "__main__":
    unittest.main()
