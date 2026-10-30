import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))

from visualise_scatterlines import PsPCas13bGuide, siRNA  # noqa: E402
from visualise_selected import nearest_sirna, plot_selected_guides, read_guide_list  # noqa: E402


def guide(start, end, score=0.95, seq="GGAAAA"):
    return PsPCas13bGuide("T1:%d-%d" % (start, end), "T1", seq, score, start, end)


class TestSelectedPlot(unittest.TestCase):
    def test_read_guide_list_skips_blanks_and_comments(self):
        with tempfile.NamedTemporaryFile("w", suffix=".list", delete=False) as f:
            f.write("# header\nA:1-30\n\nB:5-34")
        try:
            self.assertEqual(read_guide_list(f.name), ["A:1-30", "B:5-34"])
        finally:
            os.unlink(f.name)

    def test_nearest_sirna_prefers_overlap_then_closest(self):
        near, far = siRNA("near", "T1", 100, 120), siRNA("far", "T1", 500, 520)
        self.assertEqual(nearest_sirna(guide(110, 139), [far, near]), (near, "overlap", 0))
        self.assertEqual(nearest_sirna(guide(300, 329), [far, near])[:2], (far, "upstream"))
        self.assertEqual(nearest_sirna(guide(1, 30), []), (None, None, None))

    def test_plot_writes_one_trace_per_guide_spanning_its_length(self):
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "p.html")
            plot_selected_guides(
                [guide(10, 39), guide(200, 229, score=0.5)],
                [siRNA("s", "T1", 20, 40)], 300, out, target_region=(5, 250),
            )
            html = open(out).read()
        self.assertIn("nearest siRNA", html)
        self.assertIn("T1:10-39", html)
        self.assertIn("Target region", html)


if __name__ == "__main__":
    unittest.main()
