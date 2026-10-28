"""Plot PsPCas13b guides along a transcript, positioned relative to known siRNAs.

Each guide is classified as ``upstream``, ``downstream`` or ``overlap`` **with
respect to one siRNA**. A transcript carries several siRNAs, so a guide can be
upstream of one and downstream of another: the classification is therefore a
property of a (guide, siRNA) *pair*, not of a guide, and this module emits one
record per pair. Guides on a transcript with no siRNA are kept with
``classification=None`` so the guide set stays complete.

Conventions follow src/tile.py: coordinates are 1-based inclusive, so a guide
spans ``end - start + 1`` bases and ``sequence[start - 1:end]`` is the match.
Sequences are 5' -> 3'. ``spacer_seq`` is the guide, never the target.

This is the library half of the stage: it parses nothing from the command line
and reads no file by a hard-coded path. Every function takes a DataFrame or a
sequence of tuples, so the classification rules are testable without plotly and
the layout without pandas. ``scripts/plot_scatterlines.py`` is the CLI around
it -- argument parsing, file reading, the loop over transcripts and the stderr
summary live there.
"""

from __future__ import annotations

import sys
from typing import Iterable, NamedTuple, Sequence

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from annotation import Cds, Exon


SURFACE = "#fcfcfb"      # plot background, and the ring that separates stacked points
GRID = "#e8e7e2"
INK = "#52514e"          # axis and tick text
LINE = "#d9d8d3"         # the score profile joining the points, kept recessive
SIRNA_BAND = "#d6453f"   # the siRNA footprint, light red at low opacity
OVERLAP_RING = "#2f2e2c" # the highlight on guides overlapping an siRNA

# The annotation lanes, all one light grey. On an mRNA the exons abut, so a
# single row of one fill would draw an unbroken bar from base 1 to the end of
# the transcript with no boundary visible anywhere. The exons are therefore
# split across **two** rows by parity -- odd above, even below -- which makes
# every junction a vertical step between the rows. The offset does the work a
# colour change would otherwise have to do, so both rows can be the same grey.
EXON_FILL = "#b9b7b1"
CDS_FILL = "#6f6d68"     # darker, so the CDS reads as a different kind of feature
LANE_EDGE = SURFACE      # hairline in the background colour, around each bar

# Lane geometry, as fractions of the figure height. The lanes are deliberately
# thin: they are a positional reference for the scatter below, not the subject
# of the figure.
#
# The gaps are set per boundary rather than by `make_subplots`' uniform
# `vertical_spacing`, which cannot do what is wanted here: the two exon rows are
# one track split in two and should read as a pair, so they are separated by a
# hairline, while the scatter below needs a real gap to read as a different
# panel. One spacing value would have to compromise between the two.
EXON_LANE_HEIGHT = 0.022
CDS_LANE_HEIGHT = 0.022
EXON_LANE_GAP = 0.004    # between the two exon rows: just enough to see the step
LANE_GAP = 0.022         # between tracks, and above the scatter

# Marker fill, by whether a guide is worth looking at rather than by its score.
#
# A continuous scale over the efficacy was replaced by these three states
# deliberately: the score distribution is top-heavy -- 79% of guides exceed 0.7
# -- so a gradient spends most of its range on distinctions that do not change
# a decision, while the threshold does. The exact score is still in the hover.
#
# Orange is a **refinement of yellow, not a competitor**: a guide is orange only
# if it clears the threshold *and* begins with GG, so the orange points are
# always a subset of the ones worth considering. A GG guide below the threshold
# stays grey, because its GG start is no use if the guide itself is poor.
BELOW_FILL = "#b0aea8"   # below threshold: present, recessive, not competing
ABOVE_FILL = "#f4d03f"   # above threshold
GG_FILL = "#e8833a"      # above threshold and starting with GG

HIGHLIGHT_LABELS = {
    "below": "Below threshold",
    "above": "Above threshold",
    "gg": "Above threshold, GG start",
}
HIGHLIGHT_FILLS = {"below": BELOW_FILL, "above": ABOVE_FILL, "gg": GG_FILL}

#: Only a default -- the caller is expected to set this. 0.9 rather than
#: something lower because the distribution is top-heavy: at 0.7 nearly four
#: guides in five clear the bar, which highlights almost the whole transcript
#: and so highlights nothing.
DEFAULT_THRESHOLD = 0.9


# label, outline colour, marker symbol. Every guide is a circle, so position
# and the efficacy fill are the only things the shape competes with; the
# highlight on the overlaps is the outline, plus a larger marker.
#
# `upstream`, `downstream` and the no-siRNA case share the label "Guide": they
# are styled identically, and two legend entries that look the same would
# imply a distinction the plot does not draw. Which side of an siRNA a guide
# misses by is still in the hover box, and `classification` still carries it.
#
# `overlap` is last because it must be drawn last. A guide is emitted once per
# siRNA on its transcript, so its copies stack at the same point, and the trace
# drawn last is the one on top -- an overlap listed first would be buried under
# the plain copies of the very same guide.
CLASS_STYLE = {
    "upstream": ("Guide", SURFACE, "circle"),
    "downstream": ("Guide", SURFACE, "circle"),
    None: ("Guide", SURFACE, "circle"),
    "overlap": ("Overlap with siRNA", OVERLAP_RING, "circle"),
}


BLASTN_COLUMNS = ("qseqid", "sseqid", "sstart", "send", "bitscore")
PREDICTION_COLUMNS = (
    "tile_id",
    "transcript_id",
    "start",
    "end",
    "spacer_seq",
    "predicted_efficacy",
)


class siRNA(NamedTuple):
    qseqid: str                 # siRNA name, from the blastn query
    transcript: str             # accession it mapped to (blastn subject)
    start: int                  # 1-based inclusive, always <= end
    end: int


class PsPCas13bGuide(NamedTuple):
    """One guide, judged against one siRNA or against none.

    `classification` is relative to `sirna`, not an intrinsic property of the
    guide, so a guide on a transcript with several siRNAs appears once per
    siRNA with the same tile_id. The siRNA is carried whole rather than by
    name because the plot needs its span to draw the band.
    """

    tile_id: str
    transcript: str
    sequence: str               # the spacer, 5' -> 3'
    score: float
    start: int                  # 1-based inclusive
    end: int
    sirna: siRNA | None = None
    classification: str | None = None
    distance: int | None = None  # signed nt gap; negative upstream, 0 on overlap



def classify_guide(guide: PsPCas13bGuide, sirna: siRNA) -> str:
    """Where `guide` sits relative to `sirna`. Spans are 1-based inclusive."""
    if guide.end < sirna.start:
        return "upstream"
    elif guide.start > sirna.end:
        return "downstream"
    else:
        return "overlap"


def signed_distance(guide: PsPCas13bGuide, sirna: siRNA) -> int:
    """Gap in nt between guide and siRNA: negative upstream, positive
    downstream, 0 when they overlap."""
    if guide.end < sirna.start:
        return guide.end - sirna.start
    if guide.start > sirna.end:
        return guide.start - sirna.end
    return 0


def build_sirnas(df: pd.DataFrame) -> list[siRNA]:
    """One siRNA per blastn query, taking its highest-scoring HSP."""

    # Best HSP per siRNA; idxmax breaks ties by first occurrence.
    df_filt = df.loc[df.groupby("qseqid")["bitscore"].idxmax()]

    sirnas = []
    for _, row in df_filt.iterrows():
        # blastn reports sstart > send for a minus-strand hit; the span is the
        # same either way, and classification only needs it ordered.
        start, end = int(row["sstart"]), int(row["send"])
        sirnas.append(
            siRNA(
                qseqid=str(row["qseqid"]),
                transcript=str(row["sseqid"]),
                start=min(start, end),
                end=max(start, end),
            )
        )
    return sirnas


def build_guides(
    df: pd.DataFrame, sirnas: Iterable[siRNA] = ()
) -> list[PsPCas13bGuide]:
    """One guide per (prediction row, siRNA-on-the-same-transcript) pair.

    A row whose transcript carries no siRNA yields a single guide with
    `sirna`, `classification` and `distance` all None, so no guide is
    silently dropped.
    """

    by_transcript: dict[str, list[siRNA]] = {}
    for sirna in sirnas:
        by_transcript.setdefault(sirna.transcript, []).append(sirna)

    guides = []
    for _, row in df.iterrows():
        guide = PsPCas13bGuide(
            tile_id=str(row["tile_id"]),
            transcript=str(row["transcript_id"]),
            sequence=str(row["spacer_seq"]),
            score=float(row["predicted_efficacy"]),
            start=int(row["start"]),
            end=int(row["end"]),
        )
        matches = by_transcript.get(guide.transcript, ())
        if not matches:
            guides.append(guide)
            continue
        for sirna in matches:
            guides.append(
                guide._replace(
                    sirna=sirna,
                    classification=classify_guide(guide, sirna),
                    distance=signed_distance(guide, sirna),
                )
            )
    return guides


def highlight_category(guide: PsPCas13bGuide, threshold: float) -> str:
    """``below``, ``above`` or ``gg`` -- which of the three fills a guide gets.

    The threshold is **strict**: a guide exactly on it counts as below, so
    "above 0.9" never includes 0.9 itself.

    ``gg`` requires both conditions, so it is a subset of the above-threshold
    guides rather than a third, competing state -- a GG start is a property of
    how easily the guide can be transcribed, which is only worth seeing on a
    guide that is otherwise worth using. The spacer is RNA, but G is G in both
    alphabets, so no transcription step is needed to test it.
    """
    if guide.score <= threshold:
        return "below"
    return "gg" if guide.sequence.upper().startswith("GG") else "above"


def transcript_lengths(df: pd.DataFrame) -> dict[str, int]:
    """Transcript length from the *unfiltered* tiling.

    Stage 1 tiles at step 1 across the whole transcript, so the largest tile
    end equals the transcript length. Taken before any efficacy filter,
    because a filtered set would shrink the axis and make uncovered regions
    read as absent rather than uncovered.
    """
    return df.groupby("transcript_id")["end"].max().astype(int).to_dict()


def plot_scatterlines(
    guides: Iterable[PsPCas13bGuide],
    sirnas: Iterable[siRNA],
    lengths: dict[str, int],
    output: str,
    inline_plotlyjs: bool = True,
    exons: Sequence[Exon] = (),
    cds: Cds | None = None,
    threshold: float = DEFAULT_THRESHOLD,
) -> str | None:
    """One transcript per plot: the score profile with the siRNA bands on it.

    `exons` and `cds` are optional. Each one present adds a thin lane above the
    scatter, sharing its x-axis so a feature sits directly over the guides that
    target it -- the lanes are a positional reference, not the subject. Passing
    neither gives a single-panel figure, as this function gave before the lanes
    existed.
    """
    # Every record, not one per guide: `build_guides` emits a guide once per
    # siRNA on its transcript, so a guide on a transcript with three siRNAs
    # is drawn three times at the same position and score, once in each
    # class. Points therefore stack, and a class count is a count of
    # guide-siRNA pairs rather than of guides.
    guides = sorted(guides, key=lambda g: g.start)

    if not guides:
        print("skipping: no guides to plot", file=sys.stderr)
        return None

    sirnas = sorted(sirnas, key=lambda s: s.start)
    transcript = guides[0].transcript

    df = pd.DataFrame(
        [
            {
                "start": g.start,
                "end": g.end,
                "score": g.score,
                "classification": g.classification,
                "tile_id": g.tile_id,
                "sequence": g.sequence,
                "sirna": g.sirna.qseqid if g.sirna else "-",
                "distance": g.distance if g.distance is not None else "-",
                "highlight": highlight_category(g, threshold),
            }
            for g in guides
        ]
    )

    # One row per annotation lane, then the scatter, sharing a single x-axis so
    # the lanes and the guides cannot drift apart under zoom or panning. The
    # scatter keeps whatever height the lanes leave it.
    # The lane rows, top to bottom: the odd exons, the even exons, then the CDS.
    # Splitting the exons by parity is what makes their boundaries visible --
    # see EXON_FILL above.
    lane_specs: list[tuple[str, Sequence[Exon] | None]] = []
    if exons:
        lane_specs.append(("exons", [e for e in exons if e.index % 2]))
        lane_specs.append(("exons", [e for e in exons if not e.index % 2]))
    if cds is not None:
        lane_specs.append(("CDS", None))

    scatter_row = len(lane_specs) + 1

    # Each row's vertical extent, top down, in paper units. Computed here and
    # applied to the axes below, because the gap above a row depends on what
    # that row and the one before it are -- see LANE_GAP above.
    lane_domains: list[tuple[float, float]] = []
    top = 1.0
    for index, (label, _) in enumerate(lane_specs):
        height = EXON_LANE_HEIGHT if label == "exons" else CDS_LANE_HEIGHT
        lane_domains.append((top - height, top))
        following = lane_specs[index + 1][0] if index + 1 < len(lane_specs) else None
        top -= height + (EXON_LANE_GAP if label == following == "exons" else LANE_GAP)
    scatter_domain = (0.0, top)

    fig = make_subplots(
        rows=scatter_row,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0,   # every domain is set explicitly below
    )

    # The profile first and in one piece: colouring it by class would break it
    # into disconnected segments, and it is the points that carry the class.
    fig.add_trace(go.Scattergl(
        x=df["start"], y=df["score"],
        mode="lines",
        line=dict(color=LINE, width=1),
        hoverinfo="skip",       # the markers' own hover covers every point
        showlegend=False,
        name="score profile",
    ), row=scatter_row, col=1)

    # One trace per class, in CLASS_STYLE order rather than the order the
    # classes happen to appear in, so the legend reads the same on every plot
    # and the overlaps land on top. A class with no guides is left out, and the
    # three classes sharing a label contribute one legend entry between them.
    shown: set[str] = set()
    for cls, (label, colour, symbol) in CLASS_STYLE.items():
        points = df[df["classification"].isna() if cls is None
                    else df["classification"] == cls]
        if points.empty:
            continue
        highlight = cls == "overlap"
        fig.add_trace(go.Scattergl(
            x=points["start"], y=points["score"],
            mode="markers",
            name=label,
            legendgroup=label,
            legendgrouptitle_text="Position vs nearest siRNA",
            showlegend=label not in shown,
            marker=dict(
                # Fill is the highlight state, per point -- not a coloraxis.
                # Drawing order still comes from CLASS_STYLE, so an overlapping
                # guide stays on top of its own plain copies whatever its fill.
                color=[HIGHLIGHT_FILLS[state] for state in points["highlight"]],
                symbol=symbol,
                # The overlaps are bigger as well as ringed: a contrasting
                # outline alone is easy to lose among 21,000 plain points.
                size=11 if highlight else 5,
                # A ring in SURFACE keeps touching points separate -- at a 1 nt
                # step they sit one base apart. On the overlaps it is instead a
                # heavy charcoal: dark enough to read against the light end of
                # the green-to-yellow fill, where the good guides are, and
                # achromatic so it cannot be mistaken for a fill value.
                line=dict(color=colour, width=2.0 if highlight else 0.6),
            ),
            customdata=points[["end", "tile_id", "sequence", "sirna", "distance"]],
            hovertemplate=(
                f"<b>{label}</b><br>"
                "efficacy %{y:.4f}<br>"
                "%{x:,d}-%{customdata[0]:,d}<br>"
                "%{customdata[1]}<br>"
                "nearest siRNA %{customdata[3]}, gap %{customdata[4]} nt<br>"
                # Ends marked: every sequence in this pipeline is 5' -> 3'.
                "5'-%{customdata[2]}-3'<extra></extra>"
            ),
        ), row=scatter_row, col=1)
        shown.add(label)

    # Counts for the legend and subtitle are over **distinct guides**, not over
    # the rows plotted. `build_guides` emits a guide once per siRNA, so on a
    # transcript with three siRNAs every guide is drawn three times and a count
    # of plotted rows would treble -- reporting "3,870 above threshold" for
    # 1,290 guides, which no other output in this pipeline would agree with.
    # The points still stack; only the arithmetic is deduplicated.
    counts = df.drop_duplicates("tile_id")["highlight"].value_counts()
    n_guides = int(counts.sum())

    # The fill legend, as entries with no points: the fills vary *within* each
    # class trace, so no real trace can stand for one. Only the states actually
    # present are listed, and in a fixed order, so a plot with no GG guide
    # above the threshold says so by omission rather than showing a key to a
    # colour that never appears.
    present = set(df["highlight"])
    for state in ("below", "above", "gg"):
        if state not in present:
            continue
        fig.add_trace(go.Scattergl(
            x=[None], y=[None],
            mode="markers",
            name=f"{HIGHLIGHT_LABELS[state]} ({int(counts.get(state, 0)):,})",
            legendgroup="fill",
            legendgrouptitle_text=f"Fill (threshold {threshold:g})",
            marker=dict(
                color=HIGHLIGHT_FILLS[state], size=9,
                line=dict(color=SURFACE, width=0.6),
            ),
            hoverinfo="skip",
        ), row=scatter_row, col=1)

    # One band per siRNA. Light red, and the only hue in the figure besides the
    # overlap points it explains, so a band reads as the siRNA footprint rather
    # than as another guide class. A transcript with no siRNA gets no band.
    for sirna in sirnas:
        fig.add_vrect(
            x0=sirna.start - 0.5, x1=sirna.end + 0.5,  # whole bases, not centres
            fillcolor=SIRNA_BAND,
            opacity=0.15,                 # lower = more see-through
            line_width=0,                 # removes the border
            layer="below",                # keeps the band behind the data
            annotation_text=sirna.qseqid,
            annotation_position="top left",
            annotation_font_color=INK,
            row=scatter_row, col=1,
        )

    # The lanes, drawn as horizontal bars on a one-category axis so each
    # feature keeps its own hover: a shape would be silent, and the exon number
    # is most of the point of having the lane.
    for row, (label, features) in enumerate(lane_specs, start=1):
        if label == "CDS":
            fig.add_trace(go.Bar(
                x=[cds.length],
                base=[cds.start - 0.5],
                y=[label],
                orientation="h",
                marker=dict(color=CDS_FILL, line=dict(width=0)),
                showlegend=False,
                customdata=[(cds.start, cds.end, cds.length)],
                hovertemplate=(
                    "<b>CDS</b><br>"
                    "%{customdata[0]:,d}-%{customdata[1]:,d}<br>"
                    "%{customdata[2]:,d} nt<extra></extra>"
                ),
            ), row=row, col=1)
            continue

        # An exon row can legitimately be empty -- a transcript with a single
        # exon leaves the even row with nothing in it -- and an empty Bar is
        # fine: the row still reserves its height, so the lanes stay aligned.
        fig.add_trace(go.Bar(
            x=[exon.length for exon in features],
            # `base` is where each bar starts; x is its width. Half a base is
            # subtracted so a bar covers whole bases rather than centres, as
            # the siRNA bands do.
            base=[exon.start - 0.5 for exon in features],
            y=[label] * len(features),
            orientation="h",
            marker=dict(color=EXON_FILL, line=dict(color=LANE_EDGE, width=1)),
            showlegend=False,
            customdata=[(e.index, e.start, e.end, e.length) for e in features],
            hovertemplate=(
                "<b>exon %{customdata[0]}</b><br>"
                "%{customdata[1]:,d}-%{customdata[2]:,d}<br>"
                "%{customdata[3]:,d} nt<extra></extra>"
            ),
        ), row=row, col=1)

    fig.update_layout(
        title=(
            f"PsPCas13b guide efficacy along {transcript}"
            # The siRNA and exon counts that used to sit in brackets here are
            # still on the plot: the siRNAs in the legend and as labelled
            # bands, the exons in the lanes and their hover. The run summary
            # prints both.
            f"<br><span style='font-size:11px'>"
            f"{n_guides - int(counts.get('below', 0)):,} of {n_guides:,} "
            f"guides above {threshold:g}, of which "
            f"{int(counts.get('gg', 0)):,} start with GG</span>"
        ),
        # The lane bars are the full height of their lane. Without this they
        # would be drawn at a bar chart's default 80%, leaving a gap that reads
        # as if the annotation stopped short.
        bargap=0,
        plot_bgcolor=SURFACE,
        paper_bgcolor=SURFACE,
        font=dict(color=INK),
        hovermode="closest",
    )
    # Transcript coordinates, so show them in full: Plotly's default SI
    # formatting turns 2724 into "2.7k", which is useless for finding a guide.
    #
    # Applied to every row, and the range in particular: the lanes share the
    # scatter's x-axis, but `shared_xaxes` only links the *interaction*, so
    # without an explicit range on each row a lane could be autoscaled to its
    # own features and a feature would then sit over the wrong guides.
    axis_range = [0.5, lengths.get(transcript, max(df["end"])) + 0.5]
    fig.update_xaxes(
        tickformat=",d",    # tick labels: 2,724
        hoverformat=",d",   # the x value in the hover box, likewise
        # Always the whole transcript, never just the span the guides occupy:
        # uncovered regions must read as genuinely uncovered.
        range=axis_range,
    )
    fig.update_xaxes(gridcolor=GRID, row=scatter_row, col=1)
    fig.update_xaxes(title_text="Position", row=scatter_row, col=1)

    fig.update_yaxes(
        gridcolor=GRID, zeroline=False,
        title_text="Predicted efficacy",
        row=scatter_row, col=1,
    )
    # No grid and no zero line on a lane: a gridline across a one-category axis
    # reads as a feature boundary. Each lane's domain is set here rather than
    # left to `vertical_spacing`.
    for index, (label, _) in enumerate(lane_specs, start=1):
        fig.update_yaxes(
            showgrid=False, zeroline=False, ticks="",
            domain=list(lane_domains[index - 1]),
            # The two exon rows are one track, so they get one label, placed
            # between them as an annotation below. Letting each row draw its
            # own tick would print "exons" twice, and at this lane height the
            # two would collide -- and reading them as two different tracks is
            # exactly the wrong conclusion.
            showticklabels=label != "exons",
            row=index, col=1,
        )
    fig.update_yaxes(domain=list(scatter_domain), row=scatter_row, col=1)

    # The single label for the exon pair, centred on the two rows it spans.
    if exons:
        exon_domains = [
            domain for domain, (label, _) in zip(lane_domains, lane_specs)
            if label == "exons"
        ]
        fig.add_annotation(
            text="exons",
            x=0, xref="paper", xanchor="right", xshift=-8,
            y=(min(d[0] for d in exon_domains) + max(d[1] for d in exon_domains)) / 2,
            yref="paper", yanchor="middle",
            showarrow=False,
            font=dict(color=INK, size=11),
        )
    # Inlining plotly.js costs ~3 MB but opens without a network connection,
    # which is what an archived or emailed plot needs; "cdn" keeps it small.
    fig.write_html(output, include_plotlyjs=True if inline_plotlyjs else "cdn")
    return output
