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
from typing import Iterable, NamedTuple

import pandas as pd
import plotly.graph_objects as go


SURFACE = "#fcfcfb"      # plot background, and the ring that separates stacked points
GRID = "#e8e7e2"
INK = "#52514e"          # axis and tick text
LINE = "#d9d8d3"         # the score profile joining the points, kept recessive
SIRNA_BAND = "#d6453f"   # the siRNA footprint, light red at low opacity
OVERLAP_RING = "#2f2e2c" # the highlight on guides overlapping an siRNA

# Marker fill, by predicted efficacy: green at 0, yellow at 1. Lightness rises
# with the score, so the best guides are also the brightest points. It stops at
# a saturated yellow rather than running to near-white, which would make the
# highest-scoring guides -- the ones worth finding -- vanish into SURFACE.
EFFICACY_SCALE = [
    [0.0, "#14532d"],
    [0.35, "#2f8a4a"],
    [0.70, "#8ac926"],
    [1.0, "#f4d03f"],
]


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
) -> str | None:
    """One transcript per plot: the score profile with the siRNA bands on it.

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
            }
            for g in guides
        ]
    )

    fig = go.Figure()
    # The profile first and in one piece: colouring it by class would break it
    # into disconnected segments, and it is the points that carry the class.
    fig.add_trace(go.Scattergl(
        x=df["start"], y=df["score"],
        mode="lines",
        line=dict(color=LINE, width=1),
        hoverinfo="skip",       # the markers' own hover covers every point
        showlegend=False,
        name="score profile",
    ))

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
            showlegend=label not in shown,
            marker=dict(
                # Fill is the score, on one shared coloraxis so every trace
                # and every plot reads against the same colourbar.
                color=points["score"],
                coloraxis="coloraxis",
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
        ))
        shown.add(label)

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
        )

    fig.update_layout(
        title=(
            f"PsPCas13b guide efficacy along {transcript} "
            f"({len(sirnas)} siRNA site{'' if len(sirnas) == 1 else 's'})"
        ),
        xaxis_title="Position",
        yaxis_title="Predicted efficacy",
        legend_title_text="Position vs nearest siRNA",
        # Pinned to 0-1 rather than to the data: efficacy is a probability, and
        # a per-plot range would give the same colour different meanings on
        # each transcript, so two plots could not be compared.
        coloraxis=dict(
            colorscale=EFFICACY_SCALE,
            cmin=0.0,
            cmax=1.0,
            colorbar=dict(
                title="Predicted<br>efficacy",
                tickformat=".1f",
                thickness=14,
                len=0.5,
                yanchor="bottom",
                y=0,
            ),
        ),
        plot_bgcolor=SURFACE,
        paper_bgcolor=SURFACE,
        font=dict(color=INK),
        hovermode="closest",
    )
    fig.update_yaxes(gridcolor=GRID, zeroline=False)
    # Transcript coordinates, so show them in full: Plotly's default SI
    # formatting turns 2724 into "2.7k", which is useless for finding a guide.
    fig.update_xaxes(
        tickformat=",d",    # tick labels: 2,724
        hoverformat=",d",   # the x value in the hover box, likewise
        gridcolor=GRID,
        # Always the whole transcript, never just the span the guides occupy:
        # uncovered regions must read as genuinely uncovered.
        range=[0.5, lengths.get(transcript, max(df["end"])) + 0.5],
    )
    # Inlining plotly.js costs ~3 MB but opens without a network connection,
    # which is what an archived or emailed plot needs; "cdn" keeps it small.
    fig.write_html(output, include_plotlyjs=True if inline_plotlyjs else "cdn")
    return output
