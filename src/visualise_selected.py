"""Plot only the *selected* guides along a transcript, each over its full length.

Companion to ``visualise_scatterlines``: the exon and CDS lanes, the siRNA bands,
the colours and the hover content are the same, but instead of one point per
tile at its start coordinate, each selected guide is a horizontal segment from
its first to its last base, drawn at the height of its predicted efficacy.

Coordinates are 1-based inclusive, as everywhere in this pipeline, so a guide
covers ``end - start + 1`` bases and its segment runs from ``start - 0.5`` to
``end + 0.5``.
"""

from __future__ import annotations

import sys
from typing import Iterable, Sequence

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from annotation import Cds, Exon
from visualise_scatterlines import (
    CDS_FILL,
    CDS_LANE_HEIGHT,
    DEFAULT_THRESHOLD,
    EXON_FILL,
    EXON_LANE_GAP,
    EXON_LANE_HEIGHT,
    GRID,
    HIGHLIGHT_LABELS,
    INK,
    LANE_EDGE,
    LANE_GAP,
    SIRNA_BAND,
    SURFACE,
    PsPCas13bGuide,
    highlight_category,
    siRNA,
)

GUIDE_FILL = "#2f6fb0"      # every guide, one blue that contrasts with the grey lanes and red bands
TARGET_FILL = "#2f8f6f"     # the selected target region lane
CONSERVED_FILL = "#8a6fb5"  # the conserved-region lane


def read_bed_regions(path: str) -> dict[str, list[tuple[int, int]]]:
    """Regions per reference from a BED file, as 1-based inclusive spans.

    BED is 0-based half-open, so the start is shifted by one. Split on any
    whitespace: these files mix tabs and spaces and may lack a final newline.
    """
    regions: dict[str, list[tuple[int, int]]] = {}
    with open(path) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 3:
                regions.setdefault(fields[0], []).append((int(fields[1]) + 1, int(fields[2])))
    return regions


def read_guide_list(path: str) -> list[str]:
    """Tile ids from a one-per-line list; blank lines and ``#`` comments skipped."""
    with open(path) as handle:
        lines = (line.strip() for line in handle)
        return [line for line in lines if line and not line.startswith("#")]


def nearest_sirna(
    guide: PsPCas13bGuide, sirnas: Sequence[siRNA]
) -> tuple[siRNA | None, str | None, int | None]:
    """The siRNA closest to `guide` (an overlap is distance 0), with its class and gap."""
    from visualise_scatterlines import classify_guide, signed_distance

    if not sirnas:
        return None, None, None
    best = min(sirnas, key=lambda s: abs(signed_distance(guide, s)))
    return best, classify_guide(guide, best), signed_distance(guide, best)


def plot_selected_guides(
    guides: Iterable[PsPCas13bGuide],
    sirnas: Iterable[siRNA],
    length: int,
    output: str,
    inline_plotlyjs: bool = True,
    exons: Sequence[Exon] = (),
    cds: Cds | None = None,
    threshold: float = DEFAULT_THRESHOLD,
    target_region: tuple[int, int] | None = None,
    conserved: Sequence[tuple[int, int]] = (),
    gene: str | None = None,
) -> str | None:
    """One transcript per plot: the selected guides as full-length segments."""
    guides = sorted(guides, key=lambda g: (g.start, g.end))
    if not guides:
        print("skipping: no guides to plot", file=sys.stderr)
        return None

    sirnas = sorted(sirnas, key=lambda s: s.start)
    transcript = guides[0].transcript

    # Same lane layout as plot_scatterlines: exons split by parity so junctions
    # show as steps, then the CDS, then the guide panel.
    lane_specs: list[tuple[str, Sequence[Exon] | None]] = []
    if exons:
        lane_specs.append(("exons", [e for e in exons if e.index % 2]))
        lane_specs.append(("exons", [e for e in exons if not e.index % 2]))
    if cds is not None:
        lane_specs.append(("CDS", None))
    if conserved:
        lane_specs.append(("Conserved", list(conserved)))
    if target_region is not None:
        lane_specs.append(("Target region", [target_region]))
    main_row = len(lane_specs) + 1

    lane_domains: list[tuple[float, float]] = []
    top = 1.0
    for index, (label, _) in enumerate(lane_specs):
        height = EXON_LANE_HEIGHT if label == "exons" else CDS_LANE_HEIGHT
        lane_domains.append((top - height, top))
        following = lane_specs[index + 1][0] if index + 1 < len(lane_specs) else None
        top -= height + (EXON_LANE_GAP if label == following == "exons" else LANE_GAP)

    fig = make_subplots(rows=main_row, cols=1, shared_xaxes=True, vertical_spacing=0)

    states: dict[str, int] = {}
    for guide in guides:
        state = highlight_category(guide, threshold)
        states[state] = states.get(state, 0) + 1
        sirna, cls, gap = nearest_sirna(guide, sirnas)
        overlap = cls == "overlap"
        label = "Overlap with siRNA" if overlap else "Guide"

        # One marker per base along the segment: a line only answers hover at
        # its vertices, so this makes the whole guide hoverable, not its ends.
        xs = [guide.start - 0.5] + list(range(guide.start, guide.end + 1)) + [guide.end + 0.5]
        ys = [guide.score] * len(xs)

        legend_name = HIGHLIGHT_LABELS[state]
        fig.add_trace(go.Scatter(
            x=xs, y=ys,
            mode="lines+markers",
            name=legend_name,
            showlegend=False,
            line=dict(color=GUIDE_FILL, width=5),
            marker=dict(color=GUIDE_FILL, size=5, opacity=0),
            customdata=[(
                guide.end, guide.tile_id, guide.sequence,
                sirna.qseqid if sirna else "-", gap if gap is not None else "-",
            )] * len(xs),
            hovertemplate=(
                f"<b>{label}</b><br>"
                "efficacy %{y:.4f}<br>"
                f"{guide.start:,d}-{guide.end:,d}<br>"
                "%{customdata[1]}<br>"
                "nearest siRNA %{customdata[3]}, gap %{customdata[4]} nt<br>"
                "5'-%{customdata[2]}-3'<extra></extra>"
            ),
        ), row=main_row, col=1)

    for sirna in sirnas:
        fig.add_vrect(
            x0=sirna.start - 0.5, x1=sirna.end + 0.5,
            fillcolor=SIRNA_BAND, opacity=0.15, line_width=0, layer="below",
            annotation_text=sirna.qseqid, annotation_position="top left",
            annotation_font_color=INK,
            row=main_row, col=1,
        )

    for row, (label, features) in enumerate(lane_specs, start=1):
        if label in ("Target region", "Conserved"):
            fill = TARGET_FILL if label == "Target region" else CONSERVED_FILL
            fig.add_trace(go.Bar(
                x=[e - s + 1 for s, e in features],
                base=[s - 0.5 for s, _ in features],
                y=[label] * len(features), orientation="h",
                marker=dict(color=fill, line=dict(color=LANE_EDGE, width=1)),
                showlegend=False,
                customdata=[(s, e, e - s + 1) for s, e in features],
                hovertemplate=(
                    f"<b>{label}</b><br>"
                    "%{customdata[0]:,d}-%{customdata[1]:,d}<br>"
                    "%{customdata[2]:,d} nt<extra></extra>"
                ),
            ), row=row, col=1)
            continue
        if label == "CDS":
            fig.add_trace(go.Bar(
                x=[cds.length], base=[cds.start - 0.5], y=[label], orientation="h",
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
        fig.add_trace(go.Bar(
            x=[e.length for e in features],
            base=[e.start - 0.5 for e in features],
            y=[label] * len(features), orientation="h",
            marker=dict(color=EXON_FILL, line=dict(color=LANE_EDGE, width=1)),
            showlegend=False,
            customdata=[(e.index, e.start, e.end, e.length) for e in features],
            hovertemplate=(
                "<b>exon %{customdata[0]}</b><br>"
                "%{customdata[1]:,d}-%{customdata[2]:,d}<br>"
                "%{customdata[3]:,d} nt<extra></extra>"
            ),
        ), row=row, col=1)

    name = f"{transcript} ({gene})" if gene else transcript
    region = (
        f", target region {target_region[0]:,}-{target_region[1]:,}"
        if target_region else ""
    )
    fig.update_layout(
        title=(
            f"Selected PsPCas13b guides along {name}"
            f"<br><span style='font-size:11px'>{len(guides)} guides{region}</span>"
        ),
        bargap=0,
        plot_bgcolor=SURFACE, paper_bgcolor=SURFACE,
        font=dict(color=INK),
        hovermode="closest",
        showlegend=False,
    )
    fig.update_xaxes(tickformat=",d", hoverformat=",d", range=[0.5, length + 0.5])
    fig.update_xaxes(gridcolor=GRID, title_text="Position", row=main_row, col=1)
    fig.update_yaxes(
        gridcolor=GRID, zeroline=False, title_text="Predicted efficacy",
        domain=[0.0, top], row=main_row, col=1,
    )
    for index, (label, _) in enumerate(lane_specs, start=1):
        fig.update_yaxes(
            showgrid=False, zeroline=False, ticks="",
            domain=list(lane_domains[index - 1]),
            showticklabels=label != "exons",
            row=index, col=1,
        )
    if exons:
        exon_domains = [d for d, (label, _) in zip(lane_domains, lane_specs) if label == "exons"]
        fig.add_annotation(
            text="exons", x=0, xref="paper", xanchor="right", xshift=-8,
            y=(min(d[0] for d in exon_domains) + max(d[1] for d in exon_domains)) / 2,
            yref="paper", yanchor="middle", showarrow=False,
            font=dict(color=INK, size=11),
        )

    fig.write_html(output, include_plotlyjs=True if inline_plotlyjs else "cdn")
    return output
