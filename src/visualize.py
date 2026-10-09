"""Interactive maps of where guides sit along a transcript.

Takes a transcript's length and a set of already-filtered guides and builds a
Plotly figure showing each guide's position and extent against the full span of
the transcript. The point is to answer "where on the message do my surviving
guides actually land, and what did the filter leave uncovered" -- so the x axis
is always the whole transcript, 1 to ``length``, never just the range the guides
happen to occupy.

Coordinate convention
---------------------
Positions are **1-based, inclusive**, the same convention as ``src/tile.py`` and
the ``tile_id`` join key, so a guide spans ``end - start + 1`` bases and
``sequence[start - 1:end]`` is its target. Each base is drawn as a unit cell
centred on its own coordinate: a guide is rendered from ``start - 0.5`` to
``end + 0.5``, so a 30-mer at ``start=1`` occupies 0.5 to 30.5 and abuts, rather
than overlaps, a guide starting at 31.

Filtering happens elsewhere
---------------------------
Nothing here filters, scores or ranks. The caller decides which guides are worth
plotting -- by predicted efficacy, by ``class_exact``, by GC, by whatever stage
comes next -- and this module draws what it is handed. That keeps the plot honest
about being a view of one particular filter, and keeps the filter's rules out of
the plotting code.

Lanes
-----
Filtered guides from step-1 tiling overlap heavily: a run of consecutive
high-scoring tiles is 30 nt wide and 1 nt apart. Drawing those on one row would
stack them into a single opaque block and hide how many there are, so
:func:`pack_lanes` spreads overlapping guides onto separate rows (greedy
interval packing, first lane that has room). Lane number carries no meaning
beyond "these did not fit on one line"; the number of lanes at a given position
is a direct read-out of guide density there.

Plotly
------
Figure construction is Plotly's and is confined to :func:`guide_track` and
:func:`write_html`. :class:`Guide`, :func:`guides_from_rows` and
:func:`pack_lanes` are plain Python and testable without it. No function here
takes a DataFrame: rows come in as mappings, which is what both
``pandas.DataFrame.to_dict("records")`` and ``csv.DictReader`` already produce,
so pandas is not a dependency of this stage.

The ``plotly`` and ``pandas`` packages are in ``environment.yml``; use the
``pik3-guide`` env.
"""

from __future__ import annotations

from typing import Iterable, Mapping, NamedTuple, Sequence

try:
    import plotly.graph_objects as go
except ImportError as exc:  # pragma: no cover - environment problem, not logic
    raise ImportError(
        "plotly is required for src/visualize.py. Install the project env:\n"
        "    mamba env update -f environment.yml\n"
        "    conda activate pik3-guide"
    ) from exc


# Drawn under the lanes, pinned to the bottom of the plot in paper coordinates so
# it is visible in both y modes without competing for data-space room.
_TRANSCRIPT_BAND = (0.0, 0.035)

#: Annotation values longer than this are broken at ``;`` in the hover box. One
#: accession-with-coordinates locus is ~24 characters, so this keeps a single
#: locus on one line and splits a multi-reference string.
_WRAP_OVER = 40

_DEFAULT_COLOR = "#2a6f97"
_TRANSCRIPT_COLOR = "#b8c4cc"


class Guide(NamedTuple):
    """One guide to draw, as a span on a single transcript.

    ``start`` and ``end`` are 1-based inclusive coordinates on the transcript
    (see the module docstring). ``score`` is whatever continuous value the caller
    wants the colour or height to carry -- predicted efficacy in the worked
    example -- and is ``None`` when there is nothing to colour by. ``label`` is
    free text appended to the hover box, for a class or a note.

    ``sequence`` is the guide (spacer) itself, shown in the hover box so a guide
    can be read off the plot and ordered without going back to the CSV. It is the
    spacer written 5' -> 3', per the convention in ``src/tile.py`` -- not the
    target -- and nothing here checks that, since the caller names the column.

    ``annotations`` are ordered ``(name, value)`` pairs appended to the hover box
    one per line, for whatever a later stage knows about the guide -- its family
    classification, a GC fraction, an off-target count. They are deliberately
    untyped and unnamed here: this module should not have to learn the column
    names of every stage that wants to annotate a guide.
    """

    tile_id: str
    start: int
    end: int
    score: float | None = None
    label: str | None = None
    sequence: str | None = None
    annotations: tuple[tuple[str, str], ...] = ()

    @property
    def length(self) -> int:
        return self.end - self.start + 1


def guides_from_rows(
    rows: Iterable[Mapping[str, object]],
    *,
    transcript_id: str | None = None,
    score_column: str | None = "predicted_efficacy",
    label_column: str | None = None,
    sequence_column: str | None = None,
    annotation_columns: Sequence[str] = (),
) -> list[Guide]:
    """Build :class:`Guide` records from already-filtered tabular rows.

    ``rows`` is any iterable of mappings with ``tile_id``, ``start`` and ``end``
    keys -- a ``csv.DictReader``, or ``DataFrame.to_dict("records")`` after the
    caller's own filtering. Coordinates are coerced to ``int``, so the string
    values a ``DictReader`` yields are fine.

    ``transcript_id`` keeps only rows whose ``transcript_id`` column matches,
    which is how a multi-transcript predictions file is reduced to one plot. It
    is a hard error if that column is missing while a filter was requested,
    rather than silently plotting every transcript's guides on one axis.

    ``score_column`` is ignored when the column is absent, so the same call
    works on a predictions file and on a plain tiles file; a present but empty
    value becomes ``None``. Pass ``score_column=None`` to ignore scores entirely.

    ``sequence_column`` names the column holding the guide (spacer) sequence to
    show in the hover box -- ``spacer_seq`` in a predictions file, ``guide_rna``
    in a stage-1 tiles file. It is not defaulted, because guessing wrong here
    would put the *target* in a box labelled "guide".

    ``annotation_columns`` names further columns to carry into the hover box, in
    the order given. A column absent from a row is skipped rather than shown
    empty, which is what makes a left join against a partial annotation file
    usable: unmatched guides simply carry fewer lines.
    """
    guides: list[Guide] = []
    for row in rows:
        if transcript_id is not None:
            if "transcript_id" not in row:
                raise KeyError(
                    "transcript_id filter requested but rows have no "
                    f"'transcript_id' column; got {sorted(row)}"
                )
            if str(row["transcript_id"]) != transcript_id:
                continue
        guides.append(
            Guide(
                tile_id=str(row["tile_id"]),
                start=int(row["start"]),
                end=int(row["end"]),
                score=_optional_float(row.get(score_column)) if score_column else None,
                label=_optional_str(row.get(label_column)) if label_column else None,
                sequence=(
                    _optional_str(row.get(sequence_column)) if sequence_column else None
                ),
                annotations=tuple(
                    (column, value)
                    for column in annotation_columns
                    if (value := _optional_str(row.get(column))) is not None
                ),
            )
        )
    guides.sort(key=lambda guide: (guide.start, guide.end, guide.tile_id))
    return guides


def _optional_float(value: object) -> float | None:
    if value is None or value == "":
        return None
    return float(value)  # type: ignore[arg-type]


def _optional_str(value: object) -> str | None:
    if value is None or value == "":
        return None
    return str(value)


def _sequence_line(sequence: str | None) -> str:
    """A guide sequence as a hover line, or an empty string when there is none.

    Carries its own leading line break, so an absent sequence adds nothing to the
    hover box. The ends are marked 5' and 3' because every sequence in this
    pipeline is written 5' -> 3' (see ``src/tile.py``) and a bare 30-mer in a
    tooltip gives the reader no way to confirm that. Monospaced so the bases line
    up between successive hovers.
    """
    if not sequence:
        return ""
    return (
        "<br><span style=\"font-family:monospace\">"
        f"5'-{sequence}-3'</span>"
    )


def _annotation_lines(annotations: Sequence[tuple[str, str]]) -> str:
    """Annotations as hover lines, or an empty string when there are none.

    Carries its own leading line break, so a guide with no annotations adds
    nothing to the hover box.

    Long values are broken at ``;`` but **never** at ``,``. In this pipeline's
    locus strings those two separators mean different things -- ``;`` divides
    references, ``,`` divides several loci on the *same* reference, which is the
    whole distinction between cross-family reactivity and a self-repeat (see the
    stage 2 notes). Breaking at ``;`` therefore puts one reference per line and
    keeps a reference's own loci together, rather than scattering a distinction
    the reader is meant to see.

    Only values longer than :data:`_WRAP_OVER` are broken, so a short gene list
    like ``PIK3CB;PIK3CD`` stays on one line while a multi-reference locus string
    is split. Hover boxes have no width to spare, and a needless break reads as
    two findings where there is one.
    """
    if not annotations:
        return ""
    lines = []
    for name, value in annotations:
        if len(value) > _WRAP_OVER:
            value = value.replace(";", ";<br>&nbsp;&nbsp;&nbsp;&nbsp;")
        lines.append(f"<br><i>{name}</i>: {value}")
    return "".join(lines)


def check_within_transcript(guides: Iterable[Guide], length: int) -> None:
    """Raise unless every guide lies inside ``1..length``.

    A guide off the end of its transcript means the length and the guides came
    from different sequences -- the wrong FASTA, or a stale predictions file.
    That would render as a plot that looks plausible but is measured against the
    wrong ruler, so it fails here instead.
    """
    if length < 1:
        raise ValueError(f"transcript length must be >= 1, got {length}")
    for guide in guides:
        if guide.start < 1 or guide.end > length or guide.start > guide.end:
            raise ValueError(
                f"guide {guide.tile_id} spans {guide.start}-{guide.end}, which is not "
                f"within 1-{length}; are the guides and the transcript length from "
                "the same sequence?"
            )


def pack_lanes(guides: Sequence[Guide], *, gap: int = 1) -> list[int]:
    """Assign each guide a lane so no two guides in a lane overlap.

    Returns one lane index per guide, positionally matching ``guides``. Greedy:
    guides are considered in ascending ``start`` and each goes in the lowest lane
    whose last guide ended at least ``gap`` bases earlier, so a guide only moves
    up a lane when it genuinely has nowhere below to sit. The result is not the
    provably minimum number of lanes, but for interval packing greedy-by-start
    *is* optimal, and the lane count at any position equals the depth of
    overlapping guides there.

    ``gap`` is the clear space required between two guides sharing a lane, in
    bases. The default of 1 separates guides that merely abut (``end=30`` and
    ``start=31``), which would otherwise draw as one unbroken bar.
    """
    if gap < 0:
        raise ValueError(f"gap must be >= 0, got {gap}")

    lane_ends: list[int] = []
    lanes = [0] * len(guides)
    for index in sorted(range(len(guides)), key=lambda i: (guides[i].start, guides[i].end)):
        guide = guides[index]
        for lane, last_end in enumerate(lane_ends):
            if last_end + gap < guide.start:
                lanes[index] = lane
                lane_ends[lane] = guide.end
                break
        else:
            lanes[index] = len(lane_ends)
            lane_ends.append(guide.end)
    return lanes


def guide_track(
    length: int,
    guides: Sequence[Guide],
    *,
    transcript_id: str = "transcript",
    y_mode: str = "lane",
    score_label: str = "predicted efficacy",
    colorscale: str = "Viridis",
    gap: int = 1,
    title: str | None = None,
) -> "go.Figure":
    """An interactive map of ``guides`` along a transcript of ``length`` bases.

    ``y_mode`` picks what the vertical axis means:

    ``"lane"``
        Each guide is a bar at its true width, stacked into non-overlapping rows
        (:func:`pack_lanes`). Shows position, extent and local density. This is
        the default because extent is the thing a 30-mer plot is usually asked
        about.
    ``"score"``
        Each guide is a bar at its true width placed at the height of its
        ``score``. Shows how good the guides are *and* where, at the cost of
        overlapping bars. Requires every guide to carry a score.

    Guides are coloured by ``score`` when they have one, with a colourbar;
    otherwise they are drawn in one colour. The x axis is always the full
    transcript, so empty regions read as genuinely uncovered. Hovering a guide
    gives its ``tile_id``, span, length, score and label.

    The figure is returned, not written -- see :func:`write_html`.
    """
    if y_mode not in ("lane", "score"):
        raise ValueError(f"y_mode must be 'lane' or 'score', got {y_mode!r}")
    check_within_transcript(guides, length)

    scores = [guide.score for guide in guides]
    has_scores = bool(guides) and all(score is not None for score in scores)
    if y_mode == "score" and not has_scores:
        raise ValueError(
            "y_mode='score' needs a score on every guide; "
            f"{sum(score is None for score in scores)} of {len(guides)} have none. "
            "Pass score_column= to guides_from_rows, or use y_mode='lane'."
        )

    figure = go.Figure()
    _add_transcript_band(figure, length, transcript_id)

    if guides:
        lanes = pack_lanes(guides, gap=gap)
        y = scores if y_mode == "score" else lanes
        figure.add_trace(
            go.Bar(
                x=[guide.length for guide in guides],
                base=[guide.start - 0.5 for guide in guides],
                y=y,
                orientation="h",
                width=0.72,
                marker=dict(
                    color=scores if has_scores else _DEFAULT_COLOR,
                    colorscale=colorscale if has_scores else None,
                    colorbar=dict(title=score_label, thickness=12) if has_scores else None,
                    line=dict(width=0),
                ),
                customdata=[
                    (
                        guide.tile_id,
                        guide.start,
                        guide.end,
                        guide.length,
                        "n/a" if guide.score is None else f"{guide.score:.4g}",
                        # These two carry their own line break so the template
                        # needs no conditional: when absent they add nothing to
                        # the box rather than leaving a blank line.
                        f"<br>{guide.label}" if guide.label else "",
                        _sequence_line(guide.sequence),
                        _annotation_lines(guide.annotations),
                    )
                    for guide in guides
                ],
                hovertemplate=(
                    "<b>%{customdata[0]}</b><br>"
                    "%{customdata[1]}-%{customdata[2]} (%{customdata[3]} nt)<br>"
                    f"{score_label}: %{{customdata[4]}}"
                    "%{customdata[5]}%{customdata[6]}%{customdata[7]}<extra></extra>"
                ),
                showlegend=False,
            )
        )
    else:
        figure.add_annotation(
            text="no guides passed the filter for this transcript",
            xref="paper",
            yref="paper",
            x=0.5,
            y=0.5,
            showarrow=False,
            font=dict(size=13, color="#8a949c"),
        )

    if title is None:
        title = f"{len(guides)} guide(s) on {transcript_id} ({length:,} nt)"
    figure.update_layout(
        title=title,
        xaxis=dict(
            title=f"position on {transcript_id} (nt, 1-based inclusive)",
            range=[0.5, length + 0.5],
            rangeslider=dict(visible=True, thickness=0.06),
            showgrid=True,
        ),
        yaxis=(
            dict(title=score_label)
            if y_mode == "score"
            else dict(title="guides (stacked to avoid overlap)", showticklabels=False, zeroline=False)
        ),
        bargap=0.0,
        template="plotly_white",
        hovermode="closest",
        margin=dict(l=70, r=20, t=60, b=60),
    )
    return figure


def _add_transcript_band(figure: "go.Figure", length: int, transcript_id: str) -> None:
    """Draw the transcript itself as a band pinned to the bottom of the plot."""
    low, high = _TRANSCRIPT_BAND
    figure.add_shape(
        type="rect",
        xref="x",
        yref="paper",
        x0=0.5,
        x1=length + 0.5,
        y0=low,
        y1=high,
        fillcolor=_TRANSCRIPT_COLOR,
        line=dict(width=0),
        layer="below",
    )
    figure.add_annotation(
        text=f"{transcript_id} · {length:,} nt",
        xref="paper",
        yref="paper",
        x=0.0,
        y=high,
        xanchor="left",
        yanchor="bottom",
        showarrow=False,
        font=dict(size=10, color="#6b7780"),
    )


def write_html(figure: "go.Figure", path: str, *, include_plotlyjs: str = "cdn") -> str:
    """Write ``figure`` to a standalone HTML file and return the path.

    ``include_plotlyjs="cdn"`` keeps the file small (~50 kB) but needs a network
    connection to open; pass ``True`` to inline plotly.js instead (~3 MB, opens
    offline), which is the right choice for a file being shared or archived.
    """
    figure.write_html(path, include_plotlyjs=include_plotlyjs)
    return path
