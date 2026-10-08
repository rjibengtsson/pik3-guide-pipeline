#!/usr/bin/env python3
"""Map stage-1 tiles to the gene family and label each by which genes it targets.

    python scripts/classify_guides.py --tiles results/pik3/tiles.fasta \
        --input-dir data/pik3-iso --output results/pik3/classification.tsv

Stages:

1. Build a reference from every FASTA in ``--input-dir`` (the same transcripts
   that were tiled) and index it with ``bwa index``.
2. Map the tiles with ``bwa aln``/``samse``, allowing up to ``--max-mismatches``.
3. Collapse primary and ``XA`` hits into one gene set per tile, at ``NM=0`` and
   at ``NM<=max-mismatches``, keeping every hit *locus* as
   ``accession:start-end(nm)`` so alternative and repeated positions survive.
4. Label each tile ``specific`` / ``dual`` / ``pan`` / ``other`` at both cutoffs.
5. Validate the ``NM=0`` sets *and* their positions against an exhaustive
   substring oracle.

Input is ``results/{dataset}/tiles.fasta`` -- the **target** sequences, verbatim DNA. Not
``guides.fasta``: those are reverse-complemented and in the RNA alphabet, so they
would map to the opposite strand (or not at all), and a ``U`` here is rejected.

Mapping uses ``bwa aln`` rather than ``bwa mem``, which is tuned for >=70 bp and
unreliable on 30-mers. Both cutoffs come from a **single** alignment run, filtered
two ways, so the two sets of labels cannot diverge through aligner nondeterminism.

The ``NM=0`` columns are validated exhaustively. The mismatch-tolerant columns are
**not** -- they rest on bwa's heuristic and are advisory. See ``src/classify.py``.
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from classify import (  # noqa: E402
    DUAL_PAIR,
    check_self_hits,
    classify,
    exact_hit_loci,
    exact_hit_sets,
    format_loci,
    gene_symbols,
    hit_loci,
    hit_sets,
    hits_from_sam,
    source_transcript,
)
from tile import find_fasta_files, parse_fasta, read_fasta  # noqa: E402

COLUMNS = [
    "tile_id",
    "transcript_id",
    "gene",
    "n_hits_exact",
    "hit_genes_exact",
    "hit_locs_exact",
    "class_exact",
    "n_hits_nm",
    "hit_genes_nm",
    "hit_locs_nm",
    "max_nm",
    "class_nm",
]

CLASSES = ("specific", "dual", "pan", "other", "unmapped")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--tiles",
        default="results/pik3/tiles.fasta",
        help="FASTA of target tiles from stage 1, headers are tile_ids "
        "(default: %(default)s). Must be targets, not guides.",
    )
    parser.add_argument(
        "--input-dir",
        default="data/pik3-iso",
        help="folder of the transcript FASTAs to map against; the same ones that "
        "were tiled (default: %(default)s)",
    )
    parser.add_argument(
        "--output",
        default="results/pik3/classification.tsv",
        help="TSV of one row per tile, keyed on tile_id (default: %(default)s)",
    )
    parser.add_argument(
        "--max-mismatches",
        type=int,
        default=3,
        help="edit distance for the permissive columns. These are NOT validated "
        "and are advisory only (default: %(default)s)",
    )
    parser.add_argument(
        "--work-dir",
        help="where to write the reference, bwa index and parsed hits (default: a "
        "temporary directory, removed afterwards unless --keep-temp). Index files are "
        "binary; keep them out of results/. A directory given here is always kept.",
    )
    parser.add_argument(
        "--keep-temp",
        action="store_true",
        help="keep the temporary work directory instead of deleting it, and print "
        "its path. For inspecting the reference, the bwa index, the SAM and the "
        "parsed hits; see --work-dir to choose the location yourself.",
    )
    parser.add_argument(
        "--threads", type=int, default=1, help="threads for bwa aln (default: %(default)s)"
    )
    parser.add_argument(
        "--no-validate",
        action="store_true",
        help="skip the exhaustive NM=0 substring check. Only for triage -- the "
        "check is what makes the exact columns trustworthy.",
    )
    return parser.parse_args(argv)


def require_bwa() -> str:
    path = shutil.which("bwa")
    if path is None:
        sys.exit(
            "error: bwa not found on PATH. Create and activate the project env:\n"
            "    mamba env create -f environment.yml\n"
            "    conda activate pik3-guide"
        )
    return path


def run(command: list[str], stdout_path: str | None = None) -> None:
    """Run a command, failing with its stderr rather than a bare exit code."""
    with open(stdout_path, "wb") if stdout_path else _NullContext() as stdout:
        result = subprocess.run(
            command,
            stdout=stdout if stdout_path else subprocess.DEVNULL,
            stderr=subprocess.PIPE,
        )
    if result.returncode != 0:
        sys.exit(
            f"error: {' '.join(command)} failed (exit {result.returncode}):\n"
            f"{result.stderr.decode('utf-8', 'replace').strip()}"
        )


class _NullContext:
    def __enter__(self):
        return None

    def __exit__(self, *exc):
        return False


def load_references(input_dir: str) -> tuple[dict[str, str], dict[str, str]]:
    """Return ``(sequences, genes)`` keyed by transcript id, read from ``input_dir``."""
    try:
        paths = find_fasta_files(input_dir)
    except NotADirectoryError as exc:
        sys.exit(f"error: {exc}")
    if not paths:
        sys.exit(f"error: no FASTA files in {input_dir}")

    records = [record for path in paths for record in read_fasta(path)]
    if not records:
        sys.exit(f"error: no FASTA records in {input_dir}")

    sequences: dict[str, str] = {}
    for record in records:
        if record.transcript_id in sequences:
            sys.exit(
                f"error: transcript id {record.transcript_id} appears twice in "
                f"{input_dir}; hit sets would be ambiguous"
            )
        sequences[record.transcript_id] = record.sequence
    try:
        genes = gene_symbols(records)
    except ValueError as exc:
        sys.exit(f"error: {exc}")

    duplicates = {gene for gene in genes.values() if list(genes.values()).count(gene) > 1}
    if duplicates:
        sys.exit(
            f"error: gene symbol(s) {sorted(duplicates)} map to more than one transcript. "
            "Per-gene labels assume one transcript per gene; see the scope note in "
            "src/classify.py."
        )
    return sequences, genes


def load_tiles(path: str) -> dict[str, str]:
    """Read the tile FASTA, rejecting guides (RNA alphabet) and duplicate ids."""
    if not os.path.isfile(path):
        sys.exit(f"error: no such file: {path}")
    tiles: dict[str, str] = {}
    with open(path, "r", encoding="utf-8") as handle:
        for record in parse_fasta(handle, source_file=os.path.basename(path)):
            if "U" in record.sequence:
                sys.exit(
                    f"error: {path}: {record.transcript_id} contains U, so this looks "
                    "like guides.fasta. Map the targets (tiles.fasta) instead -- guides "
                    "are reverse-complemented and would map to the wrong strand."
                )
            if record.transcript_id in tiles:
                sys.exit(f"error: {path}: duplicate tile id {record.transcript_id}")
            tiles[record.transcript_id] = record.sequence
    if not tiles:
        sys.exit(f"error: {path}: no records")
    return tiles


def write_reference(sequences: dict[str, str], path: str) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        for transcript_id, sequence in sequences.items():
            handle.write(f">{transcript_id}\n{sequence}\n")


HIT_COLUMNS = [
    "tile_id",
    "reference",
    "position",
    "nm",
    "strand",
    "in_bounds",
    "counted_exact",
    "counted_nm",
]


def dump_hits(hits: list, path: str, max_nm: int) -> None:
    """Write every parsed hit, before deduplication, with why it was kept or not.

    The intermediate worth having when a hit set looks wrong: one row per primary
    or ``XA`` alignment, so ``XA`` parsing, reverse-strand drops and the mismatch
    cutoff are all visible separately. ``in_bounds`` is 0 for an alignment running
    off the end of its reference -- a concatenated-index artifact, not a match. ``hit_sets`` collapses these to one row per
    tile, which is where the detail is otherwise lost -- in particular a tile that
    hits its own transcript at several positions appears once there and several
    times here. The ``hit_locs_*`` columns of the output carry the surviving
    positions; this dump also holds the ones that were filtered out.
    """
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=HIT_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for hit in hits:
            usable = not hit.is_reverse and not hit.overhangs
            writer.writerow(
                {
                    "tile_id": hit.query_id,
                    "reference": hit.reference,
                    "position": hit.position,
                    "nm": hit.nm,
                    "strand": "-" if hit.is_reverse else "+",
                    "in_bounds": int(not hit.overhangs),
                    "counted_exact": int(usable and hit.nm == 0),
                    "counted_nm": int(usable and hit.nm <= max_nm),
                }
            )


def map_tiles(tiles_path: str, reference: str, work_dir: str, args: argparse.Namespace) -> str:
    """Index the reference, map the tiles, and return the path to the SAM."""
    sai = os.path.join(work_dir, "tiles.sai")
    sam = os.path.join(work_dir, "tiles.sam")
    mismatches = str(args.max_mismatches)

    run(["bwa", "index", reference])
    run(
        [
            "bwa", "aln",
            # Report *every* hit at this edit distance. Without -N bwa stops early
            # and the hit set is incomplete, which would mislabel pan as specific.
            "-N",
            "-n", mismatches,
            # With the seed spanning the whole read (-l below), -k caps total
            # mismatches, so its default of 2 would silently override -n.
            "-k", mismatches,
            # Seed >= read length: seeding cannot limit sensitivity.
            "-l", "30",
            # No gaps; a 30-mer either matches or does not.
            "-o", "0",
            "-t", str(args.threads),
            reference,
            tiles_path,
        ],
        stdout_path=sai,
    )
    run(
        # bwa truncates XA at 3 alternative hits by default, which would drop real
        # hits for a conserved tile. A generous cap costs nothing at this scale.
        ["bwa", "samse", "-n", "1000", reference, sai, tiles_path],
        stdout_path=sam,
    )
    return sam


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.max_mismatches < 0:
        sys.exit(f"error: --max-mismatches must be >= 0, got {args.max_mismatches}")
    if args.threads < 1:
        sys.exit(f"error: --threads must be >= 1, got {args.threads}")
    require_bwa()

    sequences, genes = load_references(args.input_dir)
    tiles = load_tiles(args.tiles)

    unknown = {source_transcript(tile_id) for tile_id in tiles} - set(sequences)
    if unknown:
        sys.exit(
            f"error: tiles reference transcript(s) absent from {args.input_dir}: "
            f"{sorted(unknown)}. Tiles and reference must come from the same input."
        )

    # An explicit --work-dir is always kept; a temporary one is kept only with
    # --keep-temp, which uses mkdtemp so nothing deletes it on the way out.
    work_dir = args.work_dir
    temporary = None
    if work_dir is not None:
        os.makedirs(work_dir, exist_ok=True)
    elif args.keep_temp:
        work_dir = tempfile.mkdtemp(prefix="classify-guides-")
    else:
        temporary = tempfile.TemporaryDirectory(prefix="classify-guides-")
        work_dir = temporary.name
    kept = temporary is None

    try:
        reference = os.path.join(work_dir, "transcriptome.fasta")
        write_reference(sequences, reference)
        sam = map_tiles(args.tiles, reference, work_dir, args)
        hits = list(hits_from_sam(sam))
        dump_hits(hits, os.path.join(work_dir, "hits.tsv"), args.max_mismatches)
    finally:
        if temporary is not None:
            temporary.cleanup()

    all_genes = frozenset(genes.values())
    exact, _, reverse_dropped = hit_sets(hits, max_nm=0)
    permissive, worst_nm, _ = hit_sets(hits, max_nm=args.max_mismatches)
    exact_loci = hit_loci(hits, max_nm=0)
    permissive_loci = hit_loci(hits, max_nm=args.max_mismatches)

    missing = check_self_hits({tile_id: exact.get(tile_id, frozenset()) for tile_id in tiles})
    if missing:
        sys.exit(
            f"error: {len(missing)} tile(s) did not match their own source transcript "
            f"at NM=0, e.g. {missing[:3]}. Every tile was cut from the reference, so "
            "bwa lost a true alignment and these hit sets are incomplete."
        )

    if not args.no_validate:
        oracle = exact_hit_sets(tiles, sequences)
        disagreements = sorted(
            tile_id
            for tile_id in tiles
            if exact.get(tile_id, frozenset()) != oracle[tile_id]
        )
        if disagreements:
            example = disagreements[0]
            sys.exit(
                f"error: bwa's NM=0 hits disagree with the substring oracle for "
                f"{len(disagreements)} tile(s). First: {example} -- bwa "
                f"{sorted(exact.get(example, frozenset()))} vs oracle "
                f"{sorted(oracle[example])}."
            )

        # The locus columns get the same treatment as the sets they sit beside:
        # compared against positions found by substring search, which is complete
        # at NM=0. Positions only, since the oracle says nothing about NM>0 loci.
        oracle_loci = exact_hit_loci(tiles, sequences)
        misplaced = sorted(
            tile_id
            for tile_id in tiles
            if {
                reference: [position for position, _nm in loci]
                for reference, loci in exact_loci.get(tile_id, {}).items()
            }
            != oracle_loci[tile_id]
        )
        if misplaced:
            example = misplaced[0]
            sys.exit(
                f"error: bwa's NM=0 hit positions disagree with the substring oracle "
                f"for {len(misplaced)} tile(s). First: {example} -- bwa "
                f"{exact_loci.get(example, {})} vs oracle {oracle_loci[example]}."
            )

    rows = []
    for tile_id, tile in tiles.items():
        transcript_id = source_transcript(tile_id)
        # Alignments are ungapped, so each locus spans exactly the tile's length.
        length = len(tile)
        exact_genes = frozenset(genes[ref] for ref in exact.get(tile_id, frozenset()))
        permissive_genes = frozenset(genes[ref] for ref in permissive.get(tile_id, frozenset()))
        rows.append(
            {
                "tile_id": tile_id,
                "transcript_id": transcript_id,
                "gene": genes[transcript_id],
                "n_hits_exact": len(exact_genes),
                "hit_genes_exact": ";".join(sorted(exact_genes)),
                "hit_locs_exact": format_loci(exact_loci.get(tile_id, {}), length),
                "class_exact": classify(exact_genes, all_genes),
                "n_hits_nm": len(permissive_genes),
                "hit_genes_nm": ";".join(sorted(permissive_genes)),
                "hit_locs_nm": format_loci(permissive_loci.get(tile_id, {}), length),
                "max_nm": worst_nm.get(tile_id, 0),
                "class_nm": classify(permissive_genes, all_genes),
            }
        )

    directory = os.path.dirname(os.path.abspath(args.output))
    os.makedirs(directory, exist_ok=True)
    with open(args.output, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    overhanging = sum(1 for hit in hits if hit.overhangs and not hit.is_reverse)
    report(
        rows, args, all_genes, reverse_dropped, work_dir if kept else None, overhanging
    )
    return 0


def report(
    rows: list[dict[str, object]],
    args: argparse.Namespace,
    all_genes: frozenset[str],
    reverse_dropped: int,
    work_dir: str | None = None,
    overhanging: int = 0,
) -> None:
    validated = "skipped (--no-validate)" if args.no_validate else "passed"
    print(
        f"reference: {len(all_genes)} genes ({', '.join(sorted(all_genes))})  "
        f"tiles: {len(rows)}",
        file=sys.stderr,
    )
    print(
        f"dual pair: {'/'.join(sorted(DUAL_PAIR))}  "
        f"NM=0 oracle check: {validated}  "
        f"permissive cutoff: NM<={args.max_mismatches} (NOT validated, advisory)",
        file=sys.stderr,
    )
    if overhanging:
        print(
            f"note: dropped {overhanging} hit(s) running off the end of their "
            "reference; bwa indexes the transcripts concatenated, so an alignment "
            "can cross into the next one and still be attributed to the previous. "
            "See in_bounds in hits.tsv.",
            file=sys.stderr,
        )
    if reverse_dropped:
        print(
            f"note: dropped {reverse_dropped} reverse-strand hit(s); the reference is "
            "mRNA sense sequence. A large count means guides were mapped, not targets.",
            file=sys.stderr,
        )

    by_gene: dict[str, dict[str, int]] = {}
    for row in rows:
        counts = by_gene.setdefault(str(row["gene"]), {name: 0 for name in CLASSES})
        counts[str(row["class_exact"])] += 1

    header = f"{'gene':<10}" + "".join(f"{name:>10}" for name in CLASSES)
    print("NM=0 (validated):", file=sys.stderr)
    print(header, file=sys.stderr)
    for gene in sorted(by_gene):
        counts = by_gene[gene]
        print(
            f"{gene:<10}" + "".join(f"{counts[name]:>10}" for name in CLASSES),
            file=sys.stderr,
        )
    totals = {name: sum(c[name] for c in by_gene.values()) for name in CLASSES}
    print(f"{'TOTAL':<10}" + "".join(f"{totals[name]:>10}" for name in CLASSES), file=sys.stderr)

    # The row that matters: looks isoform-specific exactly, but has near-matches
    # elsewhere in the family. Surfaced here so it is not left to be discovered by
    # filtering the TSV later.
    degraded = [
        row for row in rows if row["class_exact"] == "specific" and row["class_nm"] != "specific"
    ]
    print(
        f"specific at NM=0 but not at NM<={args.max_mismatches}: {len(degraded)} tile(s)"
        " - near-matches elsewhere in the family, treat as unconfirmed",
        file=sys.stderr,
    )
    print(f"wrote {args.output}", file=sys.stderr)
    if work_dir is not None:
        print(f"kept intermediates in {work_dir}", file=sys.stderr)
        for name, what in (
            ("transcriptome.fasta", "the reference that was indexed"),
            ("transcriptome.fasta.{amb,ann,bwt,pac,sa}", "bwa index"),
            ("tiles.sai", "bwa aln output"),
            ("tiles.sam", "bwa samse output, with the XA tags"),
            ("hits.tsv", "every parsed hit before dedup, with why it counted"),
        ):
            print(f"    {name:<38} {what}", file=sys.stderr)
    else:
        print(
            "intermediates discarded; pass --keep-temp to inspect them "
            "(or --work-dir to choose where)",
            file=sys.stderr,
        )


if __name__ == "__main__":
    raise SystemExit(main())
