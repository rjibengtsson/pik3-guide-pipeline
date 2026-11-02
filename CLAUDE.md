# pik3-guide-pipeline

Pipeline for designing CRISPR-Cas13 guide RNAs that target a single gene isoform, a
defined subset, or all members of a gene family. The worked example is the class IA/IB
PI3K catalytic subunits (PIK3CA/B/D/G).

Stages 1 (tiling + guide generation) and 2 (family classification) exist, plus
`scripts/plot_scatterlines.py` for visualization. Remaining stages — GC/MFE
filtering, genome-wide off-target mapping — are not written yet.

Inputs and outputs are organised **per dataset**: `data/{dataset}/` is tiled into
`results/{dataset}/`. `pik3-iso`/`pik3` is the worked example (the four PI3K
paralogs) and is what every CLI default and the real-data tests point at; `ppib` is
a second, single-gene set. When adding a dataset, pass `--input-dir` and `--output`
explicitly — the defaults name pik3, so omitting them overwrites the PI3K run.

## Commands

Everything must run inside the project conda env; see **Environment** below for why.

```bash
mamba env create -f environment.yml   # or: conda env create -f environment.yml
conda activate pik3-guide

# tests (320, all passing; needs bwa and blastn, else some skip)
python -m unittest discover -s test -v

# stage 1, exactly as the current results/pik3/ was produced (git-ignored; rerun to recreate)
python scripts/tile_transcripts.py --input-dir data/pik3-iso \
    --output results/pik3/tiles.tsv --input-type mrna -k 30 --step 1 \
    --per-transcript --fasta-out results/pik3/guides.fasta


## Layout

- [src/tile.py](src/tile.py) — the library: FASTA parsing, tiling, the guide
  convention. Holds the module-level docstring that defines the coordinate,
  direction and alphabet conventions; read it before changing anything here.
- [scripts/tile_transcripts.py](scripts/tile_transcripts.py) — the CLI around it:
  argument parsing, row assembly, TSV/FASTA writing, the stderr run summary.
- [test/test_tile.py](test/test_tile.py) — unittest suite covering both.
- [data/pik3-iso/](data/pik3-iso/) — input RefSeq mRNA FASTA, one transcript per
  file, one folder per dataset.
- [results/pik3/](results/pik3/) — generated output for that dataset, git-ignored (see **Outputs**).

`src/` and `scripts/` are not a package. Both the CLI and the tests reach `tile` by
`sys.path.insert` of the repo's `src/` (and `scripts/`) directory — see
[scripts/tile_transcripts.py:29](scripts/tile_transcripts.py#L29) and
[test/test_tile.py:15-17](test/test_tile.py#L15-L17). There is no `setup.py`,
`pyproject.toml`, linter config or CI.

## Conventions that must not drift

These are load-bearing across stages, asserted by tests, and documented in three
places (module docstring, README, CLI docstring). Changing one means changing all of
them.

- **Coordinates are 1-based inclusive**, so `sequence[start - 1:end] == tile.sequence`
  exactly. Tested by `test_slicing_round_trip` and
  `test_output_is_deterministic_and_coordinates_round_trip`.
- **`tile_id` is `{transcript_id}:{start}-{end}`** (e.g. `NM_006218.4:1-30`) and is the
  join key every downstream stage is expected to use. Duplicate transcript ids across
  input files are a hard error precisely because tile_ids would collide.
- **Every sequence in and out is written 5' -> 3'.** Guide and target are antiparallel,
  so the guide's 5' end pairs with the target's 3' end.
- **`target_seq` is verbatim from the input**, in whatever alphabet the FASTA used — a
  DNA-alphabet RefSeq mRNA keeps its `T`. Nothing is transcribed on the target side.
- **`guide_rna` is transcribed**, because a guide base-pairs as RNA: `A` in the target
  gives `U` in the guide. `--guide-alphabet dna` writes `T` instead, for aligners
  (bowtie/blast) that expect DNA; `target_seq` is unaffected either way.
- **`strand` is always `+`** — mRNA is single-stranded.
- Sequence manipulation is **Biopython's, not ours**. `transcribe`,
  `back_transcribe`, `reverse_complement` and `reverse_complement_rna` are re-exported
  from `src/tile.py` so callers need one import. `guide_rna()` is a thin wrapper over
  `reverse_complement_rna` that exists to name the biology and carry the direction
  convention. `TestSequenceHelpers.test_biopython_helpers_behave_as_this_module_assumes`
  pins the Biopython semantics so an upgrade fails loudly.
- `pairs_antiparallel` in the test suite is **deliberately not written in terms of
  `reverse_complement`**, so it fails rather than agrees if that function is wrong.
  Keep it independent. It compares in RNA space, so a DNA target passes.

## Environment

`Bio.SeqIO` needs Biopython **and** numpy: `Bio/SeqIO/__init__.py` imports `TwoBitIO`
unconditionally, so a missing numpy breaks plain FASTA parsing too. The base conda env
has no numpy. `src/tile.py` catches that `ImportError` and re-raises it with the
`environment.yml` instructions, because Biopython's own message talks about TwoBit
files and is misleading. Don't "simplify" that handler away.

## Outputs

`results/pik3/` holds the run over the four PI3K transcripts: 28,032 tiles at k=30,
step=1. It is **git-ignored** (see **Known rough edges**), so these files exist only
where a stage has been run.

- `tiles.tsv` — combined, columns `tile_id, transcript_id, source_file, start, end,
  strand, length, target_seq, guide_rna, has_ambiguous`.
- `tiles/{transcript_id}.tsv` — same columns, one file per transcript (`--per-transcript`).
- `guides.fasta` — guide sequences, headers set to `tile_id`.
- `tiles.fasta` — the same tiles written with `--fasta-field target_seq`, i.e. targets
  rather than guides. The two FASTA files are **not** interchangeable.
- `guide_qc.tsv` — per-guide homopolymer and GC flags over `guides.fasta` (see
  **Guide-level QC**). Joins onto `tiles.tsv` on `guide_id` = `tile_id`.

Output order is deterministic: input files sorted by path, records in file order, tiles
by ascending start.

## Known rough edges

- **`results/` is ignored and untracked, so no output in it is in the repo.**
  `.gitignore` line 7 is `results/`, which does match (an earlier `./results` did not),
  and `git ls-files results/` is empty. Every file named under **Outputs** has to be
  regenerated by rerunning its stage command — they are reproducible, which is the
  stated reason for ignoring them, but do not expect to find them in a fresh clone, and
  do not describe them as committed. `data/` is still tracked.
- **`--input-type` is required but never verified**, and currently cannot change the
  output: `target_seq` is verbatim and `guide_rna` goes through
  `reverse_complement_rna`, which accepts both alphabets, so `dna` and `mrna` produce
  identical rows (asserted by `test_declared_type_does_not_change_the_output`). The
  flag documents the input and is echoed in the run summary.
- **Genomic DNA is not handled.** "DNA" means the coding/sense strand of a transcript;
  splice out introns and orient the sequence before tiling.
- Sequences mixing `T` and `U`, or holding non-nucleotide characters, pass straight
  through. Only `has_ambiguous` (any non-ACGTU base in the tile) is reported.
- `find_fasta_files` accepts `.txt` alongside `.fasta/.fa/.fna`, but the CLI's
  "no FASTA files" error message lists only the latter three.

## Stage 2 — family classification

Label each tile by which members of the gene family it targets.

- [src/classify.py](src/classify.py) — pure logic: hit sets, labels, the exact-match
  oracle. No aligner in any signature, so the rules are testable without bwa.
- [scripts/classify_guides.py](scripts/classify_guides.py) — CLI: builds the reference,
  runs bwa, parses, validates, writes.
- [test/test_classify.py](test/test_classify.py) — unittest; the bwa-dependent cases
  skip when bwa is absent.

```bash
python scripts/classify_guides.py --tiles results/pik3/tiles.fasta \
    --input-dir data/pik3-iso --output results/pik3/classification.tsv \
    --work-dir /tmp/bwa-work --threads 4
```

Runs in ~3 s over 28,032 tiles.

### Inspecting the intermediates

By default the work directory is temporary and deleted. To keep it:

```bash
python scripts/classify_guides.py ... --keep-temp    # mkdtemp, path printed, not deleted
python scripts/classify_guides.py ... --work-dir /tmp/bwa-work   # your own location, always kept
```

Either way the summary names the directory and lists what is in it. Never point
`--work-dir` at `results/` — the bwa index files are binary.

| file | what it is |
|---|---|
| `transcriptome.fasta` | the reference that was indexed |
| `transcriptome.fasta.{amb,ann,bwt,pac,sa}` | bwa index |
| `tiles.sai` | `bwa aln` output |
| `tiles.sam` | `bwa samse` output, with the `XA` tags |
| `hits.tsv` | every parsed hit **before** dedup, with why it counted |

`hits.tsv` is the one to reach for first. `hit_sets` collapses hits to one row per
tile, so the dump is where `XA` parsing, reverse-strand drops and the mismatch cutoff
are each visible separately — columns `in_bounds`, `counted_exact` and `counted_nm`
say which filter excluded a hit. Columns are `tile_id, reference, position, nm, strand,
in_bounds, counted_exact, counted_nm`. It also distinguishes the two cases that
`max_nm` conflates:

```
NM_002649.3:5952-5981   NM_002649.3  5952  0  +  1  1  1
NM_002649.3:5952-5981   NM_002649.3  6004  3  +  1  0  1   <- same reference: self-repeat
NM_006218.4:2724-2753   NM_006218.4  2724  0  +  1  1  1
NM_006218.4:2724-2753   NM_006219.3  2765  3  +  1  0  1   <- different reference: real cross-family
```

Since the output now carries hit positions, `hits.tsv` is also where a *dropped*
position shows up: `in_bounds 0` marks an alignment running off the end of its
reference (see **bwa reports alignments that run past a reference's end** below).

Input is `results/{dataset}/tiles.fasta`, **not** `guides.fasta`: tiles.fasta holds `target_seq`
as verbatim DNA, so hits land on the forward strand and a target hit is a
guide-binding site. A FASTA containing `U` should be rejected with a clear error.

Mapping is `bwa aln` (not `mem` — that is tuned for >=70 bp and unreliable on 30-mers)
against a concatenated 4-transcript reference, with `-N` for all hits, `-o 0` for no
gaps, `-l 30` so the seed spans the read, and `-k` equal to `-n` (with a full-read seed,
`-k` caps total mismatches and its default of 2 would silently override a higher `-n`).
`samse -n` must be raised well above bwa's default `XA` cap of 3. The full hit set per
tile is primary alignment plus every `XA` entry, deduplicated by reference name, with
per-hit `NM` parsed from the tag rather than inherited from the primary. Hit *positions*
are kept alongside the sets, deduplicated by `(reference, position, nm)` since bwa may
repeat the primary as an `XA` entry.

Taxonomy: `specific` (1 gene), `dual` (**only** the `PIK3CD;PIK3CG` pair), `pan` (all
4), `other` (any other 2- or 3-gene combination). The raw sorted signature stays in a
`hit_genes` column on every row, so `other` is never a dead end.

### Output columns: gene symbols and accession loci side by side

`classification.tsv` reports each hit set twice, because the two serve different
purposes and must not be collapsed into one another:

- `hit_genes_exact` / `hit_genes_nm` — sorted **gene symbols**, `;`-joined. This is
  what `class_*` and `n_hits_*` are computed from, so the label can always be audited
  from its own row. The taxonomy is gene-level by definition: `specific`/`dual`/`pan`
  only mean anything over genes.
- `hit_locs_exact` / `hit_locs_nm` — the same hits as **accessions with coordinates**,
  `accession:start-end(nm)`, `;` between references and `,` between several loci on one
  reference. Coordinates are 1-based inclusive, the same convention as `tile_id`, so a
  locus string reads directly as a coordinate on the reference and
  `reference[start - 1:end]` is the matched sequence. Alignments are ungapped (`-o 0`),
  so `end` is always `start + k - 1`.

A gene symbol cannot address a position, and with more than one transcript variant per
gene it will not even address a sequence — which is why the locus columns are
accession-keyed. Example, from the real run:

```
tile_id                 hit_genes_nm     hit_locs_nm                                         max_nm  class_nm
NM_006218.4:2724-2753   PIK3CA;PIK3CB    NM_006218.4:2724-2753(0);NM_006219.3:2765-2794(3)   3       other
NM_002649.3:5952-5981   PIK3CG           NM_002649.3:5952-5981(0),NM_002649.3:6004-6033(3)    3       specific
```

The second row is the case `max_nm` alone cannot explain: `specific`, one gene, yet
`max_nm 3` — because the tile matches its *own* transcript a second time 52 nt
downstream, in the PIK3CG tandem repeat. Comma versus semicolon is the whole
distinction between a self-repeat and cross-family reactivity.

The `n_hits_*` columns count **genes**, not loci, so `n_hits_nm` is 1 on that row
while `hit_locs_nm` lists two positions. Count loci by splitting on both separators.

Both NM=0 and NM<=3 are reported, from a **single** bwa run at `-n 3 -k 3` filtered two
ways — not two runs, so the two views cannot diverge through aligner nondeterminism.

### Validation scope: NM=0 only

The NM=0 columns -- sets **and** positions -- are validated against an
**exact-substring oracle** — with four
transcripts and 30-mers, `str.find` over each reference gives provably complete ground
truth in well under a second. bwa's NM=0 hit sets must equal `exact_hit_sets`', and its
NM=0 positions must equal `exact_hit_loci`' (every occurrence, overlaps included). This is the
same independent-verification trick as `pairs_antiparallel`: written without reusing the
thing it checks, so it fails rather than agrees when that thing is wrong.

**The NM<=3 columns are deliberately not validated.** Enumerating every <=3-mismatch
match of a 30-mer is far too large for an oracle, so those columns rest on bwa's
heuristic alone and should be read as advisory, not as a complete cross-reactivity
census. Treat an absent NM<=3 hit as "not found", never as "does not exist". Coverage
there is limited to a small hand-built reference with mismatches planted at known
positions, including at the first and last base, where aligners tend to clip rather
than report a mismatch.

Two invariants are asserted at runtime, not only in tests: every tile must hit its own
source transcript at NM=0 (it came from there — a missing self-hit means bwa dropped a
true alignment, so the run fails rather than emitting confident labels over incomplete
data), and every input tile appears exactly once in the output.

### bwa sets the unmapped flag on reads that *do* have alignments

The non-obvious trap, and the reason `hits_from_sam` gates on `RNAME` being set rather
than on `read.is_unmapped`. For a read with many equally-best hits, bwa `samse` sets
flag `0x4` while still filling in `RNAME`, `POS`, `NM:i:0`, `MD:Z:30` and a full `XA`
list of perfect matches:

```
NM_000001.1:197-226  4  NM_000001.1  213  0  30M  ...  NM:i:0  X0:i:17  XA:Z:...16 perfect hits...
```

Trusting `is_unmapped` there silently discards perfect alignments — including the tile's
own source transcript — and `XA` is never reached if the check `continue`s first. This
was caught by the self-hit invariant, on a low-complexity tile. Real mRNA has poly-A
tails and tandem repeats, so it is not a synthetic-only concern;
`TestCliLowComplexity` and `test_unmapped_flag_with_a_reported_alignment_is_still_a_hit`
pin it. Accepting these records is safe only because the NM=0 oracle independently
confirms them.

### bwa reports alignments that run past a reference's end

The second non-obvious bwa trap, found by extending the NM=0 oracle to positions.
bwa indexes the transcripts as **one concatenated sequence**, so an `XA` entry can
describe an alignment that crosses the junction into the next transcript while still
being attributed to the previous one, with a full-length CIGAR and `NM:i:0`:

```
NM_000002.1:32-61 ... XA:Z:NM_000002.1,+212,30M,0;     <- that reference is 240 nt
```

A 30M at 212 ends at 241, so the last base came from a different transcript: within the
named reference the match does not exist, and the substring oracle rightly disagrees.
`hits_from_sam` flags these (`Hit.overhangs`, via `reference_span` on the CIGAR) rather
than dropping them, so `hits.tsv` shows them with `in_bounds 0`; `hit_sets` and
`hit_loci` both skip them, and the run summary counts them.

This was **invisible to the set-level oracle**: `hit_sets` deduplicates by reference
name, so a phantom locus on a reference the tile legitimately hits elsewhere changes
nothing about the set. Only the position check sees it. It can in principle inject a
reference a tile never matched into its gene set — caught loudly by the set oracle in
that case, but silently wrong at NM>0, where there is no oracle. The real four-transcript
run produces none of these; the synthetic fixture in `TestCli` does, because its leading
and trailing filler are the same 60 nt and the second copy ends flush with the
transcript.

### Independent blastn cross-check

[test/test_blast_validation.py](test/test_blast_validation.py) reproduces the NM=0 hit
sets with blastn, which shares none of this repo's assumptions — so a hit set is wrong
only if bwa, the substring oracle and blastn are all wrong the same way. It skips
cleanly when blast is absent or when `results/pik3/` has not been generated.

It is pinned to the `pik3-iso`/`pik3` dataset by `DATA_DIR`/`RESULTS_DIR` at the top
of the file, not to whatever is in `results/`, because its expectations are specific
to the four PI3K paralogs (one gene per tile). **A skip here is silent** — if those
paths move, the strongest independent check stops running without failing. Check the
skip count.

```bash
python -m unittest discover -s test -p 'test_blast_validation.py' -v   # 400-tile sample, ~2 s
PIK3_BLAST_SAMPLE=all python -m unittest discover -s test \
    -p 'test_blast_validation.py'                                      # all 28,032, ~98 s
```

Verified to fail on a deliberately corrupted `classification.tsv` row, so it is not
passing vacuously. The sample is seeded, so a failure reproduces.

Three blastn flags are load-bearing, and all three *lose* true hits when wrong:

- `-task blastn-short` — default blastn uses word size 11 and will not reliably find a
  30 bp match at all.
- `-dust no -soft_masking false` — DUST masks low-complexity query sequence by default,
  which silently drops the PIK3CG tandem repeat around nt 5950-6100 and manufactures a
  disagreement that looks like a classifier bug.
- Filter on `length == 30 and mismatch == 0 and gapopen == 0`, not on `pident` alone: a
  9 bp fragment is also "100% identical". The full run emits ~2.8M HSPs for 28k tiles,
  of which ~28k are full-length. The length filter is what makes this equivalent to
  NM=0.

### Finding: there are no pan- or dual-guides at NM=0

Over the four transcripts in `data/`, **all 28,032 tiles are `specific`** — zero `pan`,
zero `dual`, zero `other`. The paralogs share no exact 30-mer at all: every pairwise
intersection is empty, and so is the four-way intersection, down to k=10. Verified
four independent ways — bwa, the substring oracle, a direct k-mer set intersection, and
blastn — with zero disagreements across all 28,032 tiles. blastn finds 28,046
full-length exact hits for 28,032 tiles; the 14 extra are the 7 PIK3CG sequences that
occur at two positions each.

This is expected on reflection — codon degeneracy means protein-level conservation does
not imply nucleotide identity — but it has a hard consequence: **a conserved pan-family
or PIK3CD/PIK3CG dual guide cannot be obtained by exact matching at k=30.** Achieving
one needs mismatch tolerance, a shorter spacer, or protein-level/alignment-guided region
selection, not more tiling. At NM<=3 only 14 tiles pick up any cross-family near-match
(all `other`, none the designed CD/CG pair), so relaxing the cutoff does not rescue it
either.

Practical note for stage 3: the `dual`/`pan`/`other` branches of `classify` are
therefore **unexercised by the real data** — only synthetic fixtures cover them. Don't
read the real run as evidence those paths work.

Scope limitation to keep in the module docstring: the reference is four *paralogous
genes*, one RefSeq variant each, so this measures within-family cross-reactivity only.
"PIK3CA-specific" means specific among these four sequences — not genome-wide, and not
isoform-level within a gene unless the other transcript variants are added to the
reference.


## Guide-level QC — homopolymer runs and GC band

Flags individual guides on two intrinsic sequence properties. This is **guide scale**,
not insert scale, and the distinction is the organising idea of `src/seqqc.py`: the same
file also holds the whole-insert checks behind `scripts/check_target_sites.py`, and the
two share scanners but never thresholds.

- [src/seqqc.py](src/seqqc.py) — `check_guide`, `guide_homopolymer_runs`, `GuideQC`, and
  the `DEFAULT_GUIDE_*` constants, in a section at the bottom of the file under their own
  heading. The insert-scale functions above them are unchanged.
- [scripts/check_guides.py](scripts/check_guides.py) — the CLI.
- [test/test_seqqc.py](test/test_seqqc.py) — `TestGuideHomopolymerRuns`, `TestCheckGuide`,
  `TestCheckGuidesCli`, `TestCheckGuidesCliRejectsBadInput`.

```bash
python scripts/check_guides.py --guides results/pik3/guides.fasta \
    --output results/pik3/guide_qc.tsv
```

Runs in ~2 s over 28,032 guides.

Input is `results/{dataset}/guides.fasta` — the guide file, **not** `tiles.fasta`. The
two are not interchangeable (see **Outputs**): GC and homopolymer content are being
judged of the guide that gets synthesised, and `tiles.fasta` holds targets. Headers are
`tile_id`, so the output joins back onto `tiles.tsv` on `guide_id`.

### The three flags

| flag | raised when |
|---|---|
| `gc_low` | GC fraction below `--gc-min`, default 0.40 |
| `gc_high` | GC fraction above `--gc-max`, default 0.60 |
| `homopolymer` | a run of >= `--homopolymer-min` identical bases, default 4 |

`gc_low` and `gc_high` are deliberately **two flags rather than one
`gc_out_of_band`** — they have opposite remedies, and the split is what distinguishes a
badly set band from a GC-skewed transcript. They are mutually exclusive; a guide can
still raise one of them together with `homopolymer`.

Band edges are **inclusive**, matching the insert-scale check: exactly 40.0% GC passes a
0.40 minimum. The edge is reachable — 12 G/C in a 30 nt guide is exactly 0.40 — so which
side passes is a real decision rather than a rounding detail, and
`test_exactly_forty_percent_passes` pins it. It is exact in floating point too: `12/30`
and the literal `0.4` round to the same double.

The homopolymer scan covers **all four bases**. Unlike the insert-scale check there is no
separate poly-T scan to hand T off to, because a run inside a spacer is a property of
that spacer, not a U6 termination signal in the construct — so T is not excluded here
the way it is there.

### Thresholds differ by scale, scanners do not

`DEFAULT_GUIDE_GC_MIN`/`MAX` (0.40–0.60) and `DEFAULT_GUIDE_HOMOPOLYMER_MIN` (4) are
separate constants from the insert's `DEFAULT_GC_MIN`/`MAX` (0.30–0.70) and
`DEFAULT_HOMOPOLYMER_MIN` (5), and `test_guide_scale_thresholds_are_not_the_insert_scale_ones`
asserts they stay distinct. A run of 4 is unremarkable in 2 kb of insert and
disqualifying in a 30 nt spacer.

But `check_guide` calls the *same* `gc_fraction` and `homopolymer_runs` the insert checks
use, so the two scales cannot disagree about what a run or a GC fraction **is**, only
about which values are acceptable.

`test_no_sliding_window_api_remains` still stands and is not in tension with this: it
forbids scanning an insert at spacer width, which is a guess at where a guide might go.
`check_guide` takes a guide someone has committed to. Only the second is worth a verdict.

### A run of four U must be reported as U, not T

`homopolymer_runs` normalises U to T before scanning — correct for *finding* the run,
wrong for *naming* it. Reported verbatim it would put `T:1-4(x4)` on an RNA guide, and
`TTTT` appears nowhere in that guide, so grepping the column for `UUUU` would miss its
own hit. `guide_homopolymer_runs` exists for exactly this: it reads the base back out of
the input, so a guide FASTA reports `U` and a `--guide-alphabet dna` one reports `T`.
Both alphabets give the same flags and the same spans; only the name differs
(`test_both_alphabets_give_the_same_verdict`).

A run identical only *after* normalisation (`UUTT`, from a sequence mixing alphabets) is
reported under the normalised base `T`, since no single input character describes it. Such
sequences are a known rough edge of this pipeline rather than a supported input; the span
is still right, and the span is the actionable part.

### Output columns

`guide_qc.tsv` — one row per guide, every guide by default. Columns are `guide_id,
length, sequence, gc_frac, gc_percent, gc_band, n_ambiguous, n_homopolymer_runs,
longest_homopolymer, longest_homopolymer_base, longest_homopolymer_at, homopolymer_min,
homopolymer_runs, flags, n_flags`.

```
guide_id             sequence                        gc_frac  n_homopolymer_runs  homopolymer_runs       flags               n_flags
NM_002649.3:1-30     UUUCAGAUAAUCUAGCCGAGAAGGAAGUGC  0.4333   0                                                              0
NM_002649.3:121-150  AAAAAUGUGAUGCACCUGCAAAACUAUUGC  0.3667   2                   A:1-5(x5);A:20-23(x4)  gc_low;homopolymer  2
```

- **There is no pass/fail column.** `n_flags` 0 is the verdict, and both thresholds are
  echoed into every row (`gc_band`, `homopolymer_min`) so a flag can be audited from its
  own row without knowing how the CLI was invoked — the same convention as
  `target_site_qc.tsv`. `test_the_verdict_is_derivable_from_the_row` reproduces `n_flags`
  from the other columns.
- `homopolymer_runs` lists **every** run, not only the longest: a guide with three
  separate 4-runs is a different proposition from one with a single 4-run, and the
  `longest_*` columns cannot tell them apart. Format is
  `base:start-end(xlength)`, `;`-joined.
- Run offsets are **1-based inclusive and relative to the guide**, so
  `sequence[start - 1:end]` is the run. These are **not** transcript coordinates: the
  guide is the reverse complement of its target, so offset 1 in the guide is the *last*
  base of the target span named in its `tile_id`.
- Ambiguous bases (outside ACGTU) are excluded from the GC denominator and counted in
  `n_ambiguous` rather than scored as not-GC. They are **reported, not flagged** — an N
  is a property of the input FASTA, not of the guide's chemistry. The run summary notes
  them if any are present.

`--flagged-only` drops the passing rows. The default writes all of them, which is what
lets the file be joined onto `tiles.tsv`. Either way the stderr summary counts **every**
guide, not the written subset, so the passing fraction is still readable off a
`--flagged-only` run.

Hard errors: a missing or empty FASTA, a duplicate `guide_id` (the output is keyed on it,
so a duplicate would make a downstream join silently pick one of two rows — the same
reasoning as duplicate transcript ids in stage 1), an inverted GC band, and
`--homopolymer-min` below 2 (a run of 1 is every base, so every guide would be flagged).

### Finding: only a third of the PI3K guides pass both checks

Over the 28,032 guides in `results/pik3/guides.fasta`:

| | guides |
|---|---|
| passing both checks | 9,154 (32.7%) |
| `gc_low` (< 40%) | 11,700 |
| `gc_high` (> 60%) | 3,100 |
| `homopolymer` (>= 4) | 12,389 |

The flag counts do not sum to the flagged total (18,878) because a guide can raise two.
Verified against an independent oracle written without importing `seqqc` — plain
`itertools.groupby` for runs and a direct G/C count — with all four numbers matching
exactly.

Two things worth knowing before using this as a filter. **GC-low dominates**: the four
PI3K transcripts are AT-rich enough at 30 nt that nearly three times as many guides fall
below the band as above it, so tightening `--gc-min` costs far more candidates than
tightening `--gc-max`. And the **homopolymer flag is the single largest cause of
rejection** at 44% of guides — and it is not the `GGGG` the threshold is usually written
for. Of the 16,059 runs found, 12,463 (78%) are `A` or `U`:

| run base | runs |
|---|---|
| A | 6,252 |
| U | 6,211 |
| C | 1,906 |
| G | 1,690 |

So raising `--homopolymer-min` to 5 recovers mostly AT-rich guides, which are the same
ones `gc_low` is already rejecting. Whether the A/U runs are concentrated in the UTRs and
the PIK3CG tandem repeat has **not** been checked; the run offsets in
`homopolymer_runs` are guide-relative, so answering that means joining back onto
`tiles.tsv` for the transcript coordinate.

Combined with the stage 2 finding that all 28,032 tiles are `specific`, there is no
tension between the two filters here — every surviving guide is still single-gene. These
flags narrow the pool; they do not change its composition by gene.

### Dropping flagged guides from the plots

`scripts/plot_scatterlines.py --drop-qc-flagged` removes guides that fail the checks
above before plotting. The filter lives in
[src/visualise_scatterlines.py](src/visualise_scatterlines.py) as `drop_qc_flagged`,
returning a `QcDrop` with the surviving frame, per-flag counts and per-transcript drop
counts; [test/test_visualise_scatterlines.py](test/test_visualise_scatterlines.py) covers
it. Thresholds come from the same `seqqc` constants and are overridable with `--gc-min`,
`--gc-max` and `--homopolymer-min`.

```bash
python scripts/plot_scatterlines.py \
    --blastn results/blastn/siRNA_targetSeq.blastn \
    --predictions results/pik3/pik3-predictions.csv \
    --genbank-dir data/pik3-iso \
    --output-dir results/plots-qc/filtered \
    --drop-qc-flagged
```

**It is off by default.** The script's standing promise is that every row in the
predictions file is plotted, so filtering is opt-in and the stderr summary always states
what went, per flag and per transcript.

Two things make this filter safe to reuse rather than reimplement:

- It judges `spacer_seq` with `seqqc.check_guide`, **not** a join against
  `guide_qc.tsv`. The spacer is already in the predictions CSV, so there is one source of
  truth and no way to filter a plot against a stale QC run. Verified: over the real
  28,032 guides the filter drops exactly the 18,878 `guide_qc.tsv` flags, the same
  `tile_id` set, with identical per-flag counts.
- **`transcript_lengths` must be called on the unfiltered frame**, and the CLI does so
  before filtering. The axis length is the largest tile `end`, so taking it afterwards
  would stop each axis at the last surviving guide and make a QC-poor 3' end read as
  *absent* rather than as filtered out — the opposite of what drawing the full transcript
  is for. `TestFilteringDoesNotShrinkTheAxis` pins this by showing the two lengths differ
  (529 vs 30 on its fixture).

### Finding: QC filtering thins the validated siRNA sites unevenly

Filtering keeps 9,154 of 28,032 guides (32.7%), and 1,759 of the 4,820 guides above 0.9
predicted efficacy (36.5%) — so QC and efficacy are close to independent, and the filter
costs roughly two thirds of candidates at any efficacy cutoff.

The uneven part is coverage of the six siRNA sites that map into these transcripts. Each
site is overlapped by 48 guides before filtering:

| siRNA | transcript | site | surviving QC |
|---|---|---|---|
| PIK3CD_1_ss | NM_005026.5 | 543-561 | 39 / 48 (81%) |
| PIK3CD_2_ss | NM_005026.5 | 3254-3272 | 39 / 48 (81%) |
| PIK3CG_1_ss | NM_002649.3 | 1087-1105 | 30 / 48 (63%) |
| PIK3CD_3_ss | NM_005026.5 | 3286-3304 | 21 / 48 (44%) |
| PIK3CG_3_ss | NM_002649.3 | 2940-2958 | 17 / 48 (35%) |
| PIK3CG_2_ss | NM_002649.3 | 2035-2053 | 9 / 48 (19%) |

**No site is wiped out** — every one keeps at least 9 guides, which is enough to pick
from. But the range is 19% to 81%, so the filter is not a uniform thinning: a site's
local GC and homopolymer content decides how much choice is left there, and
`PIK3CG_2_ss` has about a fifth of the headroom `PIK3CD_1_ss` does. Worth checking
survival at a site before concluding a region is poorly served.

PIK3CA and PIK3CB have **no siRNA mapped** in `siRNA_targetSeq.blastn`, so this says
nothing about them; four further siRNAs in that file are PPIB and target a transcript
absent from the PI3K predictions.
