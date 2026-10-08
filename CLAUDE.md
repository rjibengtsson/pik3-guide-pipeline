# pik3-guide-pipeline

Pipeline for designing CRISPR-Cas13 guide RNAs that target a single gene isoform, a
defined subset, or all members of a gene family. The worked example is the class IA/IB
PI3K catalytic subunits (PIK3CA/B/D/G).

Stages 1 (tiling + guide generation) and 2 (family classification) exist. Remaining
stages — GC/MFE filtering, genome-wide off-target mapping — are not written yet.

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

# tests (116, all passing; needs bwa and blastn, else some skip)
python -m unittest discover -s test -v

# stage 1, exactly as the committed results/pik3/ were produced
python scripts/tile_transcripts.py --input-dir data/pik3-iso \
    --output results/pik3/tiles.tsv --input-type mrna -k 30 --step 1 \
    --per-transcript --fasta-out results/pik3/guides.fasta
```

## Layout

- [src/tile.py](src/tile.py) — the library: FASTA parsing, tiling, the guide
  convention. Holds the module-level docstring that defines the coordinate,
  direction and alphabet conventions; read it before changing anything here.
- [scripts/tile_transcripts.py](scripts/tile_transcripts.py) — the CLI around it:
  argument parsing, row assembly, TSV/FASTA writing, the stderr run summary.
- [test/test_tile.py](test/test_tile.py) — unittest suite covering both.
- [data/pik3-iso/](data/pik3-iso/) — input RefSeq mRNA FASTA, one transcript per
  file, one folder per dataset.
- [results/pik3/](results/pik3/) — committed output for that dataset (see **Outputs**).

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

`results/pik3/` holds the committed run over the four PI3K transcripts: 28,032 tiles
at k=30, step=1.

- `tiles.tsv` — combined, columns `tile_id, transcript_id, source_file, start, end,
  strand, length, target_seq, guide_rna, has_ambiguous`.
- `tiles/{transcript_id}.tsv` — same columns, one file per transcript (`--per-transcript`).
- `guides.fasta` — guide sequences, headers set to `tile_id`.
- `tiles.fasta` — the same tiles written with `--fasta-field target_seq`, i.e. targets
  rather than guides. The two FASTA files are **not** interchangeable.

Output order is deterministic: input files sorted by path, records in file order, tiles
by ascending start.

## Known rough edges

- **`.gitignore` tries to ignore `./data` and `./results` but does not.** Git ignore
  patterns cannot begin with `./`, so both lines are inert and both directories are
  tracked — `results/` is ~170k lines of the repo. If the intent was to untrack them,
  the patterns need to be `data/` and `results/`; if the intent is to keep them
  committed, the lines should go.
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
