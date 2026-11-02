# Guide Selection

This notebook outlines the process for selecting guides for PIK3 knockdown experiments.

# Coverage analysis
To narrow down regions of interest to target with guides, perform coverage analysis on the PIK3 isoforms. For each transcript, download all avaliable variants from NCBI and align them. Then, select regions of high conservation for each isoform. These regions will be used to design guides that target all variants of a given isoform.

```bash
# Download all variants for each isoform
python3 scripts/fetch_genes.py -id 5294 -o data/pik3-iso/PIK3CG.fasta
python3 scripts/fetch_genes.py -id 5293 -o data/pik3-iso/PIK3CD.fasta
python3 scripts/fetch_genes.py -id 5290 -o data/pik3-iso/PIK3CA.fasta
python3 scripts/fetch_genes.py -id 5291 -o data/pik3-iso/PIK3CB.fasta
```
Peform multiple sequence alignment across all variants for each isoform using `mafft`

```bash
python3 src/alignment.py -f data/pik3-iso/PIK3CG.fasta -o data/pik3-iso/PIK3CG_aligned.fasta
python3 src/alignment.py -f data/pik3-iso/PIK3CD.fasta -o data/pik3-iso/PIK3CD_aligned.fasta
python3 src/alignment.py -f data/pik3-iso/PIK3CA.fasta -o data/pik3-iso/PIK3CA_aligned.fasta
python3 src/alignment.py -f data/pik3-iso/PIK3CB.fasta -o data/pik3-iso/PIK3CB_aligned.fasta
```
Highly conserved regions were selected for each isoform by visual inspection of the multiple sequence alignment in UGENE. Alignment indices were then converted to reference genome coordinates using the scripts below.

```bash
python3 scripts/get_conserved_regions.py -f data/pik3-iso/PIK3CG_aligned.fasta -r NM_002649.3 -i 341,2982 -o PIK3CG_location.blastn

python3 scripts/get_conserved_regions.py -f data/pik3-iso/PIK3CD_aligned.fasta -b data/pik3-iso/PIK3CD_aligned.bed -o PIK3CD_location.blastn

python3 scripts/get_conserved_regions.py -f data/pik3-iso/PIK3CA_aligned.fasta -r NM_006218.4 -i 18165,27176 -o PIK3CA_location.blastn

python3 scripts/get_conserved_regions.py -f data/pik3-iso/PIK3CB_aligned.fasta -b data/pik3-iso/PIK3CB_aligned.bed -o PIK3CB_location.blastn
```

# QC filter

Guides were filtered based on GC content and homopolymer runs. The following thresholds were used:

| flag | raised when | rationale |
|---|---|---|
| `gc_low` | GC fraction below `--gc-min` | weak guide-target duplex |
| `gc_high` | GC fraction above `--gc-max` | secondary structure, reduced specificity |
| `homopolymer_t` | run of >=4 T/U (`--homopolymer-t-min`) | U6 premature-termination signal |
| `homopolymer_v` | run of >=4 A, C or G (`--homopolymer-v-min`) | synthesis and folding problems |

The rationale column is the conventional reasoning for these thresholds, not something
measured in this pipeline — nothing here tests activity, only sequence composition.

GC thresholds are **fractions, not percentages** — `0.40` means 40%. The two run flags
are kept separate because they are different liabilities with different remedies: a T/U
run terminates transcription, whereas an A/C/G run is a synthesis problem. A combined
count cannot say which one a guide has. `V` is IUPAC for A, C or G.

Both the measurements and the thresholds are written into every output row, so a flag can
be audited from its own row without knowing how the script was invoked.

The two runs below cover the **four PI3K transcripts** (28,032 guides); PPIB is scanned
the same way further down, and behaves quite differently.

## PI3K run 1: default band, 40-60% GC

```bash
python scripts/check_guides.py --guides results/pik3/guides.fasta \
    --output results/pik3/guide_qc.tsv
```

| | guides | of 28,032 |
|---|---|---|
| **passing all checks** | **9,154** | **32.7%** |
| `gc_low` (< 40%) | 11,700 | 41.7% |
| `gc_high` (> 60%) | 3,100 | 11.1% |
| `homopolymer_t` | 5,606 | 20.0% |
| `homopolymer_v` | 8,427 | 30.1% |

Flags overlap, so these do not sum to the 18,878 guides flagged; 1,644 guides raise both
run flags.

Two things stand out. **GC-low dominates** — the four PI3K transcripts are AT-rich enough
at 30 nt that nearly four times as many guides fall below the band as above it, so
tightening `--gc-min` costs far more candidates than tightening `--gc-max`. And
`homopolymer_v` is **driven by poly-A, not poly-G**: of the 16,059 runs found, A accounts
for 6,252 and G for only 1,690 (U 6,211, C 1,906). So most of `homopolymer_v` is AT-rich
sequence rather than the G-rich liability that flag is usually reached for.

## PI3K run 2: wide band, 10-90% GC

To separate what GC filtering costs from what the run filters cost, the same scan was
repeated with GC effectively disabled:

```bash
python scripts/check_guides.py --guides results/pik3/guides.fasta \
    --gc-min 0.1 --gc-max 0.9 \
    --output results/pik3/guide_qc_gc10-90_t4_v4.tsv
```

| | run 1 (0.40-0.60) | run 2 (0.10-0.90) |
|---|---|---|
| **passing** | 9,154 (32.7%) | **15,637 (55.8%)** |
| `gc_low` | 11,700 | 62 |
| `gc_high` | 3,100 | 15 |
| `homopolymer_t` | 5,606 | 5,606 |
| `homopolymer_v` | 8,427 | 8,427 |

The run-flag counts are **identical** between the two, which is the check that widening
the GC band did not leak into the homopolymer scans.

At 10-90%, GC is no longer doing real work: only 77 guides fall outside the band, because
the observed GC range across all 28,032 guides is 0.033 to 0.967. Of the 12,395 guides
still flagged, 12,389 are caught by a run flag — so **only 6 guides are rejected on GC
alone**. This is effectively a homopolymer-only filter, and passing is capped at 55.8% by
the run flags no matter how wide GC is opened.

So the ~23 percentage points between the two runs is the cost of the 40-60% band
specifically, and the remaining 44% rejection is structural: it cannot be recovered by
relaxing GC, only by relaxing the homopolymer thresholds or accepting the runs.

## Output naming

QC outputs are named with their parameters, since the thresholds are what make two runs
different:

```
guide_qc.tsv                      defaults: GC 0.40-0.60, T>=4, V>=4
guide_qc_gc{min}-{max}_t{N}_v{N}  explicit, GC as percent to avoid dots in the filename
guide_qc_gc10-90_t4_v4.tsv        e.g. GC 0.10-0.90, T>=4, V>=4
```

## Carrying the filter into the plots

`scripts/plot_scatterlines.py --drop-qc-flagged` drops flagged guides from the efficacy
scatters. It recomputes QC from the `spacer_seq` column with the same `seqqc.check_guide`
rather than joining `guide_qc.tsv`, so a plot can never be filtered against a stale QC
run. Thresholds are set with the same flag names.

```bash
python scripts/plot_scatterlines.py \
    --blastn results/blastn/siRNA_targetSeq.blastn \
    --predictions results/pik3/pik3-predictions.csv \
    --genbank-dir data/pik3-iso \
    --output-dir results/plots-qc/filtered \
    --drop-qc-flagged
```

The x-axis is unaffected by filtering — transcript length is taken before the filter — so
a QC-poor stretch reads as *filtered out* rather than as uncovered.

At the default band, coverage of the six siRNA sites mapped to the **PI3K** transcripts
survives **unevenly**. Each site is overlapped by 48 guides before filtering:

| siRNA | transcript | site | surviving QC |
|---|---|---|---|
| PIK3CD_1_ss | NM_005026.5 | 543-561 | 39 / 48 (81%) |
| PIK3CD_2_ss | NM_005026.5 | 3254-3272 | 39 / 48 (81%) |
| PIK3CG_1_ss | NM_002649.3 | 1087-1105 | 30 / 48 (63%) |
| PIK3CD_3_ss | NM_005026.5 | 3286-3304 | 21 / 48 (44%) |
| PIK3CG_3_ss | NM_002649.3 | 2940-2958 | 17 / 48 (35%) |
| PIK3CG_2_ss | NM_002649.3 | 2035-2053 | 9 / 48 (19%) |

No site is wiped out — every one keeps at least 9 guides — but the range is 19% to 81%,
so the filter is not a uniform thinning. Local GC and run content decides how much choice
is left at each site, and `PIK3CG_2_ss` has about a fifth of the headroom `PIK3CD_1_ss`
does. Worth checking survival at a site before concluding a region is poorly served.

PIK3CA and PIK3CB have no siRNA mapped in `siRNA_targetSeq.blastn`, so this says nothing
about them.

# QC filter: PPIB (NM_000942.5)

PPIB was scanned with the same two bands, as a single-gene control. It is a much smaller
set — 893 nt giving **864 guides** at k=30, step=1 — and the guide FASTA had to be
generated first, since `results/ppib/` held only `tiles.tsv` and the predictions:

```bash
# guides.fasta for ppib (tiles.tsv already existed and was left alone)
python scripts/tile_transcripts.py --input-dir data/ppib \
    --output /tmp/ppib-tiles-check.tsv --input-type mrna -k 30 --step 1 \
    --fasta-out results/ppib/guides.fasta

python scripts/check_guides.py --guides results/ppib/guides.fasta \
    --output results/ppib/guide_qc.tsv

python scripts/check_guides.py --guides results/ppib/guides.fasta \
    --gc-min 0.1 --gc-max 0.9 \
    --output results/ppib/guide_qc_gc10-90_t4_v4.tsv
```

The re-tiling was written to a scratch path and diffed against the committed
`results/ppib/tiles.tsv` — identical — and `guides.fasta` was checked to match the
`spacer_seq` column of `ppib-predictions.csv` tile-for-tile, so the QC applies to exactly
the guides the predictions describe.

| | run 1 (0.40-0.60) | run 2 (0.10-0.90) |
|---|---|---|
| **passing** | 367 (42.5%) | **522 (60.4%)** |
| `gc_low` | 82 | 0 |
| `gc_high` | 235 | 0 |
| `homopolymer_t` | 85 | 85 |
| `homopolymer_v` | 284 | 284 |

## PPIB is GC-rich, so the band fails it from the opposite side

This is the main result, and it inverts the PI3K finding:

| | PI3K (4 transcripts) | PPIB |
|---|---|---|
| guides | 28,032 | 864 |
| mean GC | 0.427 | **0.545** |
| GC range | 0.033-0.967 | 0.133-0.767 |
| `gc_low` : `gc_high` | 11,700 : 3,100 | **82 : 235** |
| passing at 40-60% | 32.7% | 42.5% |

PI3K loses guides for being **too AT-rich**; PPIB loses nearly three times as many to
`gc_high` as to `gc_low`. So a GC band tuned on one of these datasets is tuned in the
wrong direction for the other — tightening `--gc-min` costs PI3K dearly and PPIB almost
nothing, and `--gc-max` the reverse.

The run composition inverts too. PI3K's runs are overwhelmingly A/U (12,463 of 16,059,
78%); PPIB's 369 runs are led by **G** and are nearly flat:

| run base | PI3K | PPIB |
|---|---|---|
| A | 6,252 | 95 |
| U | 6,211 | 85 |
| C | 1,906 | 81 |
| G | 1,690 | **108** |

So `homopolymer_v` means different things in the two datasets: mostly poly-A in PI3K,
genuinely G-led in PPIB. Only 27 PPIB guides raise both run flags.

## At 10-90%, GC filters nothing at all for PPIB

PPIB's observed GC range (0.133-0.767) sits **entirely inside** the wide band, so run 2
flags zero guides on GC — not merely few, as with PI3K's 77. Run 2 is therefore a pure
homopolymer-only filter here, and the 342 guides it drops are exactly the 342 raising a
run flag.

That makes the split clean: of the 57.5% rejected at the default band, **17.9 percentage
points are the GC band** (42.5% -> 60.4%) and the remaining 39.6% is structural, from the
run flags alone.

## Plots

```bash
python scripts/plot_scatterlines.py \
    --blastn results/blastn/siRNA_targetSeq.blastn \
    --predictions results/ppib/ppib-predictions.csv \
    --genbank-dir data/ppib \
    --output-dir results/plots-qc/filtered \
    --drop-qc-flagged
```

Written to `results/plots-qc/{filtered,unfiltered}/NM_000942.5_PPIB_anno.html`. PPIB has
its own four siRNAs in `siRNA_targetSeq.blastn` (`ppib_1`-`ppib_4`). Note that table also
holds spurious **7 bp** hits of the PIK3 siRNAs against PPIB; these lose on bitscore to
the real 19 bp hits, so the best-HSP-per-siRNA rule assigns each siRNA to the right
transcript — worth knowing before reading that table by hand.

| siRNA | site | surviving QC |
|---|---|---|
| ppib_4 | 232-250 | 17 / 48 (35%) |
| ppib_1 | 305-323 | 16 / 48 (33%) |
| ppib_2 | 380-398 | 16 / 48 (33%) |
| ppib_3 | 285-303 | 13 / 48 (27%) |

Unlike PI3K, where site survival ranged from 19% to 81%, PPIB's four sites are tightly
clustered at **27-35%** — and all four sit *below* the 42.5% genome-wide pass rate, so the
siRNA-validated regions are slightly GC-richer and run-heavier than PPIB average.

**Efficacy and QC pull against each other here.** Only 38 of PPIB's 129 guides above 0.9
predicted efficacy survive QC (29.5%), against a 42.5% overall pass rate — so the
high-efficacy guides are *disproportionately* QC-poor. PI3K shows the opposite, mild
trend (36.5% of its >0.9 guides survive, against 32.7% overall). With 38 guides left
above 0.9, PPIB still has workable candidates, but the margin is much thinner than the
raw count suggests.
