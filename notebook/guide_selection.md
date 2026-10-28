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

For each conserved region, select the 1.5–2 kb sub-region with the widest coverage of high-efficacy guides. Selected sub-regions will be cloned and tested.

To determin the best sub-region, will use a sliding window scan and score each window by the number of hotspots containing high-efficacy guides. Hotspots are defined as clusters of overlapping high-efficacy guides, and the window with the highest number of hotspots will be selected for cloning.

# Hotspot window scan

Slide a window along each conserved region and score it by the guides inside it.

A good 30-mer site yields ~30 overlapping high-efficacy guides and a guide count would rank the
longest good stretch first rather than the region offering the most independent
choices. Overlapping high-efficacy guides are merged into **hotspots** — one
hotspot is one site — and a window is scored by how many whole hotspots it
contains.

```bash
python3 scripts/find_hotspot_regions.py \
    --predictions results/pik3/pik3-predictions.csv \
    --regions results/blastn/PIK3C*_location.blastn \
    --output results/pik3/hotspot_regions.tsv \
    --all-windows results/pik3/hotspot_windows.tsv \
    --hotspots-out results/pik3/hotspots.tsv \
    --max-gap 60 --sweep-thresholds 0.90,0.93,0.95
```

## Dense and evenly covered, not a few top scorers

The selection criterion is **many independent sites spread evenly**, not the
highest-scoring guides. Hotspot count gives the first half of that. `--max-gap`
gives the second: it rejects any window containing a stretch of N+ bases with no
hotspot in it, because a construct is only as useful as its worst region and a
dead zone is where no guide can be placed at all. Without it, `largest_gap` is
only a tie-break and a clumpy window can still win on count.

This matters on real data. PIK3CA's two candidate blocks are near-tied on raw
guide count but not on shape:

| PIK3CA candidate | hotspots | guides >=0.9 | guides >=0.95 | worst dead zone | core covered |
|---|---|---|---|---|---|
| 2048-4047 | **23** | 261 | 14 | **35 bp** | **81.0%** |
| 5498-7497 | 17 | **319** | **39** | 115 bp | 74.8% |

The second block has *more* guides at every cutoff and a higher ceiling, but
packs them into fewer, longer stretches with a 115 bp hole. `--max-gap 60`
excludes it.

## Choosing the two parameters that matter

`--min-efficacy` is **absolute, not a quantile**, so counts stay comparable
between transcripts — which means it has to be set against the score
distribution of the predictor, not by intuition. For `pik3-predictions.csv` the
median predicted efficacy is already **0.823**, so a 0.8 cutoff qualifies over
half of every transcript, hotspots merge into kilobase-wide blobs, and a 2 kb
window contains 0–3 of them. The default **0.9** (about the 90th percentile
here) gives hotspots 60–90 bp wide and 11–23 per 2 kb window.

`--merge-gap` — how many low-scoring bases may separate two clusters of one
hotspot — is the most sensitive knob in the script. At `--min-efficacy 0.9`:

| `--merge-gap` | PIK3CG | PIK3CD | PIK3CA | PIK3CB |
|---|---|---|---|---|
| 0 (overlap or abut) | 58 | 47 | 87 | 48 |
| 10 | 34 | 23 | 56 | 30 |
| 30 | 8 | 10 | 22 | 15 |
| 100 | 1 | 1 | 2 | 1 |

The default **0** merges only clusters that overlap or abut, so two guides
sharing no base stay two independent sites. Both parameters are echoed into the
run summary, since a hotspot count means nothing without them.

Window length is swept over 1.5–2 kb rather than fixed at 2 kb, because three of
the PIK3CD conserved regions are shorter than 2 kb. A longer window cannot
contain fewer hotspots, so ties break towards the **shorter** window: the answer
is the smallest construct achieving the best count.

## Result over the committed PI3K run

240 hotspots at `--min-efficacy 0.9 --merge-gap 0`; 55 (23%) fall outside the
conserved regions or straddle a boundary and so are not placeable in any window.
`--max-gap 60` then rejects 617 of the 747 candidate windows.

| transcript | gene | best window | length | hotspots | worst dead zone | guides |
|---|---|---|---|---|---|---|
| NM_006218.4 | PIK3CA | 2048-4047 | 2000 | 23 of 87 | 35 bp | 261 |
| NM_006219.3 | PIK3CB | 1987-3486 | 1500 | 15 of 48 | 36 bp | 228 |
| NM_005026.5 | PIK3CD | 3215-5214 | 2000 | 20 of 47 | 39 bp | 377 |
| NM_002649.3 | PIK3CG | 349-2348 | 2000 | 16 of 58 | 55 bp | 375 |

### How robust each of these is

Tested by re-running at `--merge-gap` 0, 5, 10 and 20 with everything else
fixed. Three of the four winners hold; one does not.

| gene | winner across merge-gap 0/5/10/20 | verdict |
|---|---|---|
| PIK3CG | 349-2348 at all four | solid |
| PIK3CB | 1987-3486 at all four | solid |
| PIK3CD | 3215-5214, then 2515-4514 | shifts 700 bp, still overlapping |
| PIK3CA | 2048-4047 → 2298-4297 → 3098-5097 → 6548-8297 | **drifts, then jumps** |

**PIK3CA is genuinely undetermined by this metric.** Its early block (~2048-4047)
and late block (~6000-8300) are near-tied, and which wins depends on the
parameters rather than the data — the same ambiguity shows up under the threshold
sweep. Break that tie on something else: ORF position, predicted secondary
structure, or the genome-wide off-target stage. The other three can be taken as
they stand.

### Two things not to over-read

**The gain over an average window is modest.** Hotspots are fairly evenly
distributed, so the best window beats the worst by roughly 1.5x, not 10x.

**`--max-gap 60` and a high `--min-efficacy` are incompatible here.** Above about
0.9 the hotspots thin out and no window has a 60 bp core, so the threshold sweep
returns no row for most transcripts and reports stability as UNTESTED rather than
confirmed. That is why robustness above is argued from `--merge-gap` instead. If
you want a threshold sweep to mean something, loosen `--max-gap` for it.

`hotspot_windows.tsv` holds every scored window *before* the `--max-gap` filter
and overlap suppression — deliberately, since a window rejected by the filter is
exactly what you want to see when a winner looks surprising. `hotspots.tsv` has
every hotspot with its guide count and best guide.

## Target region selection

Since the cluster of hotspot window scan analysis revealed that clusters are fairly 
even distributed across and gain over different windows is modest, we will go ahead 
and select regions within CDS of each isoform and will select the following regions for cloning:

| transcript | gene | best window | length |
|---|---|---|---|---|---|---|
| NM_006218.4 | PIK3CA | 2048-4047 | 2000 | 
| NM_006219.3 | PIK3CB | 2000-4000 | 2000 | 
| NM_005026.5 | PIK3CD | 3215-5214 | 2000 | 
| NM_002649.3 | PIK3CG | 349-2348 | 2000 |

## Checking for propertoes of selected regions for cloning

Several properties of the selected regions were checked to ensure they are suitable for cloning. The following checks were performed:

1. **GC Content**: The GC content of each selected region was calculated to ensure it falls within an optimal range for cloning and guide design.
2. **Poly-T stretches**: scanned for long stretches of thymine (T) bases.
3. ***Repeats and homopolymers**: scanned for the presence of repetitive sequences.
4. **Reading frame and stop codons**: checked for the presence of stop codons and to ensure they are in-frame with the coding sequence of the gene.
