# pik3-guide-pipeline
Pipeline for designing CRISPR-Cas13 guide RNAs that target a single gene isoform, a defined subset, or all members of a gene family.

## Environment

FASTA parsing uses `Bio.SeqIO`, which needs Biopython **and** numpy (Biopython
imports TwoBitIO unconditionally, so a missing numpy breaks plain FASTA too).
The base conda env has no numpy, so use the project env:

```bash
mamba env create -f environment.yml   # or: conda env create -f environment.yml
conda activate pik3-guide
```

Importing `src/tile.py` outside that env raises an ImportError saying exactly
this, rather than Biopython's misleading "TwoBit files" message.

## Stage 1 — tiling and guide generation

Break every input sequence in a folder into fixed-length k-mers, one tile per
starting position, tracking the source transcript and coordinates, then reverse
complement each tile into a Cas13 guide.

`--input-type` is required and must be `dna` or `mrna`. It is taken on trust
and recorded in the run summary — the script never inspects the sequences to
work out which they are, and there is no interactive prompt.

Because `target_seq` is verbatim and `guide_rna` is derived with Biopython's
`reverse_complement_rna` (which accepts both alphabets), the two choices
currently produce identical output. The flag documents what you fed in; it does
not change the result.

```bash
python scripts/tile_transcripts.py --input-dir data --output results/tiles.tsv \
    --input-type mrna -k 30 --step 1 --per-transcript --fasta-out results/guides.fasta
```

Output columns: `tile_id, transcript_id, source_file, start, end, strand,
length, target_seq, guide_rna, has_ambiguous`.

- `target_seq` is the tile **verbatim from the input**, in whatever alphabet the
  FASTA used — a DNA-alphabet input keeps its `T`. Nothing is transcribed on the
  target side.
- Positions are **1-based inclusive**, so `input_sequence[start - 1:end]`
  recovers `target_seq` exactly.
- `tile_id` (`NM_006218.4:1-30`) is the join key for downstream stages.
- `guide_rna` **is** transcribed, because a guide has to base-pair as RNA: it is
  the reverse complement of the target in the RNA alphabet, so A in the target
  gives U in the guide.
- **Both columns are written 5' -> 3'**, the direction you order an oligo in.
  Guide and target are antiparallel, so the guide's 5' end pairs with the
  target's **3'** end:

  ```
  target  5'-AAAAACCCCCGGGGGTTTTTAAGGCCTTAC-3'   (as in the input)
  guide   5'-GUAAGGCCUUAAAAACCCCCGGGGGUUUUU-3'   (RNA)
  ```

  The test suite's `pairs_antiparallel` helper asserts this on every emitted
  row, computed independently of `reverse_complement` so it fails rather than
  agrees if that function is wrong. It compares in RNA space, so a DNA-alphabet
  target is fine.
- mRNA is single-stranded, so `strand` is always `+`.
- `--guide-alphabet dna` writes the guide with T instead of U, which is what
  bowtie/blast expect for off-target mapping. `target_seq` is unaffected.
- Tiles containing ambiguous bases are kept and flagged; `--drop-ambiguous`
  omits them.
- `--fasta-out` writes the guides by default; `--fasta-field target_seq`
  writes the targets instead.

Caveats — note that none of these are checked for you:

- The declared `--input-type` is not verified against the sequences. Nothing
  inspects a sequence to decide whether it is DNA or RNA.
- RefSeq mRNA records are written in the DNA alphabet (they contain `T`).
  Answering "mRNA" for them is correct: `target_seq` keeps the `T`, and the
  T -> U substitution is applied only to the guide.
- "DNA" here means the coding/sense strand of the transcript. **Genomic** DNA
  with introns, or an antisense strand, is not handled — splice and orient it
  before tiling.
- A sequence mixing `T` and `U`, or holding anything other than nucleotides,
  is passed straight through. Only `has_ambiguous` (non-ACGTU bases per tile)
  is reported.

Tests (inside the `pik3-guide` env): `python -m unittest discover -s test -v`
