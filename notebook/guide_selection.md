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

