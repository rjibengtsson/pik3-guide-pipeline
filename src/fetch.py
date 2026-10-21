"""Download RNA (transcript) and protein sequences for an NCBI Gene ID.

Requires the NCBI Datasets CLI ('datasets') on your PATH.
"""
from __future__ import annotations

from pathlib import Path
from Bio import Entrez

Entrez.email = ""


def download_gene_id(gene_id: int, out_file: str) -> Path:
    """Download transcripts and proteins for one NCBI Gene ID.

    Parameters
    ----------
    gene_id : int
        NCBI Gene ID (e.g. 672 for human BRCA1).
    out_dir : str
        Directory where the zip and extracted files are written.

    Returns
    -------
    Path
        Directory containing the extracted files.
    """

    # 1. Link the gene ID to its transcript (RNA) records
    handle = Entrez.elink(dbfrom="gene", db="nuccore", id=gene_id,
                        linkname="gene_nuccore_refseqrna")
    links = Entrez.read(handle)
    handle.close()
    rna_ids = [l["Id"] for l in links[0]["LinkSetDb"][0]["Link"]]

    # 2. Fetch the transcripts as FASTA
    handle = Entrez.efetch(db="nuccore", id=rna_ids, 
                           rettype="fasta", retmode="text")
    text = handle.read()
    handle.close()

    with open(out_file, "w") as f:
        f.write(text)

    return Path(out_file)
