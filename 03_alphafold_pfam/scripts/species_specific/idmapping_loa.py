"""
idmapping_loa.py

Loa loa — Ensembl Metazoa BioMart ID mapping pipeline.

Loa loa gene IDs are not directly available in the UniProt mapping tool.
The workaround uses the Ensembl Metazoa BioMart (https://metazoa.ensembl.org/biomart)
to download a table linking Gene stable IDs to Transcript stable IDs.
The Gene stable IDs are then used as input for the UniProt ID mapper.

This script joins the BioMart export with the resulting UniProt mapping,
producing a Transcript ID -> UniProt Entry table for downstream use.

Usage:
    python idmapping_loa.py \
        --biomart    mart_export_loa.txt \
        --idmapping  idmapping_loa.tsv \
        --output     idmapping_loa_final.tsv
"""

import argparse

import pandas as pd


def map_loa_ids(biomart_path: str, idmapping_path: str, output_path: str) -> None:
    genes_transcripts = pd.read_csv(biomart_path, sep="\t")
    gene_uniprot      = pd.read_csv(idmapping_path, sep="\t")

    merged = genes_transcripts.merge(
        gene_uniprot,
        left_on="Gene stable ID",
        right_on="From",
        how="left",
    )

    final = (
        merged[["Transcript stable ID", "Entry"]]
        .rename(columns={"Transcript stable ID": "From"})
        .dropna(subset=["Entry"])
    )

    final.to_csv(output_path, sep="\t", index=False)
    print(f"Mapped transcripts : {len(final)}")
    print(f"Output             : {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Loa loa BioMart-based ID mapping.")
    parser.add_argument("--biomart",    required=True, help="BioMart export (Gene -> Transcript)")
    parser.add_argument("--idmapping",  required=True, help="UniProt ID mapping TSV (From, Entry)")
    parser.add_argument("--output",     required=True)
    args = parser.parse_args()

    map_loa_ids(args.biomart, args.idmapping, args.output)


if __name__ == "__main__":
    main()
