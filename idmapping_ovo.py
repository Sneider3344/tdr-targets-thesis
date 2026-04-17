"""
idmapping_ovo.py

Onchocerca volvulus — two-step ID mapping pipeline.

UniProt does not directly map OVOC gene IDs. Instead, the Ensembl Metazoa
GFF3 annotation file provides OVOC -> WBGene correspondences, and the
WBGene IDs can then be mapped to UniProt via the standard UniProt ID mapper.

Step 1 (run_gff3_to_wbgene): Parse the Ensembl GFF3, extract OVOC -> WBGene pairs.
Step 2 (run_wbgene_to_uniprot): Join the WBGene mapping with the UniProt TSV.

Download GFF3 from:
  https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-61/gff3/onchocerca_volvulus/

Usage:
    python idmapping_ovo.py gff3 \
        --gff3      Onchocerca_volvulus.ASM49940v2.61.gff3.gz \
        --gene_list ovo_gene_name.txt \
        --output    ovoc_with_wbgene.csv

    python idmapping_ovo.py uniprot \
        --wbgene_csv    ovoc_with_wbgene.csv \
        --idmapping_tsv idmapping_ovo.tsv \
        --output        idmapping_ovo_final.tsv
"""

import argparse
import gzip
import re

import pandas as pd


def run_gff3_to_wbgene(gff3_file: str, gene_list_file: str, output_csv: str) -> None:
    """Parse GFF3 and map OVOC IDs to WBGene IDs."""
    ovoc_to_wbgene: dict[str, str] = {}

    with gzip.open(gff3_file, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue
            attrs   = cols[8]
            m_ovoc  = re.search(r"Name=(OVOC\d+)", attrs)
            m_wb    = re.search(r"gene_id=(WBGene\d+)", attrs)
            if m_ovoc and m_wb:
                ovoc_id = m_ovoc.group(1)
                if ovoc_id not in ovoc_to_wbgene:
                    ovoc_to_wbgene[ovoc_id] = m_wb.group(1)

    rows = []
    with open(gene_list_file) as f:
        for line in f:
            ovoc_raw   = line.strip()
            if not ovoc_raw:
                continue
            ovoc_clean = re.sub(r"\.\d+$", "", ovoc_raw)   # strip isoform suffix
            wbgene     = ovoc_to_wbgene.get(ovoc_clean, "NA")
            rows.append({"OVOC": ovoc_raw, "WBGene": wbgene})

    df = pd.DataFrame(rows)
    df.to_csv(output_csv, index=False)
    na_count = (df["WBGene"] == "NA").sum()
    print(f"OVOC IDs processed : {len(df)}")
    print(f"Not mapped (NA)    : {na_count}")
    print(f"Output             : {output_csv}")


def run_wbgene_to_uniprot(wbgene_csv: str, idmapping_tsv: str, output_tsv: str) -> None:
    """Join WBGene IDs with UniProt entries and write final mapping."""
    ovoc_wbgene    = pd.read_csv(wbgene_csv)
    wbgene_uniprot = pd.read_csv(idmapping_tsv, sep="\t")

    # Normalise column names and whitespace
    ovoc_wbgene.columns    = ovoc_wbgene.columns.str.strip()
    wbgene_uniprot.columns = wbgene_uniprot.columns.str.strip()
    ovoc_wbgene    = ovoc_wbgene.apply(lambda c: c.str.strip() if c.dtype == object else c)
    wbgene_uniprot = wbgene_uniprot.apply(lambda c: c.str.strip() if c.dtype == object else c)

    wbgene_uniprot = wbgene_uniprot.rename(columns={"From": "WBGene", "Entry": "UniProt"})
    ovoc_wbgene    = ovoc_wbgene.dropna(subset=["WBGene"])

    merged = ovoc_wbgene.merge(wbgene_uniprot, on="WBGene", how="left")
    final  = (
        merged[["OVOC", "UniProt"]]
        .rename(columns={"OVOC": "From", "UniProt": "Entry"})
        .dropna(subset=["Entry"])
    )

    final.to_csv(output_tsv, sep="\t", index=False)
    print(f"Mapped entries     : {len(final)}")
    print(f"Output             : {output_tsv}")


def main() -> None:
    parser = argparse.ArgumentParser(description="OVO two-step ID mapping pipeline.")
    sub = parser.add_subparsers(dest="step", required=True)

    p1 = sub.add_parser("gff3", help="Step 1: GFF3 -> WBGene CSV")
    p1.add_argument("--gff3",      required=True)
    p1.add_argument("--gene_list", required=True)
    p1.add_argument("--output",    required=True)

    p2 = sub.add_parser("uniprot", help="Step 2: WBGene CSV + UniProt TSV -> final mapping")
    p2.add_argument("--wbgene_csv",    required=True)
    p2.add_argument("--idmapping_tsv", required=True)
    p2.add_argument("--output",        required=True)

    args = parser.parse_args()

    if args.step == "gff3":
        run_gff3_to_wbgene(args.gff3, args.gene_list, args.output)
    else:
        run_wbgene_to_uniprot(args.wbgene_csv, args.idmapping_tsv, args.output)


if __name__ == "__main__":
    main()
