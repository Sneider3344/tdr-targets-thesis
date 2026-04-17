"""
07_filter_by_pfam_list.py

Filters the merged all-species confidence table to keep only rows whose
PFAM domain appears in a user-supplied list of PFAMs of interest (e.g.
those associated with known drug targets).

Usage:
    python 07_filter_by_pfam_list.py \
        --all_species  all_species_confidence.tsv \
        --pfam_list    PFAM_list_clean.tsv \
        --output       PFAM_filtered.tsv
"""

import argparse

import pandas as pd


def filter_by_pfam_list(all_species_path: str, pfam_list_path: str, output_path: str) -> None:
    df = pd.read_csv(all_species_path, sep="\t")

    pfam_df   = pd.read_csv(pfam_list_path, sep="\t", header=None)
    pfam_list = pfam_df.iloc[:, 0].dropna().unique().tolist()

    filtered = df[df["PFAM"].isin(pfam_list)]
    filtered.to_csv(output_path, sep="\t", index=False)

    print(f"PFAMs in list    : {len(pfam_list)}")
    print(f"Matching rows    : {len(filtered)}")
    print(f"Output           : {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Filter confidence table by PFAM list of interest.")
    parser.add_argument("--all_species", required=True)
    parser.add_argument("--pfam_list",   required=True)
    parser.add_argument("--output",      required=True)
    args = parser.parse_args()

    filter_by_pfam_list(args.all_species, args.pfam_list, args.output)


if __name__ == "__main__":
    main()
