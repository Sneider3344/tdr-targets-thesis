"""
05_pfam_aggregate_stats.py

Groups the merged all-species table by PFAM domain and computes:
  - mean confidence score across all species and proteins
  - total annotation count
  - number of species in which the domain appears

Results are sorted by mean confidence score (descending).

Usage:
    python 05_pfam_aggregate_stats.py \
        --input  all_species_confidence.tsv \
        --output PFAM_confidence_avg.tsv
"""

import argparse

import pandas as pd


def compute_pfam_stats(input_path: str, output_path: str) -> None:
    df = pd.read_csv(input_path, sep="\t")

    stats = (
        df.groupby("PFAM")
        .agg(
            Confidence_Avg=("Confidence_Avg", "mean"),
            Count=("PFAM", "count"),
            Species=("especie", "nunique"),
        )
        .reset_index()
    )
    stats["Confidence_Avg"] = stats["Confidence_Avg"].round(2)
    stats = stats.sort_values("Confidence_Avg", ascending=False)

    stats.to_csv(output_path, sep="\t", index=False)
    print(f"Unique PFAMs     : {len(stats)}")
    print(f"Output           : {output_path}")
    print("\nTop 10 by confidence:")
    print(stats.head(10).to_string(index=False))


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute per-PFAM confidence score statistics.")
    parser.add_argument("--input",  required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    compute_pfam_stats(args.input, args.output)


if __name__ == "__main__":
    main()
