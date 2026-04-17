"""
04_merge_all_species.py

Concatenates per-species confidence score tables (produced by
03_calculate_confidence_scores.py) into a single file for cross-species
PFAM-level analysis.

Usage:
    python 04_merge_all_species.py \
        --pattern "Species/*/*_confidenceScore.txt" \
        --output  all_species_confidence.tsv
"""

import argparse
import glob

import pandas as pd


def merge_tables(pattern: str, output_path: str) -> None:
    files = glob.glob(pattern)
    if not files:
        raise FileNotFoundError(f"No files matched pattern: {pattern}")

    dfs = [pd.read_csv(f, sep="\t") for f in files]
    combined = pd.concat(dfs, ignore_index=True)
    combined.to_csv(output_path, sep="\t", index=False)

    print(f"Files merged     : {len(files)}")
    print(f"Total rows       : {len(combined)}")
    print(f"Output           : {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge per-species AlphaFold confidence tables.")
    parser.add_argument("--pattern", required=True, help="Glob pattern for input files")
    parser.add_argument("--output",  required=True)
    args = parser.parse_args()

    merge_tables(args.pattern, args.output)


if __name__ == "__main__":
    main()
