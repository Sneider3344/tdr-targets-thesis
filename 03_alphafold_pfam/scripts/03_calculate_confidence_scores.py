"""
03_calculate_confidence_scores.py

For each protein domain (PFAM) in the joined table, reads the corresponding
AlphaFold confidence JSON and computes the mean per-residue confidence score
over the domain's annotated start-end range.

Confidence scores are stored as per-residue arrays in the JSON field
'confidenceScore'. Indexing is 1-based in the annotation file, 0-based
in the JSON array.

Usage:
    python 03_calculate_confidence_scores.py \
        --table         <species_table_joined.txt> \
        --alphafold_dir <folder_with_json_files> \
        --output        <species_confidenceScore.txt>
"""

import argparse
import json
from pathlib import Path

import pandas as pd


def compute_confidence(table_path: Path, json_dir: Path, output_path: Path) -> None:
    df = pd.read_csv(table_path, sep="\t")
    df["Confidence_Avg"] = None

    missing = 0
    for i, row in df.iterrows():
        json_file = json_dir / row["AlphaFold_ID"]

        if not json_file.exists():
            print(f"[WARN] Not found: {json_file.name}")
            missing += 1
            continue

        with open(json_file) as f:
            data = json.load(f)

        start  = int(row["Inicio"])
        end    = int(row["Fin"])
        scores = data["confidenceScore"][start - 1 : end]   # convert to 0-based

        if scores:
            df.at[i, "Confidence_Avg"] = round(sum(scores) / len(scores), 2)

    df.to_csv(output_path, sep="\t", index=False)
    filled = df["Confidence_Avg"].notna().sum()
    print(f"Domains scored   : {filled}")
    print(f"Missing JSONs    : {missing}")
    print(f"Output           : {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Calculate AlphaFold confidence scores per PFAM domain.")
    parser.add_argument("--table",         required=True)
    parser.add_argument("--alphafold_dir", required=True)
    parser.add_argument("--output",        required=True)
    args = parser.parse_args()

    compute_confidence(
        Path(args.table),
        Path(args.alphafold_dir),
        Path(args.output),
    )


if __name__ == "__main__":
    main()
