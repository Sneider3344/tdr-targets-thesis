"""
02_group_by_orthogroup.py

For a given OrthoMCL output file, collapses individual sequence rows into
one row per orthogroup, listing all member Query sequence IDs as a bracketed
comma-separated string.  Sequences that could not be assigned to any group
(group id == -1) are saved to a separate file for manual inspection.

Usage:
    python 02_group_by_orthogroup.py \
        --input  <orthomcl_output.tsv> \
        --out    <grouped.tsv> \
        --ungrouped <nongrouped.tsv> \
        --prefix kpm
"""

import argparse
import pandas as pd


def strip_species_prefix(series: pd.Series, prefix: str) -> pd.Series:
    """Remove 'species|' tag from sequence IDs (e.g. 'kpm|KPM_0001' -> 'KPM_0001')."""
    return series.str.replace(rf"^{prefix}\|", "", regex=True)


def main(input_file: str, output_file: str, ungrouped_file: str, prefix: str) -> None:
    df = pd.read_csv(input_file, sep="\t", header=0)

    group_col = df.columns[2]
    seq_col   = df.columns[0]

    # Split into grouped and ungrouped entries
    ungrouped = df[df[group_col] == "-1"].copy()
    grouped   = df[df[group_col] != "-1"].copy()

    # Clean up species prefix from sequence IDs
    grouped[seq_col]   = strip_species_prefix(grouped[seq_col], prefix)
    ungrouped[seq_col] = strip_species_prefix(ungrouped[seq_col], prefix)

    # One row per orthogroup: aggregate member IDs into a bracketed list
    result = (
        grouped
        .groupby(group_col)[seq_col]
        .apply(lambda ids: "[" + ", ".join(ids) + "]")
        .reset_index()
    )
    result.columns = ["OrthoMCL group id", "Query sequence id"]

    result.to_csv(output_file, sep="\t", index=False)
    ungrouped.to_csv(ungrouped_file, sep="\t", index=False)

    print(f"Orthogroups written     : {len(result)}")
    print(f"Ungrouped sequences     : {len(ungrouped)}")
    print(f"Grouped output          : {output_file}")
    print(f"Ungrouped output        : {ungrouped_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Group OrthoMCL sequences by orthogroup.")
    parser.add_argument("--input",     required=True, help="Raw OrthoMCL TSV")
    parser.add_argument("--out",       required=True, help="Path for grouped output TSV")
    parser.add_argument("--ungrouped", required=True, help="Path for ungrouped output TSV")
    parser.add_argument("--prefix",    required=True, help="Species prefix to strip (e.g. kpm)")
    args = parser.parse_args()
    main(args.input, args.out, args.ungrouped, args.prefix)
