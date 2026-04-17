"""
03_compare_orthogroup_versions.py

Compares how sequences were grouped in OrthoMCL v6 versus v7.
For every sequence present in both versions, checks whether the numeric
orthogroup identifier is conserved across versions.  Results are split into:
  - <species>_orthomcl_iguales.tsv   : sequences that kept the same group
  - <species>_orthomcl_distintos.tsv : sequences reassigned to a different group

The numeric identifier is extracted from the full group ID string using
version-aware offsets (see extract_group_number docstring).

Usage:
    python 03_compare_orthogroup_versions.py \
        --old  <v6_grouped.tsv> \
        --new  <v7_grouped.tsv> \
        --same <iguales.tsv> \
        --diff <distintos.tsv>
"""

import argparse
import pandas as pd


def extract_group_number(group_id: str, version: int) -> str | None:
    """
    Pull the numeric core from an OrthoMCL group identifier.

    v6 IDs look like:  OG6_1_00042  -> skip first char after '_', strip zeros
    v7 IDs look like:  OG7_10_00042 -> skip first two chars after '_', strip zeros

    Returns None for IDs that don't match the expected format.
    """
    try:
        suffix = group_id.split("_")[1]
        offset = 1 if version == 6 else 2
        return suffix[offset:].lstrip("0") or "0"
    except (IndexError, AttributeError):
        return None


def expand_group_rows(df: pd.DataFrame, group_col: str, version: int) -> pd.DataFrame:
    """Explode bracketed sequence lists into individual rows and add numeric group key."""
    rows = []
    for _, row in df.iterrows():
        seq_ids = row["Query sequence id"].strip("[]").replace(" ", "").split(",")
        num = extract_group_number(row[group_col], version)
        for seq_id in seq_ids:
            rows.append({"Query sequence id": seq_id, group_col: row[group_col], "group_number": num})
    return pd.DataFrame(rows).dropna(subset=["group_number"])


def main(file_old: str, file_new: str, output_same: str, output_diff: str) -> None:
    old_df = pd.read_csv(file_old, sep="\t", header=0,
                         names=["OrthoMCL6 group id", "Query sequence id"])
    new_df = pd.read_csv(file_new, sep="\t", header=0,
                         names=["OrthoMCL7 group id", "Query sequence id"])

    old_expanded = expand_group_rows(old_df, "OrthoMCL6 group id", version=6)
    new_expanded = expand_group_rows(new_df, "OrthoMCL7 group id", version=7)

    same, diff = [], []

    for _, old_row in old_expanded.iterrows():
        seq_id = old_row["Query sequence id"]
        matches = new_expanded[new_expanded["Query sequence id"] == seq_id]

        for _, new_row in matches.iterrows():
            record = [seq_id, old_row["OrthoMCL6 group id"], new_row["OrthoMCL7 group id"]]
            if old_row["group_number"] == new_row["group_number"]:
                same.append(record)
            else:
                diff.append(record)

    cols = ["Query sequence id", "OrthoMCL6 group id", "OrthoMCL7 group id"]
    pd.DataFrame(same, columns=cols).to_csv(output_same, sep="\t", index=False)
    pd.DataFrame(diff, columns=cols).to_csv(output_diff, sep="\t", index=False)

    total = len(same) + len(diff)
    pct_same = (len(same) / total * 100) if total else 0
    print(f"Conserved groupings     : {len(same)}  ({pct_same:.1f}%)")
    print(f"Reassigned groupings    : {len(diff)}  ({100 - pct_same:.1f}%)")
    print(f"Conserved output        : {output_same}")
    print(f"Reassigned output       : {output_diff}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare OrthoMCL v6 vs v7 groupings.")
    parser.add_argument("--old",  required=True, help="Grouped v6 TSV")
    parser.add_argument("--new",  required=True, help="Grouped v7 TSV")
    parser.add_argument("--same", required=True, help="Output: conserved groupings TSV")
    parser.add_argument("--diff", required=True, help="Output: reassigned groupings TSV")
    args = parser.parse_args()
    main(args.old, args.new, args.same, args.diff)
