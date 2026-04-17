"""
01_merge_orthomcl_versions.py

Merges OrthoMCL v6 and v7 TSV outputs for a given species into a single
combined table, keyed on Query sequence ID.
Outputs matched rows and basic coverage statistics to stdout.

Usage:
    python 01_merge_orthomcl_versions.py --old <v6.tsv> --new <v7.tsv> --out <combined.tsv>
"""

import argparse
import pandas as pd


def load_orthomcl_table(path: str, group_col_name: str) -> pd.DataFrame:
    """Read an OrthoMCL TSV and return columns 0 and 2 with clean IDs."""
    df = pd.read_csv(
        path,
        sep="\t",
        usecols=[0, 2],
        names=["Query sequence id", group_col_name],
        header=0,
    )
    # Strip pipe-prefixed tags (e.g. 'kpm|KPM_0001' -> 'KPM_0001')
    df["Query sequence id"] = (
        df["Query sequence id"]
        .str.replace(r"^.*\|", "", regex=True)
        .str.strip()
        .astype(str)
    )
    return df


def main(file_old: str, file_new: str, output_file: str) -> None:
    old_df = load_orthomcl_table(file_old, "OrthoMCL6 group id")
    new_df = load_orthomcl_table(file_new, "OrthoMCL7 group id")

    merged_df = old_df.merge(new_df, on="Query sequence id", how="inner")
    merged_df.to_csv(output_file, sep="\t", index=False)

    print(f"Sequences in v6 table  : {len(old_df)}")
    print(f"Sequences in v7 table  : {len(new_df)}")
    print(f"Matched sequences      : {len(merged_df)}")
    print(f"Output written to      : {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge OrthoMCL v6 and v7 tables.")
    parser.add_argument("--old", required=True, help="Path to OrthoMCL v6 TSV")
    parser.add_argument("--new", required=True, help="Path to OrthoMCL v7 TSV")
    parser.add_argument("--out", required=True, help="Path for combined output TSV")
    args = parser.parse_args()
    main(args.old, args.new, args.out)
