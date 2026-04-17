"""
idmapping_tcru.py

Trypanosoma cruzi dm28c — UniProt ID prioritisation pipeline.

TriTrypDB gene entries can map to multiple UniProt IDs for the same protein.
This script applies a ranked selection strategy to choose the single best
UniProt ID per gene, which is then used to locate the corresponding
AlphaFold confidence JSON:

  Priority 1 — Curated UniProt IDs (those that do NOT start with 'A0A').
               These are manually reviewed entries and are preferred.
  Priority 2 — If only 'A0A' IDs exist, select the one whose AlphaFold JSON
               is largest on disk. Larger files correlate with longer, better-
               characterised sequences.
  Priority 3 — If no JSON exists for any candidate, fall back to the first
               ID in the list.

Usage:
    python idmapping_tcru.py \
        --input         idmapping_tcru.txt \
        --alphafold_dir /path/to/alphafold/tcru \
        --output        dm28c_mapped.tsv
"""

import argparse
import os

import pandas as pd


def pick_best_uniprot(uniprot_field: str, json_dir: str) -> str | None:
    """Return the best UniProt ID for a given gene according to the priority rules."""
    if not isinstance(uniprot_field, str):
        return None
    field = uniprot_field.strip()
    if not field or field.upper() == "N/A":
        return None

    ids = [x.strip() for x in field.split(",") if x.strip()]
    if not ids:
        return None

    # Priority 1: curated (non-A0A) IDs
    non_a0a = [u for u in ids if not u.upper().startswith("A0A")]
    candidates = non_a0a if non_a0a else ids

    # Priority 2: largest AlphaFold JSON among candidates
    best_id  = None
    best_len = -1
    for uid in candidates:
        json_path = os.path.join(json_dir, f"{uid}.json")
        if os.path.isfile(json_path):
            try:
                size = os.path.getsize(json_path)
                if size > best_len:
                    best_len = size
                    best_id  = uid
            except OSError:
                continue

    # Priority 3: first candidate if no JSON was found
    if best_id is None:
        best_id = candidates[0]

    return best_id


def main() -> None:
    parser = argparse.ArgumentParser(description="Prioritise UniProt IDs for T. cruzi dm28c.")
    parser.add_argument("--input",         required=True, help="TSV with 'Gene ID' and 'UniProt ID(s)' columns")
    parser.add_argument("--alphafold_dir", required=True)
    parser.add_argument("--output",        required=True)
    args = parser.parse_args()

    print(f"Reading: {args.input}")
    df = pd.read_csv(args.input, sep="\t", dtype=str)

    required = {"Gene ID", "UniProt ID(s)"}
    missing  = required - set(df.columns)
    if missing:
        raise SystemExit(f"Missing columns: {missing}. Found: {list(df.columns)}")

    out = df[["Gene ID", "UniProt ID(s)"]].rename(columns={"Gene ID": "From", "UniProt ID(s)": "Entry"}).copy()
    out["Entry"] = out["Entry"].apply(lambda x: pick_best_uniprot(x, args.alphafold_dir))
    out = out.dropna(subset=["Entry"])

    out.to_csv(args.output, sep="\t", index=False)
    print(f"Genes mapped     : {len(out)}")
    print(f"Output           : {args.output}")


if __name__ == "__main__":
    main()
