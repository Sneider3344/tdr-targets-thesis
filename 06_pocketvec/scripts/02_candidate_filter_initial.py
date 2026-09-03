#!/usr/bin/env python3
"""
Initial candidate filter ("filtrado inicial" in the thesis).

Starting from resumen_por_og7.tsv (OG7 groups with at least one BioLiP
crystal structure), applies three simultaneous filters: mean resolution
< 2.5 Å, main ligand other than RNA/DNA/peptide, and presence in at least
3 species -- counting only species with their own crystal structure.

This is the first-pass count; it was later corrected in
03_candidate_filter_corrected_coverage.py to count species present in the
network as a whole, not just species with their own crystal.
"""
import pandas as pd

OUT_DIR = "biolip_exploracion"

LIGANDS_EXCLUDE_ANALYSIS = {
    "RNA", "DNA", "PEPTIDE", "OHX", "A", "U", "G", "C",
    "LLL", "PAR", "GLC", "PTD",
}


def ligand_ok(top_ligand_str):
    """True if the most frequent ligand is NOT RNA/DNA/peptide."""
    if pd.isna(top_ligand_str):
        return False
    first_ligand = str(top_ligand_str).split(",")[0].strip().upper()
    return first_ligand not in LIGANDS_EXCLUDE_ANALYSIS


def main():
    og = pd.read_csv(f"{OUT_DIR}/resumen_por_og7.tsv", sep="\t")
    print(f"Total OG7 groups in BioLiP: {len(og):,}")

    og["ligand_ok"] = og["top_ligandos"].apply(ligand_ok)
    og["resolution_ok"] = og["resolucion_media"] < 2.5
    og["species_ok"] = og["n_especies"] >= 3

    og_filtered = og[og["ligand_ok"] & og["resolution_ok"] & og["species_ok"]].copy()
    og_filtered = og_filtered.sort_values(["n_especies", "uniprot_unicos"], ascending=False)

    print(f"Groups passing the filters (res<2.5Å, biological ligand, >=3 species): {len(og_filtered):,}")

    cols_show = ["og7", "n_especies", "uniprot_unicos", "pdb_unicos",
                 "resolucion_media", "top_ligandos", "ec_numbers", "especies_lista"]

    print("\n── TOP 20 CANDIDATES FOR POCKETVEC ─────────────────────────────────────")
    print(og_filtered[cols_show].head(20).to_string(index=False))

    og_filtered[cols_show].to_csv(f"{OUT_DIR}/candidatos_pocketvec.tsv", sep="\t", index=False)
    print(f"\nSaved: {OUT_DIR}/candidatos_pocketvec.tsv")

    print("\n── TOP 10 SUMMARY ────────────────────────────────────────────────────")
    for _, row in og_filtered.head(10).iterrows():
        print(f"\n{row['og7']}")
        print(f"  Species ({int(row['n_especies'])}): {row['especies_lista']}")
        print(f"  Unique UniProt: {int(row['uniprot_unicos'])} | Unique PDB: {int(row['pdb_unicos'])}")
        print(f"  Mean resolution: {row['resolucion_media']:.2f} Å")
        print(f"  Top ligands: {row['top_ligandos']}")
        print(f"  EC: {row['ec_numbers']}")


if __name__ == "__main__":
    main()

## Column reference
## og7 — OrthoMCL 7 group identifier. Arbitrary, says nothing about function by itself.
## n_especies — how many of the 29 species have at least one protein of this group with a BioLiP structure. The most important criterion: groups present across many species make cross-species pocket comparison meaningful.
## uniprot_unicos — how many distinct proteins of the group have a BioLiP structure. A group can have many proteins across many species but few with a crystallized structure -- this tells you how many are actually usable.
## pdb_unicos — how many distinct PDB structures. The same protein can have dozens of crystallized structures with different ligands or conditions -- always >= uniprot_unicos.
## resolucion_media — mean crystallographic resolution across the group's structures. Kept below 2.5 Å for PocketVec.
## top_ligandos — the 3 most frequent ligands in the group's structures. RNA/DNA/PEPTIDE = discarded. ATP/ADP/GTP/SAM = interesting.
## ec_numbers — EC numbers of the group's proteins. Groups with no EC are non-enzymatic (structural, RNA-binding, etc.).
## especies_lista — which specific species have structures. Useful to check whether it includes the parasites of interest (lmaj, ldon, pfal, tcru, tgon, tbrt) or only model organisms.
