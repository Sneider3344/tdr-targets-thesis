#!/usr/bin/env python3
"""
Computes the reference pocket centroid for each of the 2 pilot OG7 groups,
from their crystallographic reference structure in BioLiP.

  OG7_0006581 -> mtub (ligand DKA, decanoic acid)
  OG7_0006003 -> ecol (ligand ABF, phosphorylated arabinofuranose)

Step 1: find the best-resolution BioLiP entry for the reference
species/ligand pair, and parse its binding-site residues.
Step 2: download the corresponding PDB structure and compute the centroid
of the Cα coordinates of those residues.
"""
import gzip
import os
from pathlib import Path

import numpy as np
import pandas as pd
import requests

COLS_BIOLIP = [
    "pdb_id", "chain", "resolution", "binding_site",
    "ligand_id", "ligand_chain", "ligand_serial",
    "binding_residues_pdb", "binding_residues_renum",
    "catalytic_pdb", "catalytic_renum",
    "ec_number", "go_terms",
    "affinity_manual", "affinity_moad", "affinity_pdbbind", "affinity_bindingdb",
    "uniprot_id", "pubmed_id", "ligand_seqnum", "receptor_seq",
]

BIOLIP_PATH = "/home/ssneider/disco-big/ssneider-env/TDR_Targets7.1/PocketVec/analisis_biolip/BioLiP.txt.gz"
PDB_DIR = "pdb_referencias"

# Best-resolution PDB structures already identified as the reference for
# each pilot group's binding site
REFERENCE_STRUCTURES = {
    "OG7_0006581": {"species": "mtub", "ligand": "DKA", "pdb_id": "1w66", "chain": "A"},
    "OG7_0006003": {"species": "ecol", "ligand": "ABF", "pdb_id": "1o8b", "chain": "A"},
}


def load_biolip():
    print("Loading BioLiP...")
    with gzip.open(BIOLIP_PATH, "rt", encoding="utf-8", errors="replace") as f:
        df_bio = pd.read_csv(f, sep="\t", header=None, names=COLS_BIOLIP, low_memory=False)
    df_bio["resolution"] = pd.to_numeric(df_bio["resolution"], errors="coerce")
    df_bio["uniprot_first"] = df_bio["uniprot_id"].astype(str).str.split(",").str[0].str.strip()
    df_bio["ligand_id"] = df_bio["ligand_id"].astype(str).str.strip().str.upper()
    print(f"BioLiP loaded: {len(df_bio):,} entries")
    return df_bio


def choose_best_structure(df_candidates, name):
    if len(df_candidates) == 0:
        print(f"ALERT: no structures found for {name}")
        return None
    best = df_candidates.sort_values("resolution").iloc[0]
    print(f"\n{name} - reference structure chosen:")
    print(f"  PDB ID: {best['pdb_id']}  Chain: {best['chain']}  Resolution: {best['resolution']} Å")
    print(f"  Pocket residues (PDB numbering): {best['binding_residues_pdb']}")
    print(f"  Pocket residues (renumbered from 1): {best['binding_residues_renum']}")
    return best


def parse_residues(binding_residues_str):
    """BioLiP format: 'N180 L181 A182 V215 H218 ...' (letter + number).
    Returns a list of (letter, number) tuples."""
    residues = []
    for r in str(binding_residues_str).split():
        residues.append((r[0], int(r[1:])))
    return residues


def download_pdb(pdb_id, outdir):
    outpath = Path(outdir) / f"{pdb_id.lower()}.pdb"
    if outpath.exists():
        print(f"  {pdb_id}: already exists, skipping download")
        return outpath
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    r = requests.get(url, timeout=30)
    if r.status_code == 200:
        outpath.write_text(r.text)
        print(f"  {pdb_id}: downloaded OK ({len(r.text):,} bytes)")
        return outpath
    print(f"  {pdb_id}: ERROR {r.status_code}")
    return None


def extract_ca_coords(pdb_path, chain, wanted_residue_numbers):
    """Returns {residue_number: (x, y, z)} for Cα atoms of the given chain
    and residue numbers."""
    coords = {}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            chain_id = line[21].strip()
            res_num = int(line[22:26].strip())
            if atom_name == "CA" and chain_id == chain and res_num in wanted_residue_numbers:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coords[res_num] = (x, y, z)
    return coords


def compute_centroid(coords_dict):
    if not coords_dict:
        return None
    points = np.array(list(coords_dict.values()))
    return points.mean(axis=0)


def main():
    df_bio = load_biolip()
    os.makedirs(PDB_DIR, exist_ok=True)

    df_uniprot = pd.read_csv("og7_uniprot_por_especie.tsv", sep="\t")
    centroid_rows = []

    for og7, ref in REFERENCE_STRUCTURES.items():
        # Step 1: find the exact BioLiP entry for this reference structure
        candidates = df_bio[(df_bio["ligand_id"] == ref["ligand"]) & (df_bio["resolution"] > 0)].copy()
        uniprot_ref = df_uniprot[
            (df_uniprot["og7"] == og7) & (df_uniprot["especie"] == ref["species"])
        ]["uniprot_id"].values
        candidates = candidates[candidates["uniprot_first"].isin(uniprot_ref)]
        best = choose_best_structure(candidates, f"{og7} ({ref['species']}, {ref['ligand']})")
        if best is None:
            continue

        residues_pdb = parse_residues(best["binding_residues_pdb"])
        print(f"Pocket residues for {og7} (PDB numbering): {len(residues_pdb)} residues")

        # Step 2: download the structure and compute the Cα centroid
        pdb_path = download_pdb(ref["pdb_id"], PDB_DIR)
        if pdb_path is None:
            continue

        wanted_numbers = {num for _, num in residues_pdb}
        coords = extract_ca_coords(pdb_path, ref["chain"], wanted_numbers)
        print(f"  Residues found: {len(coords)} of {len(residues_pdb)}")

        centroid = compute_centroid(coords)
        if centroid is None:
            print(f"  ALERT: could not compute a centroid for {og7}")
            continue
        print(f"  Centroid: x={centroid[0]:.4f}  y={centroid[1]:.4f}  z={centroid[2]:.4f}")

        centroid_rows.append({
            "og7": og7,
            "especie_ref": ref["species"],
            "pdb_id": ref["pdb_id"],
            "chain": ref["chain"],
            "ligando": ref["ligand"],
            "n_residuos_pocket": len(coords),
            "centroide_x": round(float(centroid[0]), 4),
            "centroide_y": round(float(centroid[1]), 4),
            "centroide_z": round(float(centroid[2]), 4),
            "residuos_usados": " ".join(f"{letter}{num}" for letter, num in residues_pdb if num in coords),
            "uniprot_referencia": uniprot_ref[0] if len(uniprot_ref) else None,
        })

    df_centroids = pd.DataFrame(centroid_rows)
    df_centroids.to_csv("pocket_centroides.tsv", sep="\t", index=False)
    print("\nSaved: pocket_centroides.tsv")
    print(df_centroids.to_string(index=False))


if __name__ == "__main__":
    main()
