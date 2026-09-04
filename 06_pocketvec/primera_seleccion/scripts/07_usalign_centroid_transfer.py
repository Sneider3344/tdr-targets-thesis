#!/usr/bin/env python3
"""
Downloads each species' AlphaFold model, structurally aligns it against
the reference crystal structure with US-align, and transfers the pocket
centroid into the AlphaFold model's coordinate frame using the resulting
rotation matrix and translation vector.
"""
import os
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import requests

USALIGN     = "/big/USalign_/USalign"
BASE_DIR    = "/big/lab/ssneider/ssneider-env/TDR_Targets7.1/PocketVec/pLDDT_candidatos"
PDB_REF_DIR = f"{BASE_DIR}/pdb_referencias"
AF_PDB_DIR  = f"{BASE_DIR}/alphafold_pdbs"
USALIGN_DIR = f"{BASE_DIR}/usalign_outputs"

SPECIES_OF_INTEREST = ["hsap", "dmel", "cele", "scer", "kpm", "loa", "egr", "mtub"]
OG7_REFERENCE = {
    "OG7_0006581": {"pdb": "1w66", "chain": "A", "reference_species": "mtub",
                     "centroid": np.array([8.6440, 9.6843, 5.3440]),
                     "residues_pdb": [76, 78, 79, 145, 158, 159, 176]},
    "OG7_0006003": {"pdb": "1o8b", "chain": "A", "reference_species": "ecol",
                     "centroid": np.array([-1.3073, -4.9566, 1.3497]),
                     "residues_pdb": [28, 30, 31, 81, 84, 94, 95, 96, 103, 121]},
}


def download_alphafold_pdb(uniprot_id, species, outdir):
    outpath = Path(outdir) / f"AF-{uniprot_id}-F1_{species}.pdb"
    if outpath.exists():
        print(f"  {uniprot_id} ({species}): already exists")
        return outpath
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    r_api = requests.get(api_url, timeout=30)
    if r_api.status_code != 200:
        print(f"  {uniprot_id} ({species}): API error {r_api.status_code}")
        return None
    data = r_api.json()
    if not data:
        print(f"  {uniprot_id} ({species}): no data in API response")
        return None
    pdb_url = data[0].get("pdbUrl")
    if not pdb_url:
        print(f"  {uniprot_id} ({species}): no pdbUrl in API response")
        return None
    r = requests.get(pdb_url, timeout=60)
    if r.status_code == 200:
        outpath.write_text(r.text)
        print(f"  {uniprot_id} ({species}): downloaded OK from {pdb_url}")
        return outpath
    print(f"  {uniprot_id} ({species}): ERROR downloading PDB {r.status_code}")
    return None


def run_usalign(pdb_ref, chain_ref, pdb_mobile, matrix_path):
    """Aligns pdb_mobile against pdb_ref, saving the rotation matrix with -m."""
    cmd = [USALIGN, str(pdb_mobile), str(pdb_ref), "-chain2", chain_ref, "-m", str(matrix_path)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.stdout, result.stderr


def parse_matrix_file(matrix_path):
    """Reads the rotation-matrix file produced by US-align with -m.
    Format (0-indexed rows):
      0   t0   u00  u01  u02
      1   t1   u10  u11  u12
      2   t2   u20  u21  u22
    Returns t (3,) and U (3x3), or (None, None) if not found/parseable."""
    t, U = [], []
    try:
        with open(matrix_path) as f:
            for line in f:
                s = line.strip()
                if s.startswith(("0", "1", "2")):
                    parts = s.split()
                    if len(parts) >= 5:
                        t.append(float(parts[1]))
                        U.append([float(parts[2]), float(parts[3]), float(parts[4])])
    except FileNotFoundError:
        return None, None
    if len(t) == 3 and len(U) == 3:
        return np.array(t), np.array(U)
    return None, None


def parse_tm_score(stdout):
    for line in stdout.split("\n"):
        if "TM-score=" in line and "normalized by length of Structure_1" in line:
            m = re.search(r"TM-score=\s*([\d.]+)", line)
            if m:
                return float(m.group(1))
    return None


def main():
    os.makedirs(AF_PDB_DIR, exist_ok=True)
    os.makedirs(USALIGN_DIR, exist_ok=True)

    df_uniprot = pd.read_csv("og7_uniprot_por_especie.tsv", sep="\t")
    df_uniprot = df_uniprot[df_uniprot["especie"].isin(SPECIES_OF_INTEREST)]
    print("Available UniProt IDs:")
    print(df_uniprot.to_string(index=False))

    print("\nDownloading AlphaFold PDBs...")
    af_pdb_paths = {}
    for _, row in df_uniprot.iterrows():
        if pd.isna(row["uniprot_id"]):
            continue
        key = (row["og7"], row["especie"])
        path = download_alphafold_pdb(row["uniprot_id"], row["especie"], AF_PDB_DIR)
        if path:
            af_pdb_paths[key] = path
    print(f"\nAlphaFold PDBs downloaded: {len(af_pdb_paths)}")

    print("\nRunning US-align...")
    alignment_results = []
    for og7, ref in OG7_REFERENCE.items():
        pdb_ref = Path(PDB_REF_DIR) / f"{ref['pdb'].lower()}.pdb"
        reference_species = ref["reference_species"]
        reference_centroid = ref["centroid"]

        for species in SPECIES_OF_INTEREST:
            key = (og7, species)

            if species == reference_species:
                print(f"  {og7} | {species}: is the reference, centroid = {reference_centroid}")
                alignment_results.append({
                    "og7": og7, "especie": species, "es_referencia": True, "tm_score": 1.0,
                    "centroide_transferido_x": reference_centroid[0],
                    "centroide_transferido_y": reference_centroid[1],
                    "centroide_transferido_z": reference_centroid[2],
                })
                continue

            if key not in af_pdb_paths:
                print(f"  {og7} | {species}: no AlphaFold PDB, skipping")
                alignment_results.append({
                    "og7": og7, "especie": species, "es_referencia": False, "tm_score": None,
                    "centroide_transferido_x": None, "centroide_transferido_y": None,
                    "centroide_transferido_z": None,
                })
                continue

            pdb_mobile = af_pdb_paths[key]
            matrix_path = Path(USALIGN_DIR) / f"{og7}_{species}_matrix.txt"
            stdout, stderr = run_usalign(pdb_ref, ref["chain"], pdb_mobile, matrix_path)

            tm_score = parse_tm_score(stdout)
            t, U = parse_matrix_file(matrix_path)

            if t is None:
                print(f"  {og7} | {species}: could not parse the US-align matrix")
                alignment_results.append({
                    "og7": og7, "especie": species, "es_referencia": False, "tm_score": tm_score,
                    "centroide_transferido_x": None, "centroide_transferido_y": None,
                    "centroide_transferido_z": None,
                })
                continue

            # Transform the centroid: af_centroid = U @ ref_centroid + t
            af_centroid = U @ reference_centroid + t
            print(f"  {og7} | {species}: TM-score={tm_score:.3f} | "
                  f"transferred centroid = ({af_centroid[0]:.3f}, {af_centroid[1]:.3f}, {af_centroid[2]:.3f})")

            alignment_results.append({
                "og7": og7, "especie": species, "es_referencia": False, "tm_score": tm_score,
                "centroide_transferido_x": round(float(af_centroid[0]), 4),
                "centroide_transferido_y": round(float(af_centroid[1]), 4),
                "centroide_transferido_z": round(float(af_centroid[2]), 4),
            })

    df_alignment = pd.DataFrame(alignment_results)
    df_alignment = df_alignment.merge(df_uniprot[["og7", "especie", "uniprot_id"]],
                                       on=["og7", "especie"], how="left")
    df_alignment.to_csv("centroides_transferidos.tsv", sep="\t", index=False)

    print("\n── ALIGNMENT SUMMARY ────────────────────────────────────────────────")
    print(df_alignment[["og7", "especie", "uniprot_id", "tm_score",
                          "centroide_transferido_x", "centroide_transferido_y",
                          "centroide_transferido_z"]].to_string(index=False))

    print("\nSpecies with TM-score < 0.5 (weak alignment, unreliable centroid):")
    low = df_alignment[df_alignment["tm_score"] < 0.5]
    print("  None -- all alignments are reliable" if len(low) == 0
          else low[["og7", "especie", "tm_score"]].to_string(index=False))


if __name__ == "__main__":
    main()
