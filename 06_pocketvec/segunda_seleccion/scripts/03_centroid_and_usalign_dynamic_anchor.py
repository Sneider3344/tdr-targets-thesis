#!/usr/bin/env python3
"""
Centroid computation and structural alignment for the 3 finalist OG7
groups (OG7_0000433/ADP, OG7_0001020/RAP, OG7_0006854/FMN).

Unlike the first-stage pilot (one fixed crystallized species as the sole
anchor per group), this round has several crystallized species per group.
For each species without its own crystal structure, US-align is run
against every crystallized species in the group, and the crystallized
species with the highest TM-score is used as that species' anchor for
the centroid transfer -- rather than a single group-wide anchor.

GRUPOS is the per-group, per-species configuration table: which
structure to use (crystal PDB ID/chain, or AlphaFold + UniProt ID) for
each species. It's built here from BioLiP once (for the crystallized
species) and then reused by every later script in this stage; it's also
checkpointed to a pickle so long-running steps don't need to be redone.
"""
import pickle
import re
import subprocess
from pathlib import Path

import gzip
import numpy as np
import pandas as pd
import requests

BASE_DIR   = Path("/home/ssneider/disco-big/ssneider-env/TDR_Targets7.1/PocketVec")
ETAPA3_DIR = BASE_DIR / "Nueva_estrategia_pocketvec" / "etapa3_pocketvec"
USALIGN    = "/big/USalign_/USalign"

ALPHAFOLD_DIR = ETAPA3_DIR / "alphafold_pdbs"
CRYSTAL_DIR   = ETAPA3_DIR / "pdb_cristalizados"
USALIGN_DIR   = ETAPA3_DIR / "usalign_outputs"
GRUPOS_PICKLE = ETAPA3_DIR / "GRUPOS_checkpoint.pkl"

COLS_BIOLIP = [
    "pdb_id", "chain", "resolution", "binding_site",
    "ligand_id", "ligand_chain", "ligand_serial",
    "binding_residues_pdb", "binding_residues_renum",
    "catalytic_pdb", "catalytic_renum",
    "ec_number", "go_terms",
    "affinity_manual", "affinity_moad", "affinity_pdbbind", "affinity_bindingdb",
    "uniprot_id", "pubmed_id", "ligand_seqnum", "receptor_seq",
]

# Species/structure configuration per finalist group. Structural (PDB IDs,
# UniProt IDs) rather than derived from a script, since this is the fixed
# outcome of the candidate-selection stages above.
GRUPOS = {
    "OG7_0000433": {
        "ligando": "ADP",
        "especies": {
            "hsap": {"modo": "cristal", "pdb_id": "4zg4", "chain": "B"},
            "atha": {"modo": "cristal", "pdb_id": "7dhw", "chain": "A"},
            "ddis": {"modo": "cristal", "pdb_id": "1w9i", "chain": "A"},
            "pfal": {"modo": "cristal", "pdb_id": "8a12", "chain": "A"},
            "cele": {"modo": "cristal", "pdb_id": "6qdj", "chain": "A"},
            "scer": {"modo": "alphafold", "uniprot_id": "P08964"},
            "ehia": {"modo": "alphafold", "uniprot_id": "N9UR41"},
            "tcru": {"modo": "alphafold", "uniprot_id": "A0A2V2VHE6"},
            "egr":  {"modo": "alphafold", "uniprot_id": "A0A068WYU9"},
        },
    },
    "OG7_0001020": {
        "ligando": "RAP",
        "especies": {
            "hsap": {"modo": "cristal", "pdb_id": "4dri", "chain": "A"},
            "atha": {"modo": "cristal", "pdb_id": "7f2j", "chain": "A"},
            "pfal": {"modo": "cristal", "pdb_id": "4qt3", "chain": "A"},
            "ecol": {"modo": "alphafold", "uniprot_id": "P45523"},
            "scer": {"modo": "alphafold", "uniprot_id": "P38911"},
            "tcru": {"modo": "alphafold", "uniprot_id": "A0A2V2WP18"},
            "cele": {"modo": "alphafold", "uniprot_id": "O45418"},
            "egr":  {"modo": "alphafold", "uniprot_id": "A0A068X276"},
        },
    },
    "OG7_0006854": {
        "ligando": "FMN",
        "especies": {
            "hsap": {"modo": "cristal", "pdb_id": "4oqv", "chain": "A"},
            "lmaj": {"modo": "cristal", "pdb_id": "4xq6", "chain": "A"},
            "pfal": {"modo": "cristal", "pdb_id": "7l01", "chain": "A"},
            "mtub": {"modo": "cristal", "pdb_id": "3tq0", "chain": "A"},
            "ecol": {"modo": "cristal", "pdb_id": "7t5k", "chain": "A"},
            "dmel": {"modo": "alphafold", "uniprot_id": "Q9VLM9"},
            "cele": {"modo": "alphafold", "uniprot_id": "Q9XW01"},
            "calb": {"modo": "alphafold", "uniprot_id": "Q874I4"},
        },
    },
}


def parse_pdb_residues(binding_str):
    residues = []
    for token in str(binding_str).split():
        m = re.match(r"^([A-Za-z])(\d+)$", token)
        if m:
            residues.append((m.group(1), int(m.group(2))))
    return residues


def download_pdb(pdb_id, outdir):
    outpath = Path(outdir) / f"{pdb_id.lower()}.pdb"
    if outpath.exists():
        return outpath
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    r = requests.get(url, timeout=30)
    if r.status_code == 200:
        outpath.write_text(r.text)
        return outpath
    print(f"  ERROR downloading {pdb_id}: {r.status_code}")
    return None


def download_alphafold_pdb(uniprot_id, species, outdir):
    outpath = Path(outdir) / f"AF-{uniprot_id}-F1_{species}.pdb"
    if outpath.exists():
        return outpath
    r_api = requests.get(f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}", timeout=30)
    if r_api.status_code != 200:
        return None
    data = r_api.json()
    if not data:
        return None
    pdb_url = data[0].get("pdbUrl")
    if not pdb_url:
        return None
    r = requests.get(pdb_url, timeout=60)
    if r.status_code == 200:
        outpath.write_text(r.text)
        return outpath
    return None


def extract_ca_coords(pdb_path, chain, residue_numbers):
    wanted = set(residue_numbers)
    coords = {}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() == "CA" and line[21].strip() == chain:
                res_num = int(line[22:26].strip())
                if res_num in wanted:
                    coords[res_num] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    return coords


def run_usalign(pdb_mobile, chain_mobile, pdb_target, chain_target, matrix_path):
    cmd = [USALIGN, str(pdb_mobile), str(pdb_target), "-chain1", chain_mobile,
           "-chain2", chain_target, "-m", str(matrix_path)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.stdout


def parse_matrix_file(matrix_path):
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


def parse_tm_scores_and_rmsd(stdout):
    """US-align with default (asymmetric) output reports two TM-scores
    (normalized by each structure's length) and one RMSD."""
    tm1 = tm2 = rmsd = None
    for line in stdout.split("\n"):
        if "TM-score=" in line and "normalized by length of Structure_1" in line:
            m = re.search(r"TM-score=\s*([\d.]+)", line)
            if m:
                tm1 = float(m.group(1))
        elif "TM-score=" in line and "normalized by length of Structure_2" in line:
            m = re.search(r"TM-score=\s*([\d.]+)", line)
            if m:
                tm2 = float(m.group(1))
        elif line.strip().startswith("RMSD="):
            m = re.search(r"RMSD=\s*([\d.]+)", line)
            if m:
                rmsd = float(m.group(1))
    return tm1, tm2, rmsd


def compute_reference_centroids(df_biolip):
    """Cα centroid of the BioLiP binding-site residues, on each
    crystallized reference structure."""
    centroids = {}
    for og7, info in GRUPOS.items():
        for species, esp_info in info["especies"].items():
            if esp_info["modo"] != "cristal":
                continue
            row = df_biolip[
                (df_biolip["pdb_id"] == esp_info["pdb_id"]) & (df_biolip["chain"] == esp_info["chain"])
            ]
            if row.empty:
                print(f"  WARNING: no BioLiP row for {og7}|{species} ({esp_info['pdb_id']}_{esp_info['chain']})")
                continue
            residues_pdb = parse_pdb_residues(row.iloc[0]["binding_residues_pdb"])
            pdb_path = download_pdb(esp_info["pdb_id"], CRYSTAL_DIR)
            if pdb_path is None:
                continue
            coords = extract_ca_coords(pdb_path, esp_info["chain"], {n for _, n in residues_pdb})
            if not coords:
                print(f"  WARNING: no Cα coordinates found for {og7}|{species}")
                continue
            centroid = np.mean(list(coords.values()), axis=0)
            esp_info["pdb_path"] = pdb_path
            esp_info["centroide"] = centroid
            centroids[(og7, species)] = centroid
            print(f"  {og7} | {species} ({esp_info['pdb_id']}): centroid = {centroid}")
    return centroids


def find_best_anchor_and_transfer(og7, species, esp_info, group_info):
    """Downloads the AlphaFold model, aligns it against every crystallized
    species in the group, and transfers the centroid using whichever
    crystallized species gave the highest TM-score."""
    af_path = download_alphafold_pdb(esp_info["uniprot_id"], species, ALPHAFOLD_DIR)
    if af_path is None:
        print(f"  {og7} | {species}: could not download AlphaFold model")
        return
    esp_info["alphafold_pdb_path"] = af_path

    crystal_species = [(sp, info) for sp, info in group_info["especies"].items() if info["modo"] == "cristal"]
    best_anchor, best_tm, best_centroid = None, -1, None

    for anchor_species, anchor_info in crystal_species:
        matrix_path = USALIGN_DIR / f"{og7}_{species}_ANCLA_{anchor_species}_matrix.txt"
        stdout = run_usalign(af_path, "A", anchor_info["pdb_path"], anchor_info["chain"], matrix_path)
        tm1, tm2, _ = parse_tm_scores_and_rmsd(stdout)
        tm_values = [t for t in (tm1, tm2) if t is not None]
        tm_avg = np.mean(tm_values) if tm_values else None
        if tm_avg is None:
            continue
        t, U = parse_matrix_file(matrix_path)
        if t is None:
            continue
        # U/t map target -> mobile in US-align's convention; centroid is
        # expressed in the anchor's frame, so the inverse transform applies
        transferred = U.T @ (anchor_info["centroide"] - t)

        if tm_avg > best_tm:
            best_anchor, best_tm, best_centroid = anchor_species, tm_avg, transferred

    if best_anchor is None:
        print(f"  {og7} | {species}: could not align against any crystallized species")
        return

    esp_info["centroide"] = best_centroid
    esp_info["ancla_usada"] = best_anchor
    esp_info["tm_score_ancla"] = round(best_tm, 3)
    print(f"  {og7} | {species}: best anchor = {best_anchor} (TM={best_tm:.3f}) | centroid = {best_centroid}")


def main():
    ALPHAFOLD_DIR.mkdir(parents=True, exist_ok=True)
    CRYSTAL_DIR.mkdir(parents=True, exist_ok=True)
    USALIGN_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading BioLiP...")
    with gzip.open(BASE_DIR / "analisis_biolip" / "BioLiP.txt.gz", "rt", encoding="utf-8", errors="replace") as f:
        df_biolip = pd.read_csv(f, sep="\t", header=None, names=COLS_BIOLIP, low_memory=False)
    df_biolip["pdb_id"] = df_biolip["pdb_id"].astype(str).str.strip().str.lower()

    print("\nComputing reference centroids on crystallized structures...")
    compute_reference_centroids(df_biolip)

    print("\nTransferring centroids to AlphaFold species (dynamic anchor)...")
    for og7, info in GRUPOS.items():
        for species, esp_info in info["especies"].items():
            if esp_info["modo"] == "alphafold":
                find_best_anchor_and_transfer(og7, species, esp_info, info)

    GRUPOS_PICKLE.parent.mkdir(parents=True, exist_ok=True)
    with open(GRUPOS_PICKLE, "wb") as f:
        pickle.dump(GRUPOS, f)
    print(f"\nCheckpoint saved: {GRUPOS_PICKLE}")


if __name__ == "__main__":
    main()
