#!/usr/bin/env python3
"""
Two-part pLDDT evaluation for the pilot OG7 groups.

Part 1: for each species and each of the 2 OG7 groups, resolves its
UniProt ID and computes the whole-protein pLDDT statistics from the
AlphaFold confidence JSON.

Part 2: pocket-specific pLDDT -- using the AlphaFold PDB models already
downloaded (B-factor column = pLDDT) and the centroids transferred by
US-align, averages the pLDDT of residues within 8 Å of the transferred
centroid.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

# ── Part 1 config ──────────────────────────────────────────────────────
SEQS_DIR      = "/home/ssneider/disco-big/ssneider-env/TDR_Targets7.1/actualizados_short/sequences"
MAPPER_DIR    = "/home/ssneider/disco-big/ssneider-env/TDR_Targets7.1/gene_mapper/Ortho_vs_Uniprot"
NOT_ORTHO_DIR = "/big/lab/mercedesdg/TDR_Targets_7/OrthoMCL/genomes_v7/not_on_orthomcl/new_genomes"
ALPHAFOLD_DIR = "/scratch/lab/ssneider/TDR_Targets7.1/alphafold"

SPECIES_UNIPROT_DIRECT  = {"atha", "dmel", "osat", "cele", "ddis", "ecol", "mtub"}
SPECIES_NOT_ON_ORTHOMCL = {"ovo", "egr", "kpm", "loa", "sao"}
ORTHO_FILE_NAMES = {"egr": "egr_OrthoMCL_asignation_evaluation.tsv"}

OG7_OF_INTEREST = ["OG7_0006581", "OG7_0006003"]  # decanoic acid, arabinofuranose-P
SPECIES_OF_INTEREST = ["hsap", "dmel", "cele", "scer", "ecol", "kpm", "tgon", "pfal",
                        "egr", "loa", "ovo", "bmaa", "mtub"]

# ── Part 2 config ──────────────────────────────────────────────────────
BASE_DIR = "/big/lab/ssneider/ssneider-env/TDR_Targets7.1/PocketVec/pLDDT_candidatos"
AF_PDB_DIR = f"{BASE_DIR}/alphafold_pdbs"
POCKET_RADIUS_ANGSTROM = 8.0
REFERENCE_PDBS = {
    "OG7_0006581": f"{BASE_DIR}/pdb_referencias/1w66.pdb",
    "OG7_0006003": f"{BASE_DIR}/pdb_referencias/1o8b.pdb",
}


# ── Part 1: UniProt ID lookup + whole-protein pLDDT ────────────────────
def find_in_fasta_direct(species, target_og7):
    for pat in [f"{species}_aa_seqs_OrthoMCL-7.fasta", f"{species}_protein.fasta"]:
        fp = Path(SEQS_DIR) / pat
        if fp.exists():
            with open(fp) as f:
                for line in f:
                    if line.startswith(">"):
                        parts = line.strip().lstrip(">").split()
                        if len(parts) >= 2:
                            uid = parts[0].split("|")[1] if "|" in parts[0] else parts[0]
                            og = next((p for p in parts if p.startswith("OG7_")), None)
                            if og == target_og7:
                                return uid
            return None
    return None


def find_with_mapper(species, target_og7):
    ed = "tcru" if species == "tcru" else species
    mf = Path(MAPPER_DIR) / ed / "mapped_clean.csv"
    if not mf.exists():
        mf = Path(MAPPER_DIR) / "tcr" / "mapped_clean.csv"
    if not mf.exists():
        return None
    dm = pd.read_csv(mf, dtype=str)
    n2u = dict(zip(dm["OrthoMCL_ID"].str.strip(), dm["Uniprot_ID"].str.strip()))
    for pat in [f"{species}_aa_seqs_OrthoMCL-7.fasta", f"{species}_protein.fasta"]:
        fp = Path(SEQS_DIR) / pat
        if fp.exists():
            with open(fp) as f:
                for line in f:
                    if line.startswith(">"):
                        parts = line.strip().lstrip(">").split()
                        if len(parts) >= 2:
                            nid = parts[0].split("|")[1] if "|" in parts[0] else parts[0]
                            og = next((p for p in parts if p.startswith("OG7_")), None)
                            if og == target_og7:
                                uid = n2u.get(nid)
                                if uid:
                                    return uid
            return None
    return None


def find_not_on_orthomcl(species, target_og7):
    nb = ORTHO_FILE_NAMES.get(species, f"{species}_orthomcl7.1.txt")
    of = Path(NOT_ORTHO_DIR) / species / nb
    mf = Path(MAPPER_DIR) / species / "mapped_clean.csv"
    if not of.exists() or not mf.exists():
        return None
    dm = pd.read_csv(mf, dtype=str)
    n2u = dict(zip(dm["OrthoMCL_ID"].str.strip(), dm["Uniprot_ID"].str.strip()))
    do = pd.read_csv(of, sep="\t", dtype=str)
    do.columns = do.columns.str.strip()
    for _, row in do.iterrows():
        qid = str(row[do.columns[0]]).strip()
        og = str(row[do.columns[2]]).strip()
        if og != target_og7:
            continue
        nid = qid.split("|")[1] if "|" in qid else qid
        uid = n2u.get(nid)
        if uid:
            return uid
    return None


def find_alphafold_json(species, uniprot_id):
    if not uniprot_id:
        return None
    folder = Path(ALPHAFOLD_DIR) / species
    if not folder.exists():
        return None
    candidates = list(folder.glob(f"AF-{uniprot_id}-F1-confidence_v4.json"))
    return candidates[0] if candidates else None


def compute_plddt_stats(json_path):
    with open(json_path) as f:
        data = json.load(f)
    scores = data["confidenceScore"]
    n = len(scores)
    return {
        "n_residuos": n,
        "plddt_promedio": round(sum(scores) / n, 2),
        "plddt_min": min(scores),
        "plddt_max": max(scores),
        "pct_residuos_plddt>=70": round(sum(1 for s in scores if s >= 70) / n * 100, 1),
        "pct_residuos_plddt>=90": round(sum(1 for s in scores if s >= 90) / n * 100, 1),
    }


def part1_whole_protein_plddt():
    results = []
    print("Looking up UniProt IDs per species for each OG7 of interest...")
    for og7 in OG7_OF_INTEREST:
        for species in SPECIES_OF_INTEREST:
            if species in SPECIES_UNIPROT_DIRECT:
                uid = find_in_fasta_direct(species, og7)
            elif species in SPECIES_NOT_ON_ORTHOMCL:
                uid = find_not_on_orthomcl(species, og7)
            else:
                uid = find_with_mapper(species, og7)
            results.append({"og7": og7, "especie": species, "uniprot_id": uid})
            print(f"  {og7} | {species}: {uid if uid else 'NOT FOUND'}")

    df_uniprot = pd.DataFrame(results)
    df_uniprot.to_csv("og7_uniprot_por_especie.tsv", sep="\t", index=False)
    print("\nSaved: og7_uniprot_por_especie.tsv")

    print("\nComputing pLDDT for each species/OG7 pair...")
    plddt_results = []
    for _, row in df_uniprot.iterrows():
        json_path = find_alphafold_json(row["especie"], row["uniprot_id"])
        if json_path is None:
            plddt_results.append({
                "og7": row["og7"], "especie": row["especie"], "uniprot_id": row["uniprot_id"],
                "json_encontrado": False, "n_residuos": None, "plddt_promedio": None,
                "plddt_min": None, "plddt_max": None,
                "pct_residuos_plddt>=70": None, "pct_residuos_plddt>=90": None,
            })
            continue
        stats = compute_plddt_stats(json_path)
        plddt_results.append({"og7": row["og7"], "especie": row["especie"],
                               "uniprot_id": row["uniprot_id"], "json_encontrado": True, **stats})

    df_plddt = pd.DataFrame(plddt_results)
    df_plddt.to_csv("og7_plddt_por_especie.tsv", sep="\t", index=False)
    print("\n── pLDDT RESULTS BY SPECIES AND OG7 ────────────────────────────────────")
    print(df_plddt.to_string(index=False))

    pivot = df_plddt.pivot(index="especie", columns="og7", values="plddt_promedio")
    print("\nSpecies with data available for BOTH OG7 groups:")
    print(pivot.dropna().to_string())


# ── Part 2: pocket-region pLDDT ─────────────────────────────────────────
def read_ca_bfactor(pdb_path):
    """Reads an AlphaFold PDB and returns, per residue: (res_num, x, y, z,
    bfactor). In AlphaFold, B-factor = pLDDT."""
    residues = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            res_num = int(line[22:26].strip())
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            bfactor = float(line[60:66].strip())
            residues.append((res_num, x, y, z, bfactor))
    return residues


def plddt_pocket_region(pdb_path, centroid, radius=POCKET_RADIUS_ANGSTROM):
    residues = read_ca_bfactor(pdb_path)
    if not residues:
        return None, 0, []
    center = np.array(centroid)
    pocket_residues = []
    for res_num, x, y, z, bfactor in residues:
        dist = np.linalg.norm(np.array([x, y, z]) - center)
        if dist <= radius:
            pocket_residues.append((res_num, round(dist, 2), bfactor))
    if not pocket_residues:
        return None, 0, []
    plddt_mean = np.mean([r[2] for r in pocket_residues])
    return round(plddt_mean, 2), len(pocket_residues), pocket_residues


def part2_pocket_plddt():
    df_alignment = pd.read_csv(f"{BASE_DIR}/centroides_transferidos.tsv", sep="\t")
    print(f"Entries loaded: {len(df_alignment)}")

    results = []
    for _, row in df_alignment.iterrows():
        og7, species, uid = row["og7"], row["especie"], row["uniprot_id"]

        if row["es_referencia"]:
            # No AlphaFold pLDDT applies to the crystallized reference itself
            results.append({
                "og7": og7, "especie": species, "uniprot_id": uid, "tm_score": row["tm_score"],
                "plddt_pocket": "N/A (crystal)", "n_residuos_pocket": None,
                "centroide_x": row["centroide_transferido_x"],
                "centroide_y": row["centroide_transferido_y"],
                "centroide_z": row["centroide_transferido_z"],
                "residuos_detalle": None,
            })
            print(f"  {og7} | {species}: crystallized reference structure (pLDDT not applicable)")
            continue

        if pd.isna(uid) or pd.isna(row["centroide_transferido_x"]):
            print(f"  {og7} | {species}: no centroid or no UniProt ID, skipping")
            results.append({
                "og7": og7, "especie": species, "uniprot_id": uid, "tm_score": row["tm_score"],
                "plddt_pocket": None, "n_residuos_pocket": None,
                "centroide_x": None, "centroide_y": None, "centroide_z": None, "residuos_detalle": None,
            })
            continue

        pdb_path = Path(AF_PDB_DIR) / f"AF-{uid}-F1_{species}.pdb"
        if not pdb_path.exists():
            print(f"  {og7} | {species}: PDB not found at {pdb_path}")
            results.append({
                "og7": og7, "especie": species, "uniprot_id": uid, "tm_score": row["tm_score"],
                "plddt_pocket": None, "n_residuos_pocket": None,
                "centroide_x": None, "centroide_y": None, "centroide_z": None, "residuos_detalle": None,
            })
            continue

        centroid = [row["centroide_transferido_x"], row["centroide_transferido_y"], row["centroide_transferido_z"]]
        plddt, n_res, detail = plddt_pocket_region(pdb_path, centroid)

        status = f"pLDDT_pocket={plddt} ({n_res} residues within {POCKET_RADIUS_ANGSTROM}Å)"
        if plddt is None:
            status = f"ALERT: no residue within {POCKET_RADIUS_ANGSTROM}Å of the centroid"
        print(f"  {og7} | {species}: {status}")

        results.append({
            "og7": og7, "especie": species, "uniprot_id": uid, "tm_score": row["tm_score"],
            "plddt_pocket": plddt, "n_residuos_pocket": n_res,
            "centroide_x": row["centroide_transferido_x"],
            "centroide_y": row["centroide_transferido_y"],
            "centroide_z": row["centroide_transferido_z"],
            "residuos_detalle": str(detail),
        })

    df_res = pd.DataFrame(results)
    df_res.to_csv(f"{BASE_DIR}/plddt_pocket_region.tsv", sep="\t", index=False)

    print("\n── POCKET-REGION pLDDT SUMMARY ─────────────────────────────────────────")
    cols = ["og7", "especie", "tm_score", "plddt_pocket", "n_residuos_pocket"]
    print(df_res[cols].to_string(index=False))

    print("\n── ALERTS ───────────────────────────────────────────────────────────")
    alerts = df_res[
        (df_res["plddt_pocket"].apply(lambda x: isinstance(x, float) and x < 70)) |
        (df_res["n_residuos_pocket"] == 0)
    ]
    print("None -- every pocket has an acceptable pLDDT" if len(alerts) == 0
          else alerts[["og7", "especie", "plddt_pocket", "n_residuos_pocket"]].to_string(index=False))


if __name__ == "__main__":
    part1_whole_protein_plddt()
    part2_pocket_plddt()
