#!/usr/bin/env python3
"""
Sequence-divergence scoring of the 247 candidates, to pick the 3 finalist
OG7 groups for the second PocketVec round.

For each OG7 group and each crystallized ligand, the pocket residues
(BioLiP binding-site positions) are mapped onto a multiple sequence
alignment (MAFFT) built with one "anchor" sequence per distinct
crystallized structure plus one sequence per species in the group. The
pocket columns are then compared pairwise between every pair of
taxonomic lineages represented in the group, and that divergence is
scored against a null distribution built by resampling the same number
of random alignment columns 100 times -- giving a z-score and an
empirical p-value per (ligand, lineage pair).

A group's final score is its fraction of evaluated ligands with
z-score >= 2 ("frac_significativos"). The 3 groups with the highest
fraction were selected as finalists, discarding groups supported by
fewer than 4-5 evaluated ligands.
"""
import random
import subprocess
from collections import Counter
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

IN_DIR    = "biolip_exploracion"
OUT_DIR   = "biolip_exploracion2"
NUEVA_DIR = Path("Nueva_estrategia_pocketvec")
WORK_DIR  = NUEVA_DIR / "msa_temp"
MAFFT_BIN = "/home/ssneider/miniconda3/bin/mafft"

SEQS_DIR      = "/home/ssneider/disco-big/ssneider-env/TDR_Targets7.1/actualizados_short/sequences"
MAPPER_DIR    = "/home/ssneider/disco-big/ssneider-env/TDR_Targets7.1/gene_mapper/Ortho_vs_Uniprot"
NOT_ORTHO_DIR = "/big/lab/mercedesdg/TDR_Targets_7/OrthoMCL/genomes_v7/not_on_orthomcl/new_genomes"
SPECIES_UNIPROT_DIRECT  = {"atha", "dmel", "osat", "cele", "ddis", "ecol", "mtub"}
SPECIES_NOT_ON_ORTHOMCL = {"ovo", "egr", "kpm", "loa", "sao"}
ORTHO_FILE_NAMES = {"egr": "egr_OrthoMCL_asignation_evaluation.tsv"}

LINEAGE = {
    "hsap": "vertebrate", "mmus": "vertebrate",
    "dmel": "model_invertebrate", "cele": "model_invertebrate",
    "scer": "fungus", "calb": "fungus",
    "atha": "plant", "osat": "plant",
    "ddis": "model_amoeba",
    "ecol": "bacteria", "mtub": "bacteria", "kpm": "bacteria", "sao": "bacteria",
    "ldon": "kinetoplastid", "lmaj": "kinetoplastid",
    "tbrt": "kinetoplastid", "tcru": "kinetoplastid", "tcrc-l-n": "kinetoplastid",
    "pfal": "apicomplexan", "pvip": "apicomplexan", "tgon": "apicomplexan",
    "bmaa": "nematode", "loa": "nematode", "ovo": "nematode",
    "egr": "cestode",
    "gass": "other_parasite", "gmur": "other_parasite",
    "ehia": "other_parasite", "tvag": "other_parasite",
}

COLS_BIOLIP = [
    "pdb_id", "chain", "resolution", "binding_site",
    "ligand_id", "ligand_chain", "ligand_serial",
    "binding_residues_pdb", "binding_residues_renum",
    "catalytic_pdb", "catalytic_renum",
    "ec_number", "go_terms",
    "affinity_manual", "affinity_moad", "affinity_pdbbind", "affinity_bindingdb",
    "uniprot_id", "pubmed_id", "ligand_seqnum", "receptor_seq",
]
LIGANDS_EXCLUDE = {
    "HOH", "DOD", "WAT", "SO4", "PO4", "GOL", "EDO", "PEG", "BME",
    "ACT", "ACE", "FMT", "MPD", "MES", "TAR", "CIT", "TRS", "EPE",
    "MG", "ZN", "CA", "FE", "MN", "NA", "K", "CL", "CU", "NI", "CO",
    "SE", "CD", "HG", "PT", "AU", "AG", "IOD", "BR", "ION", "FLC", "DMS",
}

FINALIST_MIN_LIGANDS = 4  # discard groups supported by fewer evaluated ligands


def load_biolip_with_og7(og7_set):
    import gzip
    with gzip.open("/home/ssneider/disco-big/ssneider-env/TDR_Targets7.1/PocketVec/analisis_biolip/BioLiP.txt.gz",
                    "rt", encoding="utf-8", errors="replace") as f:
        df_bio = pd.read_csv(f, sep="\t", header=None, names=COLS_BIOLIP, low_memory=False)
    df_bio["resolution"] = pd.to_numeric(df_bio["resolution"], errors="coerce")
    df_bio["uniprot_first"] = df_bio["uniprot_id"].astype(str).str.split(",").str[0].str.strip()
    df_bio["has_uniprot"] = df_bio["uniprot_first"].str.match(r"^[A-Z][0-9][A-Z0-9]{3}[0-9]$")
    df_bio["ligand_id"] = df_bio["ligand_id"].astype(str).str.strip().str.upper()
    df_bio["ligand_bio"] = ~df_bio["ligand_id"].isin(LIGANDS_EXCLUDE)
    df_bio_clean = df_bio[(df_bio["resolution"] > 0) & df_bio["has_uniprot"] & df_bio["ligand_bio"]].copy()

    uniprot_to_og = {}
    mapper_path = Path(MAPPER_DIR)
    all_mapper_species = {p.name for p in mapper_path.iterdir() if p.is_dir()}
    species_mapper = all_mapper_species - SPECIES_UNIPROT_DIRECT - SPECIES_NOT_ON_ORTHOMCL
    if "tcr" in species_mapper:
        species_mapper.discard("tcr")
        species_mapper.add("tcru")

    for species in SPECIES_UNIPROT_DIRECT:
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
                                if og and og in og7_set:
                                    uniprot_to_og[uid] = og
                break
    for species in sorted(species_mapper):
        ed = "tcru" if species == "tcru" else species
        mf = mapper_path / ed / "mapped_clean.csv"
        if not mf.exists():
            mf = mapper_path / "tcr" / "mapped_clean.csv"
        if not mf.exists():
            continue
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
                                uid = n2u.get(nid)
                                if og and uid and og in og7_set:
                                    uniprot_to_og[uid] = og
                break
    for species in sorted(SPECIES_NOT_ON_ORTHOMCL):
        nb = ORTHO_FILE_NAMES.get(species, f"{species}_orthomcl7.1.txt")
        of = Path(NOT_ORTHO_DIR) / species / nb
        mf = mapper_path / species / "mapped_clean.csv"
        if not of.exists() or not mf.exists():
            continue
        dm = pd.read_csv(mf, dtype=str)
        n2u = dict(zip(dm["OrthoMCL_ID"].str.strip(), dm["Uniprot_ID"].str.strip()))
        do = pd.read_csv(of, sep="\t", dtype=str)
        do.columns = do.columns.str.strip()
        for _, row in do.iterrows():
            qid = str(row[do.columns[0]]).strip()
            og = str(row[do.columns[2]]).strip()
            if not og.startswith("OG7_") or og not in og7_set:
                continue
            nid = qid.split("|")[1] if "|" in qid else qid
            uid = n2u.get(nid)
            if uid:
                uniprot_to_og[uid] = og

    df_bio_clean["og7"] = df_bio_clean["uniprot_first"].map(uniprot_to_og)
    return df_bio_clean[df_bio_clean["og7"].notna()].copy()


def parse_pocket_residues(binding_str):
    """binding_residues_renum is already 1-based from BioLiP -> 0-based positions."""
    if not isinstance(binding_str, str) or not binding_str.strip():
        return []
    import re
    return [int(m.group(1)) - 1 for tok in binding_str.split()
            if (m := re.match(r"^[A-Za-z]?(\d+)$", tok))]


def list_candidate_pockets(df_with_og, og7, lig_categoria):
    """One entry per distinct ligand crystallized for this OG7 (best resolution)."""
    grp = df_with_og[df_with_og["og7"] == og7]
    if len(grp) == 0:
        return []
    pockets = []
    for lig_id, sub in grp.groupby("ligand_id"):
        row = sub.sort_values("resolution").iloc[0]
        pocket_pos = parse_pocket_residues(row["binding_residues_renum"])
        if not pocket_pos:
            continue
        pockets.append({
            "pdb_id": row["pdb_id"], "ligand_id": lig_id, "receptor_seq": row["receptor_seq"],
            "pocket_pos": pocket_pos, "categoria_ligando": lig_categoria.get(lig_id, "no_confirmado"),
        })
    return pockets


def extract_og7_sequences(og7):
    """One representative sequence per species for this OG7 group."""
    sequences = {}

    def scan_fasta(species, fp, target_og7=None, target_qids=None):
        with open(fp) as f:
            capturing, seq_lines = False, []
            for line in f:
                if line.startswith(">"):
                    if capturing and seq_lines:
                        sequences[species] = "".join(seq_lines)
                    if target_og7 is not None:
                        parts = line.strip().lstrip(">").split()
                        og = next((p for p in parts if p.startswith("OG7_")), None)
                        capturing = (og == target_og7)
                    else:
                        qid = line.strip().lstrip(">").split()[0]
                        capturing = qid in target_qids
                    seq_lines = []
                elif capturing:
                    seq_lines.append(line.strip())
            if capturing and seq_lines:
                sequences[species] = "".join(seq_lines)

    for species in SPECIES_UNIPROT_DIRECT:
        for pat in [f"{species}_aa_seqs_OrthoMCL-7.fasta", f"{species}_protein.fasta"]:
            fp = Path(SEQS_DIR) / pat
            if fp.exists():
                scan_fasta(species, fp, target_og7=og7)
                break

    mapper_path = Path(MAPPER_DIR)
    all_mapper_species = {p.name for p in mapper_path.iterdir() if p.is_dir()}
    species_mapper = all_mapper_species - SPECIES_UNIPROT_DIRECT - SPECIES_NOT_ON_ORTHOMCL
    if "tcr" in species_mapper:
        species_mapper.discard("tcr")
        species_mapper.add("tcru")
    for species in species_mapper:
        for pat in [f"{species}_aa_seqs_OrthoMCL-7.fasta", f"{species}_protein.fasta"]:
            fp = Path(SEQS_DIR) / pat
            if fp.exists():
                scan_fasta(species, fp, target_og7=og7)
                break

    for species in SPECIES_NOT_ON_ORTHOMCL:
        nb = ORTHO_FILE_NAMES.get(species, f"{species}_orthomcl7.1.txt")
        of = Path(NOT_ORTHO_DIR) / species / nb
        if not of.exists():
            continue
        do = pd.read_csv(of, sep="\t", dtype=str)
        do.columns = do.columns.str.strip()
        target_qids = set(do[do[do.columns[2]].str.strip() == og7][do.columns[0]].str.strip())
        if not target_qids:
            continue
        for pat in [f"{species}_aa_seqs_OrthoMCL-7.fasta", f"{species}_protein.fasta"]:
            fp = Path(SEQS_DIR) / pat
            if fp.exists():
                scan_fasta(species, fp, target_qids=target_qids)
                break

    return sequences


def build_msa_multi_anchor(og7, pockets):
    sequences = extract_og7_sequences(og7)
    if not sequences:
        return None

    unique_anchors = {}
    for p in pockets:
        if p["receptor_seq"] not in unique_anchors:
            unique_anchors[p["receptor_seq"]] = p["pdb_id"]

    fasta_in = WORK_DIR / f"{og7}_input.fasta"
    with open(fasta_in, "w") as f:
        for seq, pdb_id in unique_anchors.items():
            f.write(f">ANCLA_{pdb_id}\n{seq}\n")
        for species, seq in sequences.items():
            f.write(f">{species}\n{seq}\n")

    result = subprocess.run([MAFFT_BIN, "--auto", "--quiet", str(fasta_in)], capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  MAFFT ERROR for {og7}: {result.stderr[:300]}")
        return None

    aligned, current, seq_lines = {}, None, []
    for line in result.stdout.splitlines():
        if line.startswith(">"):
            if current:
                aligned[current] = "".join(seq_lines)
            current, seq_lines = line[1:].strip(), []
        else:
            seq_lines.append(line.strip())
    if current:
        aligned[current] = "".join(seq_lines)

    seq_to_anchor_id = {seq: f"ANCLA_{pdb_id}" for seq, pdb_id in unique_anchors.items()}
    return aligned, seq_to_anchor_id


def map_pocket_to_msa_columns(aligned_anchor, pocket_pos_0based):
    columns, idx, pos_set = [], -1, set(pocket_pos_0based)
    for col, char in enumerate(aligned_anchor):
        if char != "-":
            idx += 1
            if idx in pos_set:
                columns.append(col)
    return columns


def pocket_distance(seq1, seq2, cols):
    diffs, comparable = 0, 0
    for c in cols:
        if c >= len(seq1) or c >= len(seq2):
            continue
        a, b = seq1[c], seq2[c]
        if a == "-" or b == "-":
            continue
        comparable += 1
        if a != b:
            diffs += 1
    return diffs / comparable if comparable > 0 else None


def usable_columns(aligned_seqs, max_gap_frac=0.3):
    species = [e for e in aligned_seqs if not e.startswith("ANCLA_")]
    length = len(next(iter(aligned_seqs.values())))
    return [c for c in range(length)
            if sum(1 for e in species if aligned_seqs[e][c] == "-") / len(species) <= max_gap_frac]


def mean_pairwise_distance(aligned_seqs, cols, species_a, species_b):
    inter = [d for a in species_a for b in species_b
             if (d := pocket_distance(aligned_seqs[a], aligned_seqs[b], cols)) is not None]
    intra = []
    for group in (species_a, species_b):
        for a, b in combinations(group, 2):
            d = pocket_distance(aligned_seqs[a], aligned_seqs[b], cols)
            if d is not None:
                intra.append(d)
    if not inter:
        return None
    return (sum(inter) / len(inter)) - (sum(intra) / len(intra) if intra else 0.0)


def score_vs_null(aligned_seqs, pocket_cols, species_a, species_b, valid_cols, n_perm=100, seed=0):
    rng = random.Random(seed)
    real = mean_pairwise_distance(aligned_seqs, pocket_cols, species_a, species_b)
    if real is None:
        return None
    n = len(pocket_cols)
    pool = [c for c in valid_cols if c not in pocket_cols] or valid_cols
    null_scores = []
    for _ in range(n_perm):
        sample = rng.sample(pool, n) if len(pool) >= n else rng.choices(pool, k=n)
        s = mean_pairwise_distance(aligned_seqs, sample, species_a, species_b)
        if s is not None:
            null_scores.append(s)
    if len(null_scores) < 20:
        return None
    mean_null, std_null = np.mean(null_scores), (np.std(null_scores) or 1e-6)
    z = (real - mean_null) / std_null
    p_emp = sum(1 for s in null_scores if s >= real) / len(null_scores)
    return {"score_real": real, "z_score": z, "p_empirico": p_emp}


def validated_divergence_score(aligned_seqs, pocket_cols, lineage_dict, valid_cols, n_perm=100):
    groups_all = {}
    for e in aligned_seqs:
        if e.startswith("ANCLA_"):
            continue
        groups_all.setdefault(lineage_dict.get(e, "other"), []).append(e)
    valid_groups = {l: sp for l, sp in groups_all.items() if len(sp) >= 2}

    best_pair, best_result = None, None
    for l1, l2 in combinations(valid_groups.keys(), 2):
        r = score_vs_null(aligned_seqs, pocket_cols, valid_groups[l1], valid_groups[l2], valid_cols, n_perm=n_perm)
        if r and (best_result is None or r["z_score"] > best_result["z_score"]):
            best_pair, best_result = (l1, l2), r
    return best_pair, best_result


def evaluate_og7(og7, df_with_og, lig_categoria, lineage_dict, min_residues=5, n_perm=100):
    pockets = list_candidate_pockets(df_with_og, og7, lig_categoria)
    if not pockets:
        return []
    msa_result = build_msa_multi_anchor(og7, pockets)
    if msa_result is None:
        return []
    aligned, seq_to_anchor_id = msa_result
    valid_cols = usable_columns(aligned)

    results = []
    for p in pockets:
        anchor_id = seq_to_anchor_id.get(p["receptor_seq"])
        if anchor_id is None or anchor_id not in aligned:
            continue
        pocket_cols = map_pocket_to_msa_columns(aligned[anchor_id], p["pocket_pos"])
        if len(pocket_cols) < min_residues:
            continue
        pair, result = validated_divergence_score(aligned, pocket_cols, lineage_dict, valid_cols, n_perm=n_perm)
        if result is None:
            continue
        results.append({
            "og7": og7, "ligand_id": p["ligand_id"], "categoria_ligando": p["categoria_ligando"],
            "n_residuos_pocket": len(pocket_cols), "par_linajes": pair,
            "z_score": round(result["z_score"], 2), "p_empirico": round(result["p_empirico"], 3),
        })
    return results


def summarize_og7(og7, detail_rows, min_evaluated=FINALIST_MIN_LIGANDS):
    if len(detail_rows) < min_evaluated:
        return None
    z_scores = [r["z_score"] for r in detail_rows]
    n_significant = sum(1 for z in z_scores if z >= 2)
    return {
        "og7": og7,
        "n_ligandos_evaluados": len(detail_rows),
        "n_significativos": n_significant,
        "frac_significativos": round(n_significant / len(detail_rows), 3),
        "z_score_max": round(max(z_scores), 2),
    }


def main():
    NUEVA_DIR.mkdir(exist_ok=True)
    WORK_DIR.mkdir(exist_ok=True)

    candidates = pd.read_csv(f"{NUEVA_DIR}/candidatos_divergencia_stage1.tsv", sep="\t")
    lig_categoria = {}  # rebuilt from 01_candidate_filter_247.py if needed standalone
    df_with_og = load_biolip_with_og7(set(candidates["og7"]))

    print(f"Evaluating sequence divergence for {len(candidates)} candidates...")
    summaries, all_details = [], []
    for i, og7 in enumerate(candidates["og7"]):
        if i % 10 == 0:
            print(f"  {i}/{len(candidates)} OG7 processed...")
        details = evaluate_og7(og7, df_with_og, lig_categoria, LINEAGE)
        all_details.extend(details)
        summary = summarize_og7(og7, details)
        if summary is not None:
            summaries.append(summary)

    df_summary = pd.DataFrame(summaries).sort_values(
        ["frac_significativos", "z_score_max"], ascending=False
    ).reset_index(drop=True)
    df_detail = pd.DataFrame(all_details)

    df_summary.to_csv(f"{NUEVA_DIR}/resumen_divergencia_og7_v2.tsv", sep="\t", index=False)
    df_detail.to_csv(f"{NUEVA_DIR}/detalle_divergencia_og7_v2.tsv", sep="\t", index=False)

    print(f"\nOG7 with a divergence summary: {len(df_summary)}")
    print(df_summary.head(20).to_string(index=False))
    print("\nTop 3 finalists (highest fraction of ligands with z-score >= 2):")
    print(df_summary.head(3).to_string(index=False))


if __name__ == "__main__":
    main()
