#!/usr/bin/env python3
"""
Enriches the corrected candidate list (992 OG7 groups) with ligand names,
per-group ligand diversity, and a handful of overview figures used to
sanity-check the candidate pool before the unique-ligand selection step.
"""
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

IN_DIR = "biolip_exploracion"
OUT_DIR = "biolip_exploracion2"
LIGAND_FILE = "ligand.tsv.gz"

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
COMMON_COFACTORS = {
    "ATP", "ADP", "AMP", "GTP", "GDP", "GMP", "NAD", "NAP", "NDP", "FAD", "FMN",
    "SAM", "SAH", "COA", "HEM", "HEC", "FE2", "ZN", "MG", "CA", "CL", "PLP", "PMP",
    "PYR", "OXA", "CIT", "MAL", "FUM", "SUC", "AKG", "OXL", "PGA", "PEP", "FBP",
    "UMP", "CMP", "TMP", "DUT", "DUP", "PPK", "IMP", "ITT", "APC", "2PG", "G3P",
    "GSH", "DND", "FFO", "PLS", "ANP", "AD9", "HSX", "B3M", "ME8", "MET",
    "LYS", "KAA", "PRO", "RNA", "DNA", "PEPTIDE",
}

SEQS_DIR      = "/home/ssneider/disco-big/ssneider-env/TDR_Targets7.1/actualizados_short/sequences"
MAPPER_DIR    = "/home/ssneider/disco-big/ssneider-env/TDR_Targets7.1/gene_mapper/Ortho_vs_Uniprot"
NOT_ORTHO_DIR = "/big/lab/mercedesdg/TDR_Targets_7/OrthoMCL/genomes_v7/not_on_orthomcl/new_genomes"
SPECIES_UNIPROT_DIRECT  = {"atha", "dmel", "osat", "cele", "ddis", "ecol", "mtub"}
SPECIES_NOT_ON_ORTHOMCL = {"ovo", "egr", "kpm", "loa", "sao"}
ORTHO_FILE_NAMES = {"egr": "egr_OrthoMCL_asignation_evaluation.tsv"}


def load_ligand_names():
    lig_names, lig_formula = {}, {}
    if Path(LIGAND_FILE).exists():
        cols = ["lig_id", "formula", "inchi", "inchikey", "smiles", "name", "chebi", "drugbank", "zinc"]
        with gzip.open(LIGAND_FILE, "rt", encoding="utf-8", errors="replace") as f:
            df_lig = pd.read_csv(f, sep="\t", header=None, names=cols,
                                  on_bad_lines="warn", engine="python")
        df_lig["lig_id"] = df_lig["lig_id"].astype(str).str.strip().str.upper()
        df_lig["short_name"] = df_lig["name"].astype(str).str.split(";").str[0].str.strip()
        lig_names = dict(zip(df_lig["lig_id"], df_lig["short_name"]))
        lig_formula = dict(zip(df_lig["lig_id"], df_lig["formula"].astype(str)))
        print(f"Named ligands loaded: {len(lig_names):,}")
    return lig_names, lig_formula


def build_uniprot_to_og7():
    """Same three-path ID resolution used across the project (Gene Mapper +
    direct UniProt species + not-on-OrthoMCL species)."""
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
                                if og:
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
                                if og and uid:
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
            if not og.startswith("OG7_"):
                continue
            nid = qid.split("|")[1] if "|" in qid else qid
            uid = n2u.get(nid)
            if uid:
                uniprot_to_og[uid] = og

    return uniprot_to_og


def analyze_og7_ligands(df_og7):
    ligs = df_og7["ligand_id"].value_counts()
    n_unique = len(ligs)
    top3 = list(ligs.head(3).index)
    # Rare ligand: not a common cofactor, present in >= 2 structures
    rare = [l for l in ligs.index if l not in COMMON_COFACTORS and ligs[l] >= 2]
    rare_ligand = rare[0] if rare else None
    return n_unique, top3, rare_ligand


def main():
    lig_names, lig_formula = load_ligand_names()

    df_cand = pd.read_csv(f"{IN_DIR}/candidatos_pocketvec_v2.tsv", sep="\t")
    df_ligs = pd.read_csv(f"{IN_DIR}/distribucion_ligandos.tsv", sep="\t")
    df_ligs["ligand_name"] = df_ligs["ligand_id"].map(lig_names).fillna(df_ligs["ligand_id"])

    print("Loading BioLiP...")
    with gzip.open("BioLiP.txt.gz", "rt", encoding="utf-8", errors="replace") as f:
        df_bio = pd.read_csv(f, sep="\t", header=None, names=COLS_BIOLIP, low_memory=False)

    df_bio["resolution"] = pd.to_numeric(df_bio["resolution"], errors="coerce")
    df_bio["uniprot_first"] = df_bio["uniprot_id"].astype(str).str.split(",").str[0].str.strip()
    df_bio["has_uniprot"] = df_bio["uniprot_first"].str.match(r"^[A-Z][0-9][A-Z0-9]{3}[0-9]$")
    df_bio["ligand_id"] = df_bio["ligand_id"].astype(str).str.strip().str.upper()
    df_bio["ligand_bio"] = ~df_bio["ligand_id"].isin(LIGANDS_EXCLUDE)
    df_bio_clean = df_bio[(df_bio["resolution"] > 0) & df_bio["has_uniprot"] & df_bio["ligand_bio"]].copy()

    uniprot_to_og = build_uniprot_to_og7()
    df_bio_clean["og7"] = df_bio_clean["uniprot_first"].map(uniprot_to_og)
    df_with_og = df_bio_clean[df_bio_clean["og7"].notna()].copy()
    print(f"BioLiP entries cross-referenced with OG7: {len(df_with_og):,}")

    df_cand["prot_por_especie"] = (df_cand["n_prot_red"] / df_cand["n_especies_red"]).round(1)

    n_unique_map, rare_map = {}, {}
    for og7, grp in df_with_og.groupby("og7"):
        nu, _, rare = analyze_og7_ligands(grp)
        n_unique_map[og7] = nu
        rare_map[og7] = rare

    df_cand["n_ligandos_unicos"] = df_cand["og7"].map(n_unique_map)
    df_cand["ligando_raro"] = df_cand["og7"].map(rare_map)
    df_cand["nombre_ligando_raro"] = df_cand["ligando_raro"].map(lig_names)

    df_cand_sorted = df_cand.sort_values(
        ["n_especies_red", "n_uniprot_biolip"], ascending=False
    ).reset_index(drop=True)

    Path(OUT_DIR).mkdir(exist_ok=True)
    df_cand_sorted.to_csv(f"{OUT_DIR}/candidatos_final.tsv", sep="\t", index=False)

    cols_print = ["og7", "n_especies_red", "n_prot_red", "prot_por_especie",
                  "n_especies_biolip", "n_uniprot_biolip", "resolucion_media",
                  "n_ligandos_unicos", "ligando_raro", "nombre_ligando_raro", "ec_numbers"]
    print("\n── TOP 20 FINAL CANDIDATES ──────────────────────────────────────────────")
    print(df_cand_sorted[cols_print].head(20).to_string(index=False))

    # ── Ligands with >100 structures: who are they ──────────────────────────
    df_big_ligs = df_ligs[df_ligs["n_entradas"] > 100].copy()
    df_big_ligs["ligand_name"] = df_big_ligs["ligand_id"].map(lig_names).fillna("?")
    df_big_ligs["formula"] = df_big_ligs["ligand_id"].map(lig_formula).fillna("?")
    print(f"\nTotal ligands with >100 structures: {len(df_big_ligs)}")
    df_big_ligs.to_csv(f"{OUT_DIR}/ligandos_mas100.tsv", sep="\t", index=False)

    # ── Plot 1: log-log distribution, annotated with the top 15 candidates ──
    top15 = df_cand_sorted.head(15)
    og7_lig_freq = {}
    for _, row in top15.iterrows():
        if pd.isna(row["top_ligandos"]):
            continue
        code = str(row["top_ligandos"]).split(",")[0].strip()
        freq_row = df_ligs[df_ligs["ligand_id"] == code]
        if not freq_row.empty:
            og7_lig_freq[row["og7"]] = (code, int(freq_row.iloc[0]["n_entradas"]),
                                         lig_names.get(code, code))

    freq_dist = df_ligs["n_entradas"].value_counts().sort_index()

    fig, ax = plt.subplots(figsize=(14, 8))
    ax.loglog(freq_dist.index, freq_dist.values, "o", color="#AACCE8",
              markersize=3, alpha=0.6, label="All ligands")
    for _, bl_row in df_big_ligs.iterrows():
        freq_val = bl_row["n_entradas"]
        n_ligs_with_freq = freq_dist.get(freq_val, 1)
        ax.loglog(freq_val, n_ligs_with_freq, "o", color="#E84040", markersize=5, alpha=0.8)
    ax.axvspan(100, df_ligs["n_entradas"].max() * 1.5, alpha=0.06, color="red",
               label="Zone >100 structures")
    ax.axvline(100, color="red", linestyle="--", linewidth=1, alpha=0.5)

    colors_top15 = plt.cm.tab20(np.linspace(0, 1, 15))
    for i, (_, row) in enumerate(top15.iterrows()):
        og7 = row["og7"]
        if og7 not in og7_lig_freq:
            continue
        code, freq, name = og7_lig_freq[og7]
        n_ligs_freq = freq_dist.get(freq, 1)
        x_pt, y_pt = freq, n_ligs_freq
        x_txt = x_pt * (1.8 if i % 2 == 0 else 2.5)
        y_txt = y_pt * (3 if i % 3 == 0 else 0.4 if i % 3 == 1 else 1.5)
        ax.annotate(
            f"{og7}\n{code} ({name[:18]})\n{row['n_especies_red']} sp.",
            xy=(x_pt, y_pt), xytext=(x_txt, y_txt), fontsize=7, color=colors_top15[i],
            arrowprops=dict(arrowstyle="->", color=colors_top15[i], lw=1.2),
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=colors_top15[i], alpha=0.85),
        )
        ax.plot(x_pt, y_pt, "o", color=colors_top15[i], markersize=8, zorder=5)

    ax.set_xlabel("Number of structures per ligand (log scale)", fontsize=12)
    ax.set_ylabel("Number of ligands with that frequency (log scale)", fontsize=12)
    ax.set_title(f"Full distribution of {len(df_ligs):,} unique biological ligands in BioLiP\n"
                 "with the top 15 OG7 candidates located on it", fontweight="bold", fontsize=12)
    ax.grid(True, alpha=0.3, which="both")
    ax.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/06_loglog_anotado.png", dpi=150, bbox_inches="tight")
    plt.close()

    # ── Plot 2: top 20 OG7, network vs. BioLiP species coverage ────────────
    top20 = df_cand_sorted.head(20).sort_values("n_especies_red", ascending=True)
    fig, ax = plt.subplots(figsize=(14, 9))
    y = np.arange(len(top20))
    h = 0.38
    bars_red = ax.barh(y + h / 2, top20["n_especies_red"], h, color="#2E75B6",
                        label="Species in the network (with or without a crystal)")
    bars_biolip = ax.barh(y - h / 2, top20["n_especies_biolip"], h, color="#70AD47",
                           label="Species with a BioLiP structure")
    for bar, row in zip(bars_red, top20.itertuples()):
        ax.text(bar.get_width() + 0.2, bar.get_y() + bar.get_height() / 2,
                 f"{int(row.n_especies_red)} sp. | {int(row.n_prot_red)} prot. | {row.prot_por_especie:.1f}x",
                 va="center", fontsize=7.5, color="#2E75B6")
    for bar, row in zip(bars_biolip, top20.itertuples()):
        ax.text(bar.get_width() + 0.2, bar.get_y() + bar.get_height() / 2,
                 f"{int(row.n_uniprot_biolip)} struct. | res={row.resolucion_media:.2f}Å",
                 va="center", fontsize=7.5, color="#70AD47")
    labels = [f"{row.og7}\n{int(row.n_ligandos_unicos) if not pd.isna(row.n_ligandos_unicos) else 0} distinct lig. | "
              f"rare: {str(row.ligando_raro or '')[:8]} ({str(row.nombre_ligando_raro or '')[:15]})\n"
              f"EC: {str(row.ec_numbers)[:25]}" for row in top20.itertuples()]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.set_xlabel("Number of species", fontsize=11)
    ax.set_title("Top 20 OG7 candidate groups for PocketVec\n"
                 "Blue = presence in the network | Green = BioLiP structures",
                 fontweight="bold", fontsize=12)
    ax.legend(fontsize=9)
    ax.set_xlim(0, top20["n_especies_red"].max() * 1.55)
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/07_top20_og7_final.png", dpi=150, bbox_inches="tight")
    plt.close()

    # ── Plot 3: ligands with >100 structures ─────────────────────────────
    df_big_sorted = df_big_ligs.sort_values("n_entradas", ascending=True)
    labels_lig = [f"{row['ligand_id']} — {str(row['ligand_name'])[:35]}" for _, row in df_big_sorted.iterrows()]
    values = list(df_big_sorted["n_entradas"])
    fig, ax = plt.subplots(figsize=(13, max(6, len(df_big_sorted) * 0.35)))
    bars = ax.barh(labels_lig, values, color=plt.cm.Blues(np.linspace(0.35, 0.85, len(labels_lig))))
    for bar, val in zip(bars, values):
        ax.text(bar.get_width() + 10, bar.get_y() + bar.get_height() / 2, f"{val:,}", va="center", fontsize=8)
    ax.set_xlabel("Number of structures in BioLiP", fontsize=11)
    ax.set_title("Biological ligands with more than 100 structures in BioLiP\n"
                 "(the ones that dominate the dataset)", fontweight="bold", fontsize=12)
    ax.tick_params(axis="y", labelsize=8)
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/08_ligandos_mas100.png", dpi=150, bbox_inches="tight")
    plt.close()

    print(f"\nAll outputs saved to {OUT_DIR}/")


if __name__ == "__main__":
    main()
