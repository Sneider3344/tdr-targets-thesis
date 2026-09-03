#!/usr/bin/env python3
"""
Cross-references BioLiP ligand-binding entries with the network's OG7
orthology groups.

Builds a UniProt -> OG7 map from the three ID-resolution paths used
throughout the project (direct UniProt species, Gene Mapper output for
OrthoMCL-native species, and the not-on-OrthoMCL species), then uses that
map to tag each BioLiP entry with its OG7 group.

Outputs (used by the next steps in the candidate-selection pipeline):
  - resumen_por_og7.tsv     : one row per OG7 group with a BioLiP structure
  - resumen_por_ec4.tsv     : same data grouped by EC number (4 levels)
  - distribucion_ligandos.tsv
  - plots 01-05 (species/EC/ligand/resolution overview)
"""
import os
import gzip
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ── Config ──────────────────────────────────────────────────────────────
BIOLIP_FILE = "BioLiP.txt.gz"
LIGAND_FILE = "ligand_summary.txt.gz"
OUT_DIR     = "biolip_exploracion"

SEQS_DIR      = "/home/ssneider/disco-big/ssneider-env/TDR_Targets7.1/actualizados_short/sequences"
MAPPER_DIR    = "/home/ssneider/disco-big/ssneider-env/TDR_Targets7.1/gene_mapper/Ortho_vs_Uniprot"
NOT_ORTHO_DIR = "/big/lab/mercedesdg/TDR_Targets_7/OrthoMCL/genomes_v7/not_on_orthomcl/new_genomes"

SPECIES_UNIPROT_DIRECT = {"atha", "dmel", "osat", "cele", "ddis", "ecol", "mtub"}
SPECIES_NOT_ON_ORTHOMCL = {"ovo", "egr", "kpm", "loa", "sao"}
ORTHO_FILE_NAMES = {"egr": "egr_OrthoMCL_asignation_evaluation.tsv"}

LIGANDS_EXCLUDE = {
    "HOH", "DOD", "WAT", "SO4", "PO4", "GOL", "EDO", "PEG", "BME",
    "ACT", "ACE", "FMT", "MPD", "MES", "TAR", "CIT", "TRS", "EPE",
    "MG", "ZN", "CA", "FE", "MN", "NA", "K", "CL", "CU", "NI", "CO",
    "SE", "CD", "HG", "PT", "AU", "AG", "IOD", "BR", "ION", "FLC", "DMS",
}

COLS = [
    "pdb_id", "chain", "resolution", "binding_site",
    "ligand_id", "ligand_chain", "ligand_serial",
    "binding_residues_pdb", "binding_residues_renum",
    "catalytic_pdb", "catalytic_renum",
    "ec_number", "go_terms",
    "affinity_manual", "affinity_moad", "affinity_pdbbind", "affinity_bindingdb",
    "uniprot_id", "pubmed_id", "ligand_seqnum", "receptor_seq",
]

os.makedirs(OUT_DIR, exist_ok=True)


# ── Block 1: build the UniProt -> OG7 map ─────────────────────────────────
def build_uniprot_to_og7():
    uniprot_to_og = {}

    # 1A: species with direct UniProt IDs in the header
    print("[1A] Species with direct UniProt IDs...")
    for species in SPECIES_UNIPROT_DIRECT:
        fasta_candidates = list(Path(SEQS_DIR).glob(f"{species}_aa_seqs_OrthoMCL-7.fasta")) + \
                            list(Path(SEQS_DIR).glob(f"{species}_protein.fasta"))
        if not fasta_candidates:
            print(f"  WARNING: no FASTA found for {species}")
            continue
        count = 0
        with open(fasta_candidates[0]) as f:
            for line in f:
                if line.startswith(">"):
                    parts = line.strip().lstrip(">").split()
                    if len(parts) >= 2:
                        raw_id = parts[0]
                        uniprot = raw_id.split("|")[1] if "|" in raw_id else raw_id
                        og = next((p for p in parts if p.startswith("OG7_")), None)
                        if og and uniprot:
                            uniprot_to_og[uniprot] = og
                            count += 1
        print(f"  {species}: {count:,} proteins mapped")

    # 1B: species with a native ID -> mapped via Gene Mapper (mapped_clean.csv)
    print("\n[1B] Species mapped through Gene Mapper...")
    mapper_path = Path(MAPPER_DIR)
    all_mapper_species = {p.name for p in mapper_path.iterdir() if p.is_dir()}
    species_mapper = all_mapper_species - SPECIES_UNIPROT_DIRECT - SPECIES_NOT_ON_ORTHOMCL
    if "tcr" in species_mapper:
        species_mapper.discard("tcr")
        species_mapper.add("tcru")

    for species in sorted(species_mapper):
        species_dir = "tcru" if species == "tcru" else species
        mapped_file = mapper_path / species_dir / "mapped_clean.csv"
        if not mapped_file.exists():
            mapped_file = mapper_path / "tcr" / "mapped_clean.csv"
        if not mapped_file.exists():
            print(f"  WARNING: no mapped_clean.csv found for {species}")
            continue

        df_map = pd.read_csv(mapped_file, dtype=str)
        native_to_uniprot = dict(zip(df_map["OrthoMCL_ID"].str.strip(),
                                      df_map["Uniprot_ID"].str.strip()))

        fasta_candidates = list(Path(SEQS_DIR).glob(f"{species}_aa_seqs_OrthoMCL-7.fasta")) + \
                            list(Path(SEQS_DIR).glob(f"{species}_protein.fasta"))
        if not fasta_candidates:
            print(f"  WARNING: no FASTA found for {species}")
            continue

        count, not_found = 0, 0
        with open(fasta_candidates[0]) as f:
            for line in f:
                if line.startswith(">"):
                    parts = line.strip().lstrip(">").split()
                    if len(parts) >= 2:
                        raw_id = parts[0]
                        native_id = raw_id.split("|")[1] if "|" in raw_id else raw_id
                        og = next((p for p in parts if p.startswith("OG7_")), None)
                        uniprot = native_to_uniprot.get(native_id)
                        if og and uniprot:
                            uniprot_to_og[uniprot] = og
                            count += 1
                        elif og:
                            not_found += 1
        print(f"  {species}: {count:,} mapped, {not_found:,} without a UniProt ID in mapped_clean")

    # 1C: species not on OrthoMCL
    print("\n[1C] Species not on OrthoMCL...")
    for species in sorted(SPECIES_NOT_ON_ORTHOMCL):
        filename = ORTHO_FILE_NAMES.get(species, f"{species}_orthomcl7.1.txt")
        ortho_file = Path(NOT_ORTHO_DIR) / species / filename
        if not ortho_file.exists():
            print(f"  WARNING: {ortho_file} not found")
            continue
        mapped_file = mapper_path / species / "mapped_clean.csv"
        if not mapped_file.exists():
            print(f"  WARNING: no mapped_clean.csv found for {species}")
            continue

        df_map = pd.read_csv(mapped_file, dtype=str)
        native_to_uniprot = dict(zip(df_map["OrthoMCL_ID"].str.strip(),
                                      df_map["Uniprot_ID"].str.strip()))

        df_ortho = pd.read_csv(ortho_file, sep="\t", dtype=str)
        df_ortho.columns = df_ortho.columns.str.strip()
        col_query, col_og = df_ortho.columns[0], df_ortho.columns[2]

        count, not_found = 0, 0
        for _, row in df_ortho.iterrows():
            query_id = str(row[col_query]).strip()
            og = str(row[col_og]).strip()
            if not og.startswith("OG7_"):
                continue
            # egr has no species|ID prefix, the ID is used directly
            native_id = query_id.split("|")[1] if "|" in query_id else query_id
            uniprot = native_to_uniprot.get(native_id)
            if uniprot:
                uniprot_to_og[uniprot] = og
                count += 1
            else:
                not_found += 1
        print(f"  {species}: {count:,} mapped, {not_found:,} without a UniProt ID in mapped_clean")

    print(f"\nTotal UniProt -> OG7 entries mapped: {len(uniprot_to_og):,}")
    return uniprot_to_og, species_mapper, mapper_path


def build_uniprot_to_species(species_mapper, mapper_path):
    uniprot_to_species = {}
    for species in SPECIES_UNIPROT_DIRECT:
        fasta_candidates = list(Path(SEQS_DIR).glob(f"{species}_aa_seqs_OrthoMCL-7.fasta")) + \
                            list(Path(SEQS_DIR).glob(f"{species}_protein.fasta"))
        if not fasta_candidates:
            continue
        with open(fasta_candidates[0]) as f:
            for line in f:
                if line.startswith(">"):
                    parts = line.strip().lstrip(">").split()
                    if parts:
                        raw_id = parts[0]
                        uniprot = raw_id.split("|")[1] if "|" in raw_id else raw_id
                        uniprot_to_species[uniprot] = species

    for species in sorted(species_mapper):
        species_dir = "tcru" if species == "tcru" else species
        mapped_file = mapper_path / species_dir / "mapped_clean.csv"
        if not mapped_file.exists():
            mapped_file = mapper_path / "tcr" / "mapped_clean.csv"
        if not mapped_file.exists():
            continue
        df_map = pd.read_csv(mapped_file, dtype=str)
        for uniprot in df_map["Uniprot_ID"].dropna().str.strip():
            uniprot_to_species[uniprot] = species

    for species in sorted(SPECIES_NOT_ON_ORTHOMCL):
        mapped_file = mapper_path / species / "mapped_clean.csv"
        if not mapped_file.exists():
            continue
        df_map = pd.read_csv(mapped_file, dtype=str)
        for uniprot in df_map["Uniprot_ID"].dropna().str.strip():
            uniprot_to_species[uniprot] = species

    return uniprot_to_species


def main():
    uniprot_to_og, species_mapper, mapper_path = build_uniprot_to_og7()

    # ── Block 2: load and clean BioLiP ────────────────────────────────────
    print("\nLoading BioLiP...")
    opener = gzip.open if BIOLIP_FILE.endswith(".gz") else open
    with opener(BIOLIP_FILE, "rt", encoding="utf-8", errors="replace") as f:
        df = pd.read_csv(f, sep="\t", header=None, names=COLS, low_memory=False)
    print(f"Total BioLiP entries: {len(df):,}")

    df["resolution"] = pd.to_numeric(df["resolution"], errors="coerce")
    df["is_xray"] = df["resolution"] > 0
    df["uniprot_id"] = df["uniprot_id"].astype(str).str.strip()
    df["uniprot_first"] = df["uniprot_id"].str.split(",").str[0].str.strip()
    df["has_uniprot"] = df["uniprot_first"].str.match(r"^[A-Z][0-9][A-Z0-9]{3}[0-9]$")
    df["ligand_id"] = df["ligand_id"].astype(str).str.strip().str.upper()
    df["ligand_bio"] = ~df["ligand_id"].isin(LIGANDS_EXCLUDE)
    df["ec_number"] = df["ec_number"].astype(str).str.strip()

    df_clean = df[df["is_xray"] & df["has_uniprot"] & df["ligand_bio"]].copy()
    print(f"With X-ray + UniProt + biological ligand: {len(df_clean):,}")
    print(f"Unique UniProt IDs: {df_clean['uniprot_first'].nunique():,}")

    lig_names = {}
    if os.path.exists(LIGAND_FILE):
        lig_cols = ["lig_id", "formula", "inchi", "inchikey", "smiles", "name", "chebi", "drugbank", "zinc"]
        opener2 = gzip.open if LIGAND_FILE.endswith(".gz") else open
        with opener2(LIGAND_FILE, "rt", encoding="utf-8", errors="replace") as f:
            df_lig = pd.read_csv(f, sep="\t", header=None, names=lig_cols, low_memory=False,
                                  on_bad_lines="warn", engine="python")
        df_lig["lig_id"] = df_lig["lig_id"].str.strip().str.upper()
        df_lig["short_name"] = df_lig["name"].astype(str).str.split(";").str[0].str.strip()
        lig_names = dict(zip(df_lig["lig_id"], df_lig["short_name"]))

    df_clean["ligand_name"] = df_clean["ligand_id"].map(lig_names).fillna(df_clean["ligand_id"])

    # ── Block 3: cross-reference with OrthoMCL ────────────────────────────
    print("\nCross-referencing BioLiP with OrthoMCL...")
    df_clean["og7"] = df_clean["uniprot_first"].map(uniprot_to_og)
    df_with_og = df_clean[df_clean["og7"].notna()].copy()

    print(f"Entries with an OG7 group (in the network): {len(df_with_og):,}")
    print(f"Entries without an OG7 group (outside the network): {df_clean['og7'].isna().sum():,}")
    print(f"Unique OG7 groups in BioLiP: {df_with_og['og7'].nunique():,}")

    # ── Block 4: summary tables ────────────────────────────────────────────
    total_entries = len(df_clean)
    uniprot_to_species = build_uniprot_to_species(species_mapper, mapper_path)
    df_with_og["species"] = df_with_og["uniprot_first"].map(uniprot_to_species)

    og_summary = (
        df_with_og.groupby("og7")
        .agg(
            n_entries=("pdb_id", "count"),
            n_unique_pdb=("pdb_id", "nunique"),
            n_unique_uniprot=("uniprot_first", "nunique"),
            n_species=("species", "nunique"),
            species_list=("species", lambda x: ", ".join(sorted(x.dropna().unique()))),
            n_unique_ligands=("ligand_id", "nunique"),
            mean_resolution=("resolution", "mean"),
            top_ligands=("ligand_id", lambda x: ", ".join(x.value_counts().head(3).index)),
            ec_numbers=("ec_number", lambda x: ", ".join(x[x.str.match(r"^\d")].unique()[:3])),
        )
        .reset_index()
        .sort_values(["n_species", "n_unique_uniprot"], ascending=False)
    )
    og_summary["pct_entries"] = (og_summary["n_entries"] / total_entries * 100).round(2)
    # rename to the column names used by the downstream candidate-selection scripts
    og_summary_out = og_summary.rename(columns={
        "n_unique_pdb": "pdb_unicos", "n_unique_uniprot": "uniprot_unicos",
        "n_species": "n_especies", "species_list": "especies_lista",
        "n_unique_ligands": "ligandos_unicos", "mean_resolution": "resolucion_media",
        "top_ligands": "top_ligandos", "n_entries": "n_entradas",
    })
    og_summary_out.to_csv(f"{OUT_DIR}/resumen_por_og7.tsv", sep="\t", index=False)
    print(f"\nTop 20 OG7 groups (by number of species):")
    print(og_summary[["og7", "n_species", "n_unique_uniprot", "n_unique_pdb", "mean_resolution",
                       "top_ligands", "ec_numbers"]].head(20).to_string(index=False))

    df_clean["ec_4lvl"] = df_clean["ec_number"].str.extract(r"^(\d+\.\d+\.\d+\.\d+)")[0]
    ec_summary = (
        df_clean[df_clean["ec_4lvl"].notna()]
        .groupby("ec_4lvl")
        .agg(n_entries=("pdb_id", "count"), n_unique_pdb=("pdb_id", "nunique"),
             n_unique_uniprot=("uniprot_first", "nunique"))
        .reset_index()
        .sort_values("n_unique_uniprot", ascending=False)
    )
    ec_summary["pct_entries"] = (ec_summary["n_entries"] / total_entries * 100).round(2)
    ec_summary["pct_cumulative"] = ec_summary["pct_entries"].cumsum().round(2)
    ec_summary.to_csv(f"{OUT_DIR}/resumen_por_ec4.tsv", sep="\t", index=False)
    print(f"\nTop 25 EC families account for: {ec_summary.head(25)['pct_entries'].sum():.1f}% of the dataset")

    lig_dist = df_clean["ligand_id"].value_counts().reset_index()
    lig_dist.columns = ["ligand_id", "n_entries"]
    lig_dist["ligand_name"] = lig_dist["ligand_id"].map(lig_names).fillna(lig_dist["ligand_id"])
    lig_dist["pct"] = (lig_dist["n_entries"] / total_entries * 100).round(3)
    lig_dist.to_csv(f"{OUT_DIR}/distribucion_ligandos.tsv", sep="\t", index=False)
    print(f"\nUnique biological ligands: {len(lig_dist):,}")

    # ── Block 5: overview plots ─────────────────────────────────────────────
    sns.set_theme(style="whitegrid", font_scale=1.1)
    palette = sns.color_palette("Blues_r", 10)

    fig, ax = plt.subplots(figsize=(13, 10))
    top30_og = og_summary.head(30).sort_values("n_species", ascending=True)
    norm = plt.Normalize(top30_og["n_species"].min(), top30_og["n_species"].max())
    colors_og = plt.cm.Blues(norm(top30_og["n_species"]))
    bars = ax.barh(top30_og["og7"], top30_og["n_unique_uniprot"], color=colors_og)
    for bar, row in zip(bars, top30_og.itertuples()):
        ax.text(bar.get_width() + 0.3, bar.get_y() + bar.get_height() / 2,
                 f"{row.n_unique_uniprot} UniProt | {row.n_species} sp. | {row.top_ligands}",
                 va="center", fontsize=7.5)
    ax.set_xlabel("Unique UniProt IDs")
    ax.set_title("Top 30 orthology groups (OG7) with the most BioLiP structures\n"
                 "(darker bars = more species represented)", fontweight="bold")
    sm = plt.cm.ScalarMappable(cmap="Blues", norm=norm)
    sm.set_array([])
    plt.colorbar(sm, ax=ax, label="Number of species", shrink=0.6)
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/01_top30_og7.png", dpi=150)
    plt.close()

    fig, ax = plt.subplots(figsize=(12, 9))
    top30_ec = ec_summary.head(30).sort_values("n_unique_uniprot", ascending=True)
    bars = ax.barh(top30_ec["ec_4lvl"], top30_ec["n_unique_uniprot"], color=palette[3])
    for bar, row in zip(bars, top30_ec.itertuples()):
        ax.text(bar.get_width() + 0.3, bar.get_y() + bar.get_height() / 2,
                 f"{row.n_unique_uniprot}  ({row.pct_entries}%)", va="center", fontsize=8)
    ax.set_xlabel("Unique UniProt IDs")
    ax.set_title("Top 30 enzyme families (EC, 4 levels) in BioLiP\n(% = fraction of total entries)",
                 fontweight="bold")
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/02_top30_ec4_con_pct.png", dpi=150)
    plt.close()

    fig, axes = plt.subplots(1, 2, figsize=(15, 7))
    top40 = lig_dist.head(40).sort_values("n_entries", ascending=True)
    labels = [f"{row.ligand_id} ({row.ligand_name[:18]})" for row in top40.itertuples()]
    axes[0].barh(labels, top40["n_entries"], color=palette[2])
    axes[0].set_xlabel("Number of entries")
    axes[0].set_title("Top 40 biological ligands", fontweight="bold")
    axes[0].tick_params(axis="y", labelsize=7)

    freq_dist = lig_dist["n_entries"].value_counts().sort_index()
    axes[1].loglog(freq_dist.index, freq_dist.values, "o", color=palette[1], markersize=4, alpha=0.7)
    axes[1].set_xlabel("Structures per ligand (log)")
    axes[1].set_ylabel("Number of ligands (log)")
    axes[1].set_title(f"Full distribution — {len(lig_dist):,} unique ligands\n(log-log scale)",
                       fontweight="bold")
    axes[1].grid(True, alpha=0.3)
    plt.suptitle("Distribution of biological ligands in BioLiP", fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/03_distribucion_ligandos.png", dpi=150, bbox_inches="tight")
    plt.close()

    fig, ax = plt.subplots(figsize=(9, 4))
    res_data = df_clean[df_clean["resolution"].between(0.5, 4.0)]["resolution"]
    ax.hist(res_data, bins=60, color=palette[2], edgecolor="white")
    ax.axvline(2.0, color="red", linestyle="--", label="2.0 Å (very good)")
    ax.axvline(2.5, color="orange", linestyle="--", label="2.5 Å (acceptable)")
    ax.axvline(3.5, color="gray", linestyle="--", label="3.5 Å (limit)")
    ax.set_xlabel("Resolution (Å)")
    ax.set_ylabel("Entries")
    ax.set_title("Crystallographic resolution distribution — BioLiP", fontweight="bold")
    ax.legend()
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/04_resolucion.png", dpi=150)
    plt.close()

    n_og_with_structure = df_with_og["og7"].nunique()
    n_og_in_network = len(set(uniprot_to_og.values()))
    fig, ax = plt.subplots(figsize=(6, 6))
    sizes = [n_og_with_structure, max(0, n_og_in_network - n_og_with_structure)]
    labels = [f"With a BioLiP\nstructure\n({n_og_with_structure:,})",
              f"Without a BioLiP\nstructure\n({max(0, n_og_in_network - n_og_with_structure):,})"]
    ax.pie(sizes, labels=labels, colors=[palette[1], "#DDDDDD"],
           autopct="%1.1f%%", startangle=90, textprops={"fontsize": 11})
    ax.set_title("Network OrthoMCL groups\nwith/without a BioLiP structure", fontweight="bold")
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/05_og7_cobertura.png", dpi=150)
    plt.close()

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Total BioLiP entries (X-ray, biological ligand, UniProt): {total_entries:,}")
    print(f"Unique biological ligands: {len(lig_dist):,}")
    print(f"OG7 groups with at least 1 BioLiP structure: {n_og_with_structure:,}")
    print(f"\nOutputs saved to: {OUT_DIR}/")


if __name__ == "__main__":
    main()
