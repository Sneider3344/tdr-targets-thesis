#!/usr/bin/env python3
"""
Corrected candidate count ("segunda selección" in the thesis).

The initial filter counted, per OG7 group, how many species had their OWN
crystal structure in BioLiP. This undercounts real candidates: species
without a crystal can still be included in the analysis by transferring
the reference pocket's centroid via structural homology (US-align) from a
crystallized species in the same group.

This script instead counts, for each OG7 group, how many species of the
FULL network are present (reading the FASTA/orthology files directly,
independent of whether they have a BioLiP structure), and re-applies the
same three filters on that corrected count. At least one protein of the
group still has to have a crystal structure, to make the centroid
transfer possible.

Also includes a diagnostic step comparing network coverage vs. BioLiP
coverage for a handful of OG7 groups of interest, as a sanity check
before trusting the corrected numbers.
"""
from collections import defaultdict
from pathlib import Path

import pandas as pd

# ── Config ──────────────────────────────────────────────────────────────
SEQS_DIR      = "/home/ssneider/disco-big/ssneider-env/TDR_Targets7.1/actualizados_short/sequences"
MAPPER_DIR    = "/home/ssneider/disco-big/ssneider-env/TDR_Targets7.1/gene_mapper/Ortho_vs_Uniprot"
NOT_ORTHO_DIR = "/big/lab/mercedesdg/TDR_Targets_7/OrthoMCL/genomes_v7/not_on_orthomcl/new_genomes"
OUT_DIR       = "/home/ssneider/disco-big/ssneider-env/TDR_Targets7.1/PocketVec/analisis_biolip/biolip_exploracion"

SPECIES_UNIPROT_DIRECT = {"atha", "dmel", "osat", "cele", "ddis", "ecol", "mtub"}
SPECIES_NOT_ON_ORTHOMCL = {"ovo", "egr", "kpm", "loa", "sao"}
ORTHO_FILE_NAMES = {"egr": "egr_OrthoMCL_asignation_evaluation.tsv"}

# OG7 groups used for the sanity-check diagnostic below
OG7_DIAGNOSTIC = ["OG7_0001626", "OG7_0000433", "OG7_0000396", "OG7_0005325"]

LIGANDS_EXCLUDE_ANALYSIS = {
    "RNA", "DNA", "PEPTIDE", "OHX", "A", "U", "G", "C",
    "LLL", "PAR", "GLC", "PTD",
}


def build_og7_to_species():
    """og7 -> set of species present in the network (source of truth,
    independent of BioLiP), read directly from the FASTA / orthology
    assignment files."""
    og7_to_species = defaultdict(set)
    og7_to_n_proteins = defaultdict(int)

    print("Reading FASTAs to build og7 -> species (source of truth)...")

    mapper_path = Path(MAPPER_DIR)
    all_mapper_species = {p.name for p in mapper_path.iterdir() if p.is_dir()}
    species_mapper = all_mapper_species - SPECIES_UNIPROT_DIRECT - SPECIES_NOT_ON_ORTHOMCL
    if "tcr" in species_mapper:
        species_mapper.discard("tcr")
        species_mapper.add("tcru")

    for species in sorted(SPECIES_UNIPROT_DIRECT | species_mapper):
        fasta_candidates = list(Path(SEQS_DIR).glob(f"{species}_aa_seqs_OrthoMCL-7.fasta")) + \
                            list(Path(SEQS_DIR).glob(f"{species}_protein.fasta"))
        if not fasta_candidates:
            print(f"  WARNING: no FASTA found for {species}")
            continue
        with open(fasta_candidates[0]) as f:
            for line in f:
                if line.startswith(">"):
                    parts = line.strip().lstrip(">").split()
                    og = next((p for p in parts if p.startswith("OG7_")), None)
                    if og:
                        og7_to_species[og].add(species)
                        og7_to_n_proteins[og] += 1

    for species in sorted(SPECIES_NOT_ON_ORTHOMCL):
        filename = ORTHO_FILE_NAMES.get(species, f"{species}_orthomcl7.1.txt")
        ortho_file = Path(NOT_ORTHO_DIR) / species / filename
        if not ortho_file.exists():
            print(f"  WARNING: {ortho_file} not found")
            continue
        df_ortho = pd.read_csv(ortho_file, sep="\t", dtype=str)
        df_ortho.columns = df_ortho.columns.str.strip()
        col_og = df_ortho.columns[2]
        for og in df_ortho[col_og].dropna():
            og = str(og).strip()
            if og.startswith("OG7_"):
                og7_to_species[og].add(species)
                og7_to_n_proteins[og] += 1

    print(f"Total OG7 groups in the network: {len(og7_to_species):,}")
    return og7_to_species, og7_to_n_proteins


def ligand_ok(s):
    if not s or str(s) == "nan":
        return False
    return str(s).split(",")[0].strip().upper() not in LIGANDS_EXCLUDE_ANALYSIS


def main():
    og7_to_species, og7_to_n_proteins = build_og7_to_species()

    og_df = pd.read_csv(f"{OUT_DIR}/resumen_por_og7.tsv", sep="\t")
    og_dict = {row["og7"]: row for _, row in og_df.iterrows()}

    # ── Diagnostic: network coverage vs. BioLiP coverage for a few OG7 ────
    print("\n" + "=" * 70)
    print("DIAGNOSTIC: network species vs. BioLiP species for key OG7 groups")
    print("=" * 70)

    for og7 in OG7_DIAGNOSTIC:
        species_network = og7_to_species.get(og7, set())
        n_proteins_network = og7_to_n_proteins.get(og7, 0)

        if og7 in og_dict:
            row = og_dict[og7]
            species_biolip = set(str(row["especies_lista"]).split(", "))
            n_biolip = int(row["uniprot_unicos"])
            mean_res = float(row["resolucion_media"])
        else:
            species_biolip, n_biolip, mean_res = set(), 0, 0

        species_network_only = species_network - species_biolip
        species_biolip_only = species_biolip - species_network

        print(f"\n{og7}")
        print(f"  In the NETWORK ({len(species_network):>2} species, {n_proteins_network} proteins): "
              f"{', '.join(sorted(species_network))}")
        print(f"  In BioLiP      ({len(species_biolip):>2} species, {n_biolip} UniProt, res={mean_res:.2f}Å): "
              f"{', '.join(sorted(species_biolip))}")
        print(f"  In the network but WITHOUT a BioLiP structure: "
              f"{', '.join(sorted(species_network_only)) or 'none'}")
        if species_biolip_only:
            print(f"  ALERT: in BioLiP but not in the network: {', '.join(sorted(species_biolip_only))}")

    # ── Global table: og7 -> network species vs. BioLiP species ───────────
    print("\n" + "=" * 70)
    print("GLOBAL TABLE: real species coverage per OG7 (from FASTAs)")
    print("=" * 70)

    rows = []
    for og7, species_network in og7_to_species.items():
        if og7 in og_dict:
            row = og_dict[og7]
            species_biolip = set(str(row["especies_lista"]).split(", ")) if pd.notna(row["especies_lista"]) else set()
            n_uniprot = int(row["uniprot_unicos"])
            res = float(row["resolucion_media"])
            ligs = str(row["top_ligandos"])
            ec = str(row["ec_numbers"])
        else:
            species_biolip, n_uniprot, res, ligs, ec = set(), 0, 999, "", ""

        rows.append({
            "og7": og7,
            "n_especies_red": len(species_network),
            "n_prot_red": og7_to_n_proteins[og7],
            "n_especies_biolip": len(species_biolip),
            "n_uniprot_biolip": n_uniprot,
            "resolucion_media": res,
            "top_ligandos": ligs,
            "ec_numbers": ec,
            "especies_red": ", ".join(sorted(species_network)),
            "especies_biolip": ", ".join(sorted(species_biolip)),
        })

    df_global = pd.DataFrame(rows)
    df_global.to_csv(f"{OUT_DIR}/cobertura_especies_og7.tsv", sep="\t", index=False)
    print(f"Saved: {OUT_DIR}/cobertura_especies_og7.tsv")
    print(f"Total OG7 with at least 1 species in the network: {len(df_global):,}")

    # ── Corrected candidate list, using n_especies_red as the criterion ───
    df_cand = df_global[
        (df_global["resolucion_media"] < 2.5) &
        (df_global["resolucion_media"] > 0) &
        (df_global["top_ligandos"].apply(ligand_ok)) &
        (df_global["n_especies_red"] >= 3)
    ].copy()

    df_cand = df_cand.sort_values(["n_especies_red", "n_uniprot_biolip"], ascending=False)
    df_cand.to_csv(f"{OUT_DIR}/candidatos_pocketvec_v2.tsv", sep="\t", index=False)

    print(f"\nCandidates with the corrected criterion (species from FASTAs): {len(df_cand):,}")
    print("\nTop 20 candidates:")
    cols = ["og7", "n_especies_red", "n_prot_red", "n_especies_biolip",
            "n_uniprot_biolip", "resolucion_media", "top_ligandos", "ec_numbers"]
    print(df_cand[cols].head(20).to_string(index=False))


if __name__ == "__main__":
    main()
