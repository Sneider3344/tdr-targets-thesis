#!/usr/bin/env python3
"""
Pillars 1 + 2: correspondence between InterPro Family and InterPro Domain.

Pillar 1: coverage -- of the proteins annotated with InterPro Family,
          how many also have an associated InterPro Domain.

Pillar 2a: heterogeneity BETWEEN members of the same Family -- within a
           Family, do the members (potentially from different species)
           share the same Domain, have different Domains, or none at all?

Pillar 2b: positional overlap WITHIN the same protein -- when a protein
           has both a Family and a Domain, do the coordinates match, is
           the Domain nested inside the Family, or do they not overlap?

Inputs (already downloaded, nothing new to fetch here):
  - informacion_interpro.txt    : type of each InterPro entry (Domain/Family/etc.)
  - interproscan_corridas/*.tsv : per-species InterProScan TSVs (all 29)
"""
import os
import glob
from collections import defaultdict, Counter
import pandas as pd

# ── Paths ─────────────────────────────────────────────────────────────────
INTERPRO_INFO      = "/home/ssneider/disco-big/ssneider-env/TDR_Targets7.1/interpro/informacion_interpro.txt"
INTERPROSCAN_DIR   = "/home/ssneider/disco-big/ssneider-env/TDR_Targets7.1/interproscan_corridas"
OUTPUT_DIR         = "/home/ssneider/disco-big/ssneider-env/TDR_Targets7.1/interpro/family_domain_analysis"
GENE_MAPPER_DIR    = "/home/ssneider/disco-big/ssneider-env/TDR_Targets7.1/gene_mapper/Ortho_vs_Uniprot"
PFAM_AF_DIR        = "/home/ssneider/disco-big/ssneider-env/TDR_Targets7.1/PFAM_vs_alphafold/Species"

# Same as in find_lost_proteins.py: species extraction from the filename
FILENAME_TO_SPECIES = {
    "echinococcus": "egr",
}

# Same ID-mapping groups used across the rest of the pipeline (already fixed)
SPECIES_MAPPED_CLEAN = {
    "bmaa", "egr", "hsap", "kpm", "ldon", "lmaj", "loa", "mmus", "ovo", "pvip", "sao", "tcru",
    "calb", "ehia", "gass", "gmur", "pfal", "scer", "tbrt", "tgon", "tvag"
}
SPECIES_NO_MAPPING = {"atha", "dmel", "osat", "cele", "ddis", "ecol", "mtub"}
SPECIES_TRY_BOTH   = {"tcru", "tcr"}

HOMOGENEOUS_THRESHOLD = 0.90  # fraction of members that must share the majority domain


# ── 0. ID mapping (same as the rest of the pipeline) ──────────────────────
def load_mapped_clean(species):
    path = os.path.join(GENE_MAPPER_DIR, species, "mapped_clean.csv")
    mapping = {}
    if not os.path.exists(path):
        return mapping
    with open(path, newline="") as fh:
        for i, line in enumerate(fh):
            cols = line.rstrip("\n").split(",")
            if len(cols) < 2 or i == 0:
                continue
            mapping[cols[0].strip()] = cols[1].strip()
    return mapping


def load_idmapping(species):
    path = os.path.join(PFAM_AF_DIR, species, f"idmapping_{species}.tsv")
    mapping = {}
    if not os.path.exists(path):
        return mapping
    with open(path, newline="") as fh:
        for i, line in enumerate(fh):
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 2 or i == 0:
                continue
            mapping[cols[0].strip()] = cols[1].strip()
    return mapping


def load_tsv_file(path):
    mapping = {}
    if not os.path.exists(path):
        return mapping
    delim = "," if path.endswith(".csv") else "\t"
    with open(path, newline="") as fh:
        for i, line in enumerate(fh):
            cols = line.rstrip("\n").split(delim)
            if len(cols) < 2 or i == 0:
                continue
            mapping[cols[0].strip()] = cols[1].strip()
    return mapping


def load_combined(species):
    TCRU_FILES = {
        "tcru": os.path.join(PFAM_AF_DIR, "tcru",     "dm28c_mapped.tsv"),
        "tcr":  os.path.join(PFAM_AF_DIR, "tcrc-l-n", "CL-brenner_mapped.tsv"),
    }
    mapping = load_idmapping(species)
    mapping.update(load_mapped_clean(species))
    mapping.update(load_tsv_file(TCRU_FILES[species]))
    return mapping


def get_id_mapping(species):
    if species in SPECIES_NO_MAPPING:
        return {}
    elif species in SPECIES_TRY_BOTH:
        return load_combined(species)
    elif species in SPECIES_MAPPED_CLEAN:
        return load_mapped_clean(species)
    else:
        return load_idmapping(species)


def resolve_uniprot(raw_id, mapping):
    clean_id = raw_id.split("|", 1)[1] if "|" in raw_id else raw_id
    if not mapping:
        return clean_id
    return mapping.get(clean_id, clean_id)  # no mapping found, keep the raw id


# ── 1. InterPro entry types ─────────────────────────────────────────────
def load_interpro_types(path):
    types = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("ENTRY_AC"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                types[parts[0]] = parts[1]
    return types


# ── 2. Parse a single species TSV ───────────────────────────────────────
def parse_species_tsv(tsv_path, interpro_types, global_type_counter):
    """
    Returns:
      proteins : dict protein_id -> {
          "family": [(ipr_acc, start, end), ...],
          "domain": [(ipr_acc, start, end), ...],
      }
      protein_lengths : dict protein_id -> full sequence length (cols[2])
    Only considers Pfam hits (same as the rest of the pipeline). Along the
    way, tallies in global_type_counter how many hits fall under each
    InterPro entry type (Domain, Family, Repeat, Homologous_superfamily, etc.)
    """
    proteins = defaultdict(lambda: {"family": [], "domain": []})
    protein_lengths = {}
    with open(tsv_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 12:
                continue
            if cols[3] != "Pfam":
                continue
            prot_id = cols[0]
            protein_lengths[prot_id] = int(cols[2])
            start = int(cols[6])
            end = int(cols[7])
            ipr_acc = cols[11]
            if not ipr_acc or ipr_acc == "-":
                global_type_counter["No associated IPR"] += 1
                continue
            itype = interpro_types.get(ipr_acc, "Unknown")
            global_type_counter[itype] += 1
            if itype == "Family":
                proteins[prot_id]["family"].append((ipr_acc, start, end))
            elif itype == "Domain":
                proteins[prot_id]["domain"].append((ipr_acc, start, end))
    return proteins, protein_lengths


# ── 3. Main ───────────────────────────────────────────────────────────────
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("Loading InterPro entry types...")
    interpro_types = load_interpro_types(INTERPRO_INFO)

    tsv_files = sorted(glob.glob(os.path.join(INTERPROSCAN_DIR, "*ips_data.tsv")))
    print(f"{len(tsv_files)} files found in {INTERPROSCAN_DIR}\n")

    # family_acc -> list of (species, protein_id, set_of_domain_accs_for_that_protein)
    family_members = defaultdict(list)
    # positional overlap rows (Pillar 2b)
    overlap_rows = []

    # ── Pillar 1: coverage, accumulated per species ───────────────────────
    coverage_rows = []
    global_type_counter = Counter()
    all_windows_rows = []  # every Family/Domain window, needed for Pillar 3
    protein_lengths_rows = []  # protein length, needed for the enrichment test

    for tsv_path in tsv_files:
        filename        = os.path.basename(tsv_path)
        filename_prefix = filename.split("_")[0]
        species         = FILENAME_TO_SPECIES.get(filename_prefix, filename_prefix)

        proteins, protein_lengths = parse_species_tsv(tsv_path, interpro_types, global_type_counter)
        id_mapping = get_id_mapping(species)

        n_with_family = 0
        n_family_and_domain = 0

        for prot_id, d in proteins.items():
            fam_accs = {acc for (acc, s, e) in d["family"]}
            dom_accs = {acc for (acc, s, e) in d["domain"]}
            uniprot_id = resolve_uniprot(prot_id, id_mapping)

            protein_lengths_rows.append({
                "species": species, "protein_id": prot_id, "uniprot_id": uniprot_id,
                "length": protein_lengths.get(prot_id),
            })

            # full table of windows (every Family and every Domain, whether
            # or not the protein has both types) -- used by Pillar 3 / BioLiP
            for (acc, s, e) in d["family"]:
                all_windows_rows.append({
                    "species": species, "protein_id": prot_id, "uniprot_id": uniprot_id,
                    "entry_type": "Family", "ipr_acc": acc, "start": s, "end": e,
                })
            for (acc, s, e) in d["domain"]:
                all_windows_rows.append({
                    "species": species, "protein_id": prot_id, "uniprot_id": uniprot_id,
                    "entry_type": "Domain", "ipr_acc": acc, "start": s, "end": e,
                })

            if fam_accs:
                n_with_family += 1
                if dom_accs:
                    n_family_and_domain += 1

            # accumulate for Pillar 2a (once per protein x family, not per hit)
            for fam_acc in fam_accs:
                family_members[fam_acc].append((species, prot_id, dom_accs))

            # Pillar 2b: positional overlap Family <-> Domain, same protein
            for (fam_acc, fs, fe) in d["family"]:
                for (dom_acc, ds, de) in d["domain"]:
                    inter_start = max(fs, ds)
                    inter_end = min(fe, de)
                    inter_len = max(0, inter_end - inter_start + 1)
                    dom_len = de - ds + 1
                    fam_len = fe - fs + 1
                    overlap_rows.append({
                        "species": species, "protein_id": prot_id, "uniprot_id": uniprot_id,
                        "family_acc": fam_acc, "family_start": fs, "family_end": fe,
                        "domain_acc": dom_acc, "domain_start": ds, "domain_end": de,
                        "overlap_len": inter_len,
                        "pct_domain_covered": round(100 * inter_len / dom_len, 1) if dom_len else 0,
                        "pct_family_covered": round(100 * inter_len / fam_len, 1) if fam_len else 0,
                    })

        pct = 100 * n_family_and_domain / n_with_family if n_with_family else 0
        coverage_rows.append({
            "species": species,
            "n_with_family": n_with_family,
            "n_family_and_domain": n_family_and_domain,
            "pct_family_and_domain": round(pct, 1),
        })
        print(f"[{species}] proteins with Family: {n_with_family}  "
              f"with Family+Domain: {n_family_and_domain} ({pct:.1f}%)")

    # ── Save Pillar 1 ────────────────────────────────────────────────────
    cov_df = pd.DataFrame(coverage_rows)
    tot_fam = cov_df["n_with_family"].sum()
    tot_both = cov_df["n_family_and_domain"].sum()
    cov_df.loc[len(cov_df)] = {
        "species": "TOTAL",
        "n_with_family": tot_fam,
        "n_family_and_domain": tot_both,
        "pct_family_and_domain": round(100 * tot_both / tot_fam, 1) if tot_fam else 0,
    }
    cov_df.to_csv(os.path.join(OUTPUT_DIR, "pillar1_family_domain_coverage.tsv"), sep="\t", index=False)
    print("\n=== Pillar 1: Family -> Family+Domain coverage ===")
    print(cov_df.to_string(index=False))

    # ── Pillar 2a: heterogeneity per Family ───────────────────────────────
    print("\nComputing per-Family heterogeneity (Pillar 2a)...")
    fam_rows = []
    for fam_acc, members in family_members.items():
        n_members = len(members)
        n_with_domain = sum(1 for (_, _, doms) in members if doms)

        all_domains = Counter()
        for (_, _, doms) in members:
            all_domains.update(doms)

        if all_domains:
            majority_domain, n_majority = all_domains.most_common(1)[0]
        else:
            majority_domain, n_majority = None, 0

        n_distinct_domains = len(all_domains)
        pct_with_domain = 100 * n_with_domain / n_members if n_members else 0
        pct_majority = 100 * n_majority / n_members if n_members else 0

        if n_with_domain == 0:
            category = "No Domain"
        elif n_distinct_domains == 1 and pct_with_domain >= HOMOGENEOUS_THRESHOLD * 100:
            category = "Homogeneous"
        else:
            category = "Heterogeneous"

        fam_rows.append({
            "family_acc": fam_acc,
            "n_members": n_members,
            "n_with_domain": n_with_domain,
            "pct_with_domain": round(pct_with_domain, 1),
            "n_distinct_domains": n_distinct_domains,
            "majority_domain": majority_domain,
            "pct_majority_domain": round(pct_majority, 1),
            "category": category,
        })

    fam_df = pd.DataFrame(fam_rows)
    fam_df.to_csv(os.path.join(OUTPUT_DIR, "pillar2a_family_heterogeneity.tsv"), sep="\t", index=False)

    print("\n=== Pillar 2a: Family classification ===")
    counts = fam_df["category"].value_counts()
    pcts = (100 * fam_df["category"].value_counts(normalize=True)).round(1)
    for cat in counts.index:
        print(f"  {cat:<14} {counts[cat]:>6}  ({pcts[cat]}%)")

    # ── Pillar 2b: positional overlap ──────────────────────────────────────
    overlap_df = pd.DataFrame(overlap_rows)
    overlap_df.to_csv(os.path.join(OUTPUT_DIR, "pillar2b_family_domain_overlap.tsv"), sep="\t", index=False)

    print("\n=== Pillar 2b: positional overlap Family vs Domain (same protein) ===")
    if len(overlap_df):
        exact_match = ((overlap_df["family_start"] == overlap_df["domain_start"]) &
                       (overlap_df["family_end"] == overlap_df["domain_end"])).sum()
        nested = (overlap_df["pct_domain_covered"] >= 99).sum()
        no_overlap = (overlap_df["overlap_len"] == 0).sum()
        n = len(overlap_df)
        print(f"  Total Family-Domain pairs evaluated: {n}")
        print(f"  Exact coordinate match:              {exact_match} ({100*exact_match/n:.1f}%)")
        print(f"  Domain nested in Family (>=99% cov):  {nested} ({100*nested/n:.1f}%)")
        print(f"  No overlap (distinct regions):        {no_overlap} ({100*no_overlap/n:.1f}%)")

    print(f"\nOutputs saved to: {OUTPUT_DIR}")
    print("  - pillar1_family_domain_coverage.tsv")
    print("  - pillar2a_family_heterogeneity.tsv")
    print("  - pillar2b_family_domain_overlap.tsv")
    print("  - interpro_type_distribution.tsv")
    print("  - all_family_domain_windows.tsv")
    print("  - protein_lengths.tsv")

    # ── Full window table (for Pillar 3 / BioLiP) ─────────────────────────
    windows_df = pd.DataFrame(all_windows_rows)
    windows_df.to_csv(os.path.join(OUTPUT_DIR, "all_family_domain_windows.tsv"), sep="\t", index=False)

    lengths_df = pd.DataFrame(protein_lengths_rows)
    lengths_df.to_csv(os.path.join(OUTPUT_DIR, "protein_lengths.tsv"), sep="\t", index=False)

    # ── InterPro type distribution (justifies focusing on Domain+Family) ──
    print("\n=== Pfam hits by associated InterPro entry type ===")
    total_hits = sum(global_type_counter.values())
    type_rows = []
    for itype, count in global_type_counter.most_common():
        pct = 100 * count / total_hits if total_hits else 0
        type_rows.append({"interpro_type": itype, "n_hits": count, "pct": round(pct, 2)})
        print(f"  {itype:<20} {count:>10}  ({pct:.2f}%)")

    dom_fam_pct = 100 * (global_type_counter.get("Domain", 0) + global_type_counter.get("Family", 0)) / total_hits if total_hits else 0
    print(f"\n  --> Domain + Family account for {dom_fam_pct:.2f}% of Pfam hits with an associated IPR")

    type_df = pd.DataFrame(type_rows)
    type_df.to_csv(os.path.join(OUTPUT_DIR, "interpro_type_distribution.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    main()
