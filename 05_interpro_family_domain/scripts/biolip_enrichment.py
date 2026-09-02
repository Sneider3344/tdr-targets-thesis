#!/usr/bin/env python3
"""
Ligand-binding enrichment of Domain vs. Family, corrected for ID-mapping
collisions.

Two related checks are run here:

1. General enrichment: across all proteins, are BioLiP contact residues
   over- or under-represented inside a Domain window, a Family window, or
   an unannotated region, relative to what window size alone would predict?

2. Restricted enrichment: same question, but limited to proteins that have
   BOTH a Domain and a Family annotation at once -- and only using
   uniprot_id values that map cleanly to a single raw protein_id. Some
   uniprot_id values collapse several distinct raw proteins into one
   (the same collision issue found earlier in the Gene Mapper step), so
   residues tied to an ambiguous uniprot_id are excluded from this second
   check: there is no way to know which of the colliding proteins the
   residue actually belongs to.

Inputs:
  - protein_lengths.tsv           : protein length per (species, protein_id)
  - all_family_domain_windows.tsv : Family/Domain windows per protein
  - biolip_domain_mapping.tsv     : BioLiP contact residues mapped to
                                     UniProt coordinates, with their
                                     Domain/Family/unannotated category
"""
import os
from collections import defaultdict
import pandas as pd
from scipy.stats import chisquare

BASE_DIR = "/home/ssneider/disco-big/ssneider-env/TDR_Targets7.1/interpro"
LENGTHS_TSV = os.path.join(BASE_DIR, "nuevo_analisis_tesis", "protein_lengths.tsv")
WINDOWS_TSV = os.path.join(BASE_DIR, "nuevo_analisis_tesis", "all_family_domain_windows.tsv")
BIOLIP_MAPPING_TSV = os.path.join(BASE_DIR, "nuevo_analisis_tesis/pilar3_biolip", "biolip_domain_mapping.tsv")
OUTPUT_DIR = os.path.join(BASE_DIR, "nuevo_analisis_tesis/pilar3_biolip")


def merge_intervals(intervals):
    """Collapses overlapping/adjacent (start, end) intervals and returns
    the total covered length."""
    if not intervals:
        return 0
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return sum(e - s + 1 for s, e in merged)


def main():
    print("Loading Family/Domain windows...")
    windows_df = pd.read_csv(WINDOWS_TSV, sep="\t")

    # ── 1. Flag which uniprot_id values are "clean" vs "colliding" ────────
    # A uniprot_id is clean when it maps 1-to-1 to a single (species,
    # protein_id) pair; it's colliding when several raw proteins share it.
    uid_to_proteins = windows_df.groupby("uniprot_id")[["species", "protein_id"]].apply(
        lambda d: set(map(tuple, d.values))
    )
    clean_uids = {uid for uid, prots in uid_to_proteins.items() if len(prots) == 1}
    colliding_uids = {uid for uid, prots in uid_to_proteins.items() if len(prots) > 1}
    print(f"Clean uniprot_id values (1-to-1): {len(clean_uids)}")
    print(f"Colliding uniprot_id values (excluded from the analysis): {len(colliding_uids)}")

    # ── 2. Domain/Family coverage per (species, protein_id) -- unambiguous ─
    print("Loading protein lengths...")
    lengths_df = pd.read_csv(LENGTHS_TSV, sep="\t")
    length_by_protein = {}
    for row in lengths_df.itertuples(index=False):
        if pd.notna(row.length):
            length_by_protein[(row.species, row.protein_id)] = int(row.length)

    fam_intervals = defaultdict(list)
    dom_intervals = defaultdict(list)
    for row in windows_df.itertuples(index=False):
        key = (row.species, row.protein_id)
        if row.entry_type == "Family":
            fam_intervals[key].append((int(row.start), int(row.end)))
        elif row.entry_type == "Domain":
            dom_intervals[key].append((int(row.start), int(row.end)))

    coverage_by_protein = {}  # (species, protein_id) -> (dom_cov, fam_cov, length, uniprot_id)
    uid_by_protein = windows_df.drop_duplicates(["species", "protein_id"]).set_index(
        ["species", "protein_id"])["uniprot_id"].to_dict()

    all_keys = set(fam_intervals) | set(dom_intervals)
    for key in all_keys:
        plen = length_by_protein.get(key)
        if not plen:
            continue
        dom_cov = min(merge_intervals(dom_intervals.get(key, [])), plen)
        fam_cov = min(merge_intervals(fam_intervals.get(key, [])), plen)
        coverage_by_protein[key] = (dom_cov, fam_cov, plen, uid_by_protein.get(key))

    # ── 3. Restrict to proteins with Domain AND Family, AND a clean uniprot_id ─
    coverage_both_clean = {
        key: (d, f, p) for key, (d, f, p, uid) in coverage_by_protein.items()
        if d > 0 and f > 0 and uid in clean_uids
    }
    print(f"\nProteins with Domain AND Family AND a clean uniprot_id: {len(coverage_both_clean)}")

    # uniprot_id -> (species, protein_id), only for the clean ones used here
    uid_to_key_clean = {
        uid_by_protein[key]: key for key in coverage_both_clean.keys()
    }

    # ── 4. Load BioLiP and keep only residues tied to clean uniprot_id values ─
    print("Loading BioLiP...")
    biolip_df = pd.read_csv(BIOLIP_MAPPING_TSV, sep="\t")

    n_total = len(biolip_df)
    biolip_clean = biolip_df[biolip_df["uniprot_id"].isin(uid_to_key_clean.keys())]
    n_excluded = n_total - len(biolip_clean)
    print(f"Residues excluded due to a colliding/ambiguous uniprot_id: {n_excluded}")

    # ── 5. Enrichment, same logic as before but over the clean subset ─────
    expected_domain = expected_family = 0.0
    observed_domain = observed_family = 0

    for row in biolip_clean.itertuples(index=False):
        if row.uniprot_id not in uid_to_key_clean:
            continue
        dom_cov, fam_cov, plen = coverage_both_clean[uid_to_key_clean[row.uniprot_id]]
        if row.category not in ("Domain", "Family"):
            continue  # "unannotated" doesn't apply within the both-categories universe
        expected_domain += dom_cov / plen
        expected_family += fam_cov / plen
        if row.category == "Domain":
            observed_domain += 1
        else:
            observed_family += 1

    print(f"\nBinding residues used (Domain+Family, clean uniprot_id): "
          f"{observed_domain + observed_family}")

    summary = pd.DataFrame([
        {"category": "Domain", "observed": observed_domain, "expected": round(expected_domain, 1),
         "enrichment": round(observed_domain / expected_domain, 2) if expected_domain else float("nan")},
        {"category": "Family", "observed": observed_family, "expected": round(expected_family, 1),
         "enrichment": round(observed_family / expected_family, 2) if expected_family else float("nan")},
    ])
    print("\n=== Restricted enrichment (corrected for mapping collisions) ===")
    print(summary.to_string(index=False))

    observed = [observed_domain, observed_family]
    expected = [expected_domain, expected_family]
    total_obs, total_exp = sum(observed), sum(expected)
    expected_scaled = [e * total_obs / total_exp for e in expected]
    chi2, pval = chisquare(observed, f_exp=expected_scaled)
    print(f"\nChi-square: chi2={chi2:.2f}  p-value={pval:.3e}")

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    summary.to_csv(os.path.join(OUTPUT_DIR, "enrichment_domain_vs_family_both.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    main()
