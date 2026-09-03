#!/usr/bin/env python3
"""
Compares PocketVec descriptors across species and OG7 groups for the
first-stage pilot round.

Computes pairwise cosine similarity between all descriptors, then breaks
it down into intra-OG7 (same group, different species -- expected to be
high if the pocket is conserved) and inter-OG7 (same species, different
group -- expected to be low, since the two proteins are unrelated).
"""
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import cosine
import matplotlib.pyplot as plt
import seaborn as sns

OUTPUT_DIR = "/big/lab/ssneider/ssneider-env/TDR_Targets7.1/PocketVec/pLDDT_candidatos/pocketvec_outputs"
BASE_DIR   = "/big/lab/ssneider/ssneider-env/TDR_Targets7.1/PocketVec/pLDDT_candidatos"

OG7_GROUPS = ["OG7_0006581", "OG7_0006003"]
OG7_LABELS = {
    "OG7_0006581": "OG7_0006581\n(Decanoic acid, mtub)",
    "OG7_0006003": "OG7_0006003\n(Arabinofuranose-P, ecol)",
}


def load_descriptors():
    descriptors = {}
    for folder in sorted(Path(OUTPUT_DIR).iterdir()):
        pkl_path = folder / "PocketVec_fp.pkl"
        if not pkl_path.exists():
            print(f"  WARNING: {pkl_path} not found")
            continue
        name = folder.name  # e.g. OG7_0006581_hsap
        og7, species = None, name
        for og7_id in OG7_GROUPS:
            if name.startswith(og7_id):
                species = name[len(og7_id) + 1:]
                og7 = og7_id
                break
        vec = pickle.load(open(pkl_path, "rb"))
        descriptors[(og7, species)] = np.array(vec)
        print(f"  Loaded: {og7} | {species} -> vector dim={len(vec)}")
    print(f"\nTotal descriptors loaded: {len(descriptors)}")
    return descriptors


def cosine_similarity_matrix(descriptors):
    keys = sorted(descriptors.keys())
    labels = [f"{sp}\n({og.replace('OG7_000', 'OG7...')})" for og, sp in keys]
    n = len(keys)

    sim_matrix = np.zeros((n, n))
    for i, k1 in enumerate(keys):
        for j, k2 in enumerate(keys):
            sim_matrix[i, j] = 1 - cosine(descriptors[k1], descriptors[k2])

    df_sim = pd.DataFrame(sim_matrix, index=labels, columns=labels)
    df_sim.to_csv(f"{BASE_DIR}/pocketvec_similitud_coseno.tsv", sep="\t")
    print("Saved: pocketvec_similitud_coseno.tsv")
    print(df_sim.round(3).to_string())
    return keys, labels, df_sim


def plot_heatmap(keys, df_sim):
    fig, ax = plt.subplots(figsize=(13, 11))
    n_group1 = sum(1 for og, _ in keys if og == OG7_GROUPS[0])

    sns.heatmap(
        df_sim, annot=True, fmt=".2f", cmap="YlOrRd", vmin=0, vmax=1, ax=ax,
        linewidths=0.5, linecolor="white",
        cbar_kws={"label": "Cosine similarity (1=identical, 0=no similarity)"},
        annot_kws={"size": 8},
    )
    ax.axhline(n_group1, color="black", linewidth=2.5)
    ax.axvline(n_group1, color="black", linewidth=2.5)
    ax.set_title(
        "PocketVec pocket similarity between species and OG7 groups\n"
        "Black line separates OG7_0006581 (top/left) from OG7_0006003 (bottom/right)",
        fontweight="bold", fontsize=11,
    )
    ax.tick_params(axis="x", rotation=45, labelsize=8)
    ax.tick_params(axis="y", rotation=0, labelsize=8)
    plt.tight_layout()
    plt.savefig(f"{BASE_DIR}/pocketvec_heatmap_general.png", dpi=150, bbox_inches="tight")
    plt.show()
    print("Saved: pocketvec_heatmap_general.png")


def intra_og7_similarity(descriptors, keys):
    print("\n── INTRA-OG7 SIMILARITY (same group, across species) ──────────────────")
    all_sims = {}
    for og7 in OG7_GROUPS:
        keys_og7 = [(og, sp) for og, sp in keys if og == og7]
        species = [sp for _, sp in keys_og7]
        n_sp = len(keys_og7)
        print(f"\n{og7} ({n_sp} species):")

        sims = []
        for i in range(n_sp):
            for j in range(i + 1, n_sp):
                s = 1 - cosine(descriptors[keys_og7[i]], descriptors[keys_og7[j]])
                sims.append(s)
                print(f"  {species[i]} vs {species[j]}: {s:.3f}")

        print(f"  Intra-OG7 mean: {np.mean(sims):.3f} ± {np.std(sims):.3f}")
        print(f"  Min: {np.min(sims):.3f} | Max: {np.max(sims):.3f}")
        all_sims[og7] = sims
    return all_sims


def inter_og7_similarity(descriptors, keys):
    print("\n── INTER-OG7 SIMILARITY (same species, different groups) ──────────────")
    common_species = (set(sp for og, sp in keys if og == OG7_GROUPS[0]) &
                       set(sp for og, sp in keys if og == OG7_GROUPS[1]))

    sims_inter = []
    for sp in sorted(common_species):
        k1, k2 = (OG7_GROUPS[0], sp), (OG7_GROUPS[1], sp)
        if k1 in descriptors and k2 in descriptors:
            s = 1 - cosine(descriptors[k1], descriptors[k2])
            sims_inter.append(s)
            print(f"  {sp}: {s:.3f}")

    print(f"\nInter-OG7 mean (same species): {np.mean(sims_inter):.3f} ± {np.std(sims_inter):.3f}")
    print("(If this is low, the two groups have chemically distinct pockets -- expected)")
    return sims_inter


def plot_intra_vs_inter(intra_sims, sims_inter):
    fig, ax = plt.subplots(figsize=(9, 5))
    data_box = [intra_sims[OG7_GROUPS[0]], intra_sims[OG7_GROUPS[1]], sims_inter]
    labels_box = [
        f"Intra-OG7\n{OG7_GROUPS[0]}\n(Decanoic acid)",
        f"Intra-OG7\n{OG7_GROUPS[1]}\n(Arabinofuranose-P)",
        "Inter-OG7\n(same species,\ndifferent group)",
    ]
    bp = ax.boxplot(data_box, labels=labels_box, patch_artist=True,
                     medianprops=dict(color="black", linewidth=2))
    colors = ["#2E75B6", "#70AD47", "#E74C3C"]
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax.set_ylabel("Cosine similarity", fontsize=11)
    ax.set_title(
        "Pocket similarity distribution by category\n"
        "Intra-OG7 = between species of the same group | Inter-OG7 = same species, different groups",
        fontweight="bold", fontsize=10,
    )
    ax.set_ylim(0, 1.05)
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{BASE_DIR}/pocketvec_boxplot_similitud.png", dpi=150, bbox_inches="tight")
    plt.show()
    print("Saved: pocketvec_boxplot_similitud.png")


def main():
    descriptors = load_descriptors()
    keys, labels, df_sim = cosine_similarity_matrix(descriptors)
    plot_heatmap(keys, df_sim)
    intra_sims = intra_og7_similarity(descriptors, keys)
    sims_inter = inter_og7_similarity(descriptors, keys)
    plot_intra_vs_inter(intra_sims, sims_inter)


if __name__ == "__main__":
    main()
