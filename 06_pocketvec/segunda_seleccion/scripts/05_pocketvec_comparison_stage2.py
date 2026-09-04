#!/usr/bin/env python3
"""
Final intra-/inter-OG7 comparison of PocketVec descriptors for the 3
finalist groups.

Reproduces the thesis's Table 4 (per-group pair statistics against the
0.83 similarity threshold) and Figure 21 (boxplot of intra-OG7 vs.
inter-OG7 cosine similarity). Inter-OG7 comparisons are only possible for
species present in more than one of the 3 groups.
"""
import pickle
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import cosine

ETAPA3_DIR = Path("/home/ssneider/disco-big/ssneider-env/TDR_Targets7.1/PocketVec/Nueva_estrategia_pocketvec/etapa3_pocketvec")
OUTPUT_DIR = ETAPA3_DIR / "pocketvec_outputs"
GRUPOS_PICKLE = ETAPA3_DIR / "GRUPOS_checkpoint.pkl"

OG7_GROUPS = ["OG7_0000433", "OG7_0001020", "OG7_0006854"]
OG7_LABELS = {
    "OG7_0000433": "OG7_0000433\n(ADP)",
    "OG7_0001020": "OG7_0001020\n(RAP)",
    "OG7_0006854": "OG7_0006854\n(FMN)",
}
SIMILARITY_THRESHOLD = 0.83


def load_descriptors(GRUPOS):
    descriptors = {}
    for og7 in OG7_GROUPS:
        for species in GRUPOS[og7]["especies"]:
            pkl_path = OUTPUT_DIR / f"{og7}_{species}_DEFINITIVO" / "PocketVec_fp.pkl"
            if not pkl_path.exists():
                print(f"  WARNING: missing {pkl_path}")
                continue
            descriptors[(og7, species)] = np.array(pickle.load(open(pkl_path, "rb")))
    print(f"Total descriptors loaded: {len(descriptors)}")
    return descriptors


def intra_og7_table(descriptors, GRUPOS):
    """Table 4: per-group pairwise similarity statistics."""
    rows = []
    intra_by_group = {}
    for og7 in OG7_GROUPS:
        species = [sp for sp in GRUPOS[og7]["especies"] if (og7, sp) in descriptors]
        sims = []
        for a, b in combinations(species, 2):
            s = 1 - cosine(descriptors[(og7, a)], descriptors[(og7, b)])
            sims.append((a, b, s))
        intra_by_group[og7] = sims

        values = [s for _, _, s in sims]
        above = sum(1 for v in values if v > SIMILARITY_THRESHOLD)
        best = max(sims, key=lambda x: x[2])
        worst = min(sims, key=lambda x: x[2])
        rows.append({
            "grupo": og7, "ligando": GRUPOS[og7]["ligando"], "especies": len(species),
            "pares": len(sims), "promedio": round(np.mean(values), 3),
            f"pares_>{SIMILARITY_THRESHOLD}": f"{above}/{len(sims)} ({100*above/len(sims):.0f}%)",
            "ejemplo_mas_alto": f"{best[0]}-{best[1]} ({best[2]:.3f})",
            "ejemplo_mas_bajo": f"{worst[0]}-{worst[1]} ({worst[2]:.3f})",
        })
    df = pd.DataFrame(rows)
    print("\n=== Table 4: intra-OG7 statistics ===")
    print(df.to_string(index=False))
    return df, intra_by_group


def inter_og7_pairs(descriptors, GRUPOS):
    """Inter-OG7: same species, compared across pairs of distinct groups."""
    rows = []
    for og7_a, og7_b in combinations(OG7_GROUPS, 2):
        species_a = set(GRUPOS[og7_a]["especies"])
        species_b = set(GRUPOS[og7_b]["especies"])
        common = species_a & species_b
        for sp in sorted(common):
            if (og7_a, sp) not in descriptors or (og7_b, sp) not in descriptors:
                continue
            s = 1 - cosine(descriptors[(og7_a, sp)], descriptors[(og7_b, sp)])
            rows.append({"og7_a": og7_a, "og7_b": og7_b, "especie": sp, "similitud": round(s, 3)})
    df_inter = pd.DataFrame(rows)
    n_above = (df_inter["similitud"] > SIMILARITY_THRESHOLD).sum()
    print(f"\n=== Inter-OG7 comparisons (same species, unrelated proteins) ===")
    print(df_inter.to_string(index=False))
    print(f"\n{n_above} of {len(df_inter)} pairs ({100*n_above/len(df_inter):.0f}%) "
          f"above the {SIMILARITY_THRESHOLD} threshold")
    return df_inter


def plot_intra_vs_inter(intra_by_group, df_inter):
    """Figure 21: boxplot per group, all 3 groups plus the pooled
    inter-OG7 distribution."""
    fig, ax = plt.subplots(figsize=(9, 5))
    data = [[s for _, _, s in intra_by_group[og7]] for og7 in OG7_GROUPS]
    data.append(list(df_inter["similitud"]))
    labels = [OG7_LABELS[og7] for og7 in OG7_GROUPS] + ["Inter-OG7\n(unrelated proteins)"]

    bp = ax.boxplot(data, labels=labels, patch_artist=True, medianprops=dict(color="black", linewidth=2))
    colors = ["#2E75B6", "#70AD47", "#9E5FA3", "#E74C3C"]
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax.axhline(SIMILARITY_THRESHOLD, color="red", linestyle="--", linewidth=1,
               label=f"Threshold ({SIMILARITY_THRESHOLD})")
    ax.set_ylabel("PocketVec cosine similarity", fontsize=11)
    ax.set_title("Cosine similarity distribution for three binding-pocket groups\n"
                 "(ADP, RAP, FMN), intra- and inter-OG7", fontweight="bold", fontsize=11)
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(ETAPA3_DIR / "figura21_intra_vs_inter_og7.png", dpi=150, bbox_inches="tight")
    print(f"\nSaved: figura21_intra_vs_inter_og7.png")


def main():
    with open(GRUPOS_PICKLE, "rb") as f:
        GRUPOS = pickle.load(f)

    descriptors = load_descriptors(GRUPOS)
    df_table4, intra_by_group = intra_og7_table(descriptors, GRUPOS)
    df_inter = inter_og7_pairs(descriptors, GRUPOS)
    plot_intra_vs_inter(intra_by_group, df_inter)

    df_table4.to_csv(ETAPA3_DIR / "tabla4_estadisticas_por_grupo.tsv", sep="\t", index=False)
    df_inter.to_csv(ETAPA3_DIR / "pares_inter_og7.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
