"""
06_plot_confidence_distributions.py

Generates the following figures from the AlphaFold confidence score tables:

  1. Per-species histogram  : binned confidence averages for each species
  2. PFAM global histogram  : distribution of mean confidence across all PFAMs
  3. High-confidence bar    : % of PFAM annotations with avg score > 80, per species

All plots are saved as PNG files in the specified output directory.

Usage:
    python 06_plot_confidence_distributions.py \
        --all_species  all_species_confidence.tsv \
        --pfam_stats   PFAM_confidence_avg.tsv \
        --outdir       figures/
"""

import argparse
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


# ── 1. Per-species binned histogram ─────────────────────────────────────────

def plot_per_species_histogram(df: pd.DataFrame, outdir: str) -> None:
    bins   = [0, 40, 50, 60, 65, 70, 75, 80, 85, 90, 95, 100]
    labels = ["<40","40-50","50-60","60-65","65-70","70-75","75-80","80-85","85-90","90-95","95-100"]

    for species in df["especie"].unique():
        sub = df[df["especie"] == species].dropna(subset=["Confidence_Avg"])
        sub = sub.copy()
        sub["grupo"] = pd.cut(sub["Confidence_Avg"], bins=bins, labels=labels,
                              include_lowest=True, right=False)
        counts = sub["grupo"].value_counts().sort_index()

        plt.figure(figsize=(10, 6))
        counts.plot(kind="bar", color="orange", edgecolor="black")
        plt.title(f"AlphaFold confidence distribution — {species}")
        plt.xlabel("Confidence average (%)")
        plt.ylabel("Number of proteins")
        plt.xticks(rotation=45)
        plt.tight_layout()
        out_path = os.path.join(outdir, f"{species}_confidence_histogram.png")
        plt.savefig(out_path, dpi=300)
        plt.close()
        print(f"Saved: {out_path}")


# ── 2. PFAM-level global distribution ───────────────────────────────────────

def plot_pfam_global_distribution(pfam_stats: pd.DataFrame, outdir: str) -> None:
    plt.figure(figsize=(8, 6))
    sns.histplot(pfam_stats["Confidence_Avg"].dropna(), bins=30, kde=True, color="skyblue")
    plt.title("Confidence average distribution across all PFAMs")
    plt.xlabel("Confidence average")
    plt.ylabel("Number of PFAMs")
    plt.tight_layout()
    out_path = os.path.join(outdir, "PFAM_global_distribution.png")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Saved: {out_path}")


# ── 3. High-confidence PFAM percentage per species ───────────────────────────

def plot_high_confidence_by_species(df: pd.DataFrame, outdir: str, threshold: float = 80.0) -> None:
    scored = df.dropna(subset=["Confidence_Avg"])

    total = (
        scored.groupby("especie")["PFAM"]
        .count()
        .rename("Total")
        .reset_index()
    )
    above = (
        scored[scored["Confidence_Avg"] > threshold]
        .groupby("especie")["PFAM"]
        .count()
        .rename("Above_threshold")
        .reset_index()
    )
    summary = total.merge(above, on="especie", how="left").fillna(0)
    summary["Pct"] = (summary["Above_threshold"] / summary["Total"] * 100).round(1)
    summary = summary.sort_values("Pct", ascending=False)

    plt.figure(figsize=(10, 6))
    plt.bar(summary["especie"], summary["Pct"], color="darkorange")
    plt.xlabel("Species")
    plt.ylabel(f"PFAM with confidence avg > {threshold} (%)")
    plt.title(f"Proportion of high-confidence PFAM annotations per species (>{threshold}%)")
    plt.xticks(rotation=90, ha="center")
    plt.tight_layout()
    out_path = os.path.join(outdir, "species_high_confidence_pct.png")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Saved: {out_path}")

    print("\nSummary table:")
    print(summary.to_string(index=False))


# ── 4. Per-species smooth distribution (husl palette) ───────────────────────

def plot_per_species_smooth(df: pd.DataFrame, outdir: str) -> None:
    species_list = df["especie"].unique()
    palette = sns.color_palette("husl", len(species_list))

    subdir = os.path.join(outdir, "per_species_smooth")
    os.makedirs(subdir, exist_ok=True)

    for species, color in zip(species_list, palette):
        sub = df[(df["especie"] == species) & df["Confidence_Avg"].notna()]
        plt.figure(figsize=(8, 6))
        sns.histplot(sub["Confidence_Avg"], bins=30, kde=True, color=color, edgecolor="black")
        plt.title(f"Confidence distribution — {species}")
        plt.xlabel("Confidence average")
        plt.ylabel("Number of PFAMs")
        plt.xlim(0, 100)
        plt.tight_layout()
        out_path = os.path.join(subdir, f"{species}_smooth_distribution.png")
        plt.savefig(out_path, dpi=300)
        plt.close()

    print(f"Per-species smooth plots saved to: {subdir}")


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description="Generate AlphaFold confidence distribution plots.")
    parser.add_argument("--all_species", required=True, help="all_species_confidence.tsv")
    parser.add_argument("--pfam_stats",  required=True, help="PFAM_confidence_avg.tsv")
    parser.add_argument("--outdir",      required=True, help="Output directory for figures")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df_all   = pd.read_csv(args.all_species, sep="\t")
    df_pfam  = pd.read_csv(args.pfam_stats,  sep="\t")

    plot_per_species_histogram(df_all, args.outdir)
    plot_pfam_global_distribution(df_pfam, args.outdir)
    plot_high_confidence_by_species(df_all, args.outdir)
    plot_per_species_smooth(df_all, args.outdir)


if __name__ == "__main__":
    main()
