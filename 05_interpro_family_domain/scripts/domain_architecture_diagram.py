#!/usr/bin/env python3
"""
Generates schematic 2D domain architecture diagrams (Pfam-graphics style)
to illustrate specific Family vs Domain cases, from
pillar2b_family_domain_overlap.tsv.

Jupyter-friendly version: paste the whole cell and run. To switch cases,
edit the "EDIT THESE VALUES" block below and re-run the cell (no need to
restart the kernel).
"""
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

OVERLAP_TSV = "/home/ssneider/disco-big/ssneider-env/TDR_Targets7.1/interpro/nuevo_analisis_tesis/pillar2b_family_domain_overlap.tsv"

FAMILY_COLOR = "#4C72B0"  # blue
DOMAIN_COLOR = "#DD8452"  # orange
BACKBONE_COLOR = "#B0B0B0"


def load_cases():
    return pd.read_csv(OVERLAP_TSV, sep="\t")


def list_cases(df, overlap_filter=None, species=None, n=15):
    d = df.copy()
    if overlap_filter == "no_overlap":
        d = d[d["overlap_len"] == 0]
    elif overlap_filter == "nested":
        d = d[d["pct_domain_covered"] >= 99]
    elif overlap_filter == "partial":
        d = d[(d["overlap_len"] > 0) & (d["pct_domain_covered"] < 99)]

    if species:
        d = d[d["species"] == species]

    d = d.sort_values("protein_id")
    cols_to_show = ["species", "protein_id"]
    if "uniprot_id" in d.columns:
        cols_to_show.append("uniprot_id")
    cols_to_show += ["family_acc", "family_start", "family_end",
                      "domain_acc", "domain_start", "domain_end", "pct_domain_covered"]
    print(d[cols_to_show].head(n).to_string(index=False))
    print(f"\n({len(d)} total cases matching this filter, showing {min(n, len(d))})")


def draw_architecture(df, species, protein_id, protein_length=None, out_path="domain_architecture.png",
                       display_name=None):
    sub = df[(df["species"] == species) & (df["protein_id"] == protein_id)]
    if sub.empty:
        raise ValueError(f"No rows found for {species}/{protein_id} in {OVERLAP_TSV}")

    fam_regions = sub[["family_acc", "family_start", "family_end"]].drop_duplicates()
    dom_regions = sub[["domain_acc", "domain_start", "domain_end"]].drop_duplicates()

    if protein_length is None:
        protein_length = int(max(sub["family_end"].max(), sub["domain_end"].max()) * 1.05)
        length_note = ""
    else:
        length_note = ""

    fig, ax = plt.subplots(figsize=(10, 2.8))

    ax.add_patch(mpatches.FancyBboxPatch(
        (0, 0.42), protein_length, 0.16,
        boxstyle="round,pad=0,rounding_size=2",
        linewidth=1, edgecolor="black", facecolor=BACKBONE_COLOR, zorder=1))

    box_height = 0.34
    fam_y = 0.55
    dom_y = 0.05

    legend_handles = []
    seen_fam, seen_dom = set(), set()

    for _, row in fam_regions.iterrows():
        s, e, acc = row["family_start"], row["family_end"], row["family_acc"]
        ax.add_patch(mpatches.FancyBboxPatch(
            (s, fam_y), e - s, box_height,
            boxstyle="round,pad=0,rounding_size=1.5",
            linewidth=1.2, edgecolor="black", facecolor=FAMILY_COLOR, zorder=2))
        ax.text((s + e) / 2, fam_y + box_height / 2, acc,
                 ha="center", va="center", fontsize=8, color="white", fontweight="bold")
        if acc not in seen_fam:
            legend_handles.append(mpatches.Patch(color=FAMILY_COLOR, label=f"InterPro Family: {acc}"))
            seen_fam.add(acc)

    for _, row in dom_regions.iterrows():
        s, e, acc = row["domain_start"], row["domain_end"], row["domain_acc"]
        ax.add_patch(mpatches.FancyBboxPatch(
            (s, dom_y), e - s, box_height,
            boxstyle="round,pad=0,rounding_size=1.5",
            linewidth=1.2, edgecolor="black", facecolor=DOMAIN_COLOR, zorder=2))
        ax.text((s + e) / 2, dom_y + box_height / 2, acc,
                 ha="center", va="center", fontsize=8, color="white", fontweight="bold")
        if acc not in seen_dom:
            legend_handles.append(mpatches.Patch(color=DOMAIN_COLOR, label=f"InterPro Domain: {acc}"))
            seen_dom.add(acc)

    step = max(1, protein_length // 10)
    for pos in range(0, protein_length + 1, step):
        ax.plot([pos, pos], [0.40, 0.42], color="black", linewidth=0.8)
        ax.text(pos, 0.36, str(pos), ha="center", va="top", fontsize=6.5)

    ax.set_xlim(-protein_length * 0.02, protein_length * 1.02)
    ax.set_ylim(-0.05, 1.0)
    ax.axis("off")
    uniprot_id = sub["uniprot_id"].iloc[0] if "uniprot_id" in sub.columns else protein_id
    title_id = display_name if display_name else (
        f"{uniprot_id} ({protein_id})" if uniprot_id != protein_id else protein_id
    )
    ax.set_title(f"{species} | {title_id}  (length: {protein_length} aa{length_note})", fontsize=10)
    ax.legend(handles=legend_handles, loc="upper center", bbox_to_anchor=(0.5, -0.08),
              ncol=2, frameon=False, fontsize=8)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"Saved: {out_path}")
    plt.show()


# ══════════════════════════════════════════════════════════════════════
# EDIT THESE VALUES AND RE-RUN THE CELL (no need for %run)
# ══════════════════════════════════════════════════════════════════════
MODE = "plot"   # "list" to explore cases, "plot" to generate a figure

# --- parameters for MODE = "list" ---
LIST_OVERLAP_FILTER = "no_overlap"   # "no_overlap" | "nested" | "partial" | None
LIST_SPECIES = "tgon"                # species code, or None for all

# --- parameters for MODE = "plot" ---
PLOT_SPECIES = "tgon"
PLOT_PROTEIN = "tgon|TGME49_326800"
PLOT_LENGTH = None                          # true length in aa, or None to approximate
PLOT_DISPLAY_NAME = "S8EU04"                # <-- SET the name/ID you want in the title
                                             #     e.g. "A0A068WHV0 (HSP70)" -- leave None to use the automatic one
PLOT_OUT = "domain_architecture_tgon.png"
# ══════════════════════════════════════════════════════════════════════

df = load_cases()

if MODE == "list":
    list_cases(df, overlap_filter=LIST_OVERLAP_FILTER, species=LIST_SPECIES)
elif MODE == "plot":
    draw_architecture(df, PLOT_SPECIES, PLOT_PROTEIN, PLOT_LENGTH, PLOT_OUT, display_name=PLOT_DISPLAY_NAME)
