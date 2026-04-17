"""
04_plot_genes_per_orthogroup.py

Generates a bar chart showing how many orthogroups contain 1, 2, 3, ...
genes (i.e. the distribution of group sizes), for a single species/version CSV.

Usage:
    python 04_plot_genes_per_orthogroup.py \
        --input  <species_orthogroups.csv> \
        --title  "S. cerevisiae - OrthoMCL v7" \
        --output <plot.png>
"""

import argparse
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def plot_distribution(input_path: str, title: str, output_path: str) -> None:
    df = pd.read_csv(input_path)

    # Count genes per orthogroup, then count how many groups have each size
    group_sizes = df["orthology_group"].value_counts()
    size_distribution = group_sizes.value_counts().sort_index()

    plt.figure(figsize=(12, 6))
    sns.barplot(x=size_distribution.index, y=size_distribution.values, color="skyblue")
    plt.title(title, fontsize=25)
    plt.xlabel("Number of genes per orthogroup")
    plt.ylabel("Number of orthogroups")
    plt.xticks(rotation=45)
    plt.ylim(0, size_distribution.max() * 1.1)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Plot saved: {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot gene count distribution per orthogroup.")
    parser.add_argument("--input",  required=True, help="CSV with orthology_group column")
    parser.add_argument("--title",  default="Gene distribution per orthogroup")
    parser.add_argument("--output", required=True, help="Output PNG path")
    args = parser.parse_args()

    plot_distribution(args.input, args.title, args.output)


if __name__ == "__main__":
    main()
