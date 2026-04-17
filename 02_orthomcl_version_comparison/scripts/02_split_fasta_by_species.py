"""
02_split_fasta_by_species.py

Splits a multi-species FASTA (headers only or full sequences) into one
file per species, based on the header prefix (e.g. '>scer', '>tbrt').

Usage:
    python 02_split_fasta_by_species.py \
        --input   <filtered.fasta> \
        --species scer tbrt \
        --outdir  <output_directory>
"""

import argparse
import os


def split_fasta(input_path: str, species: list[str], outdir: str) -> dict:
    os.makedirs(outdir, exist_ok=True)
    handles = {
        s: open(os.path.join(outdir, f"{s}_orthogroups.fasta"), "w", encoding="utf-8")
        for s in species
    }
    counts = {s: 0 for s in species}
    current = None

    with open(input_path, encoding="utf-8") as infile:
        for line in infile:
            if line.startswith(">"):
                current = None
                for s in species:
                    if line.startswith(f">{s}"):
                        current = s
                        counts[s] += 1
                        break
            if current:
                handles[current].write(line)

    for fh in handles.values():
        fh.close()

    return counts


def main() -> None:
    parser = argparse.ArgumentParser(description="Split FASTA by species prefix.")
    parser.add_argument("--input",   required=True, help="Filtered multi-species FASTA")
    parser.add_argument("--species", nargs="+", default=["scer", "tbrt"])
    parser.add_argument("--outdir",  required=True, help="Directory for per-species output files")
    args = parser.parse_args()

    counts = split_fasta(args.input, args.species, args.outdir)
    for species, n in counts.items():
        print(f"{species}: {n} sequences written")


if __name__ == "__main__":
    main()
