"""
01_filter_fasta_by_species.py

Extracts sequences for target species (scer, tbrt) from a compressed
full-database OrthoMCL FASTA file. Reads the gzipped input line by line
to avoid loading the entire file into memory.

Usage:
    python 01_filter_fasta_by_species.py \
        --input  <orthomcl_full.fasta.gz> \
        --output <filtered_scer_tbrt.fasta> \
        --species scer tbrt
"""

import argparse
import gzip


def filter_fasta(input_path: str, output_path: str, species: list[str]) -> int:
    """Write only sequences whose header starts with one of the target species prefixes."""
    prefixes = tuple(f">{s}" for s in species)
    kept = 0
    write = False

    with gzip.open(input_path, "rt", encoding="utf-8") as infile, \
         open(output_path, "w", encoding="utf-8") as outfile:
        for line in infile:
            if line.startswith(">"):
                write = line.startswith(prefixes)
                if write:
                    kept += 1
            if write:
                outfile.write(line)

    return kept


def main() -> None:
    parser = argparse.ArgumentParser(description="Filter OrthoMCL FASTA by species prefix.")
    parser.add_argument("--input",   required=True, help="Compressed input FASTA (.fasta.gz)")
    parser.add_argument("--output",  required=True, help="Filtered output FASTA")
    parser.add_argument("--species", nargs="+", default=["scer", "tbrt"],
                        help="Species prefixes to keep (default: scer tbrt)")
    args = parser.parse_args()

    kept = filter_fasta(args.input, args.output, args.species)
    print(f"Sequences kept : {kept}")
    print(f"Output written : {args.output}")


if __name__ == "__main__":
    main()
