"""
03_fasta_headers_to_csv.py

Parses FASTA header lines from OrthoMCL v6 and v7 outputs and writes a
structured CSV with columns: gene_id, protein, orthology_group.

Header formats differ between versions:
  v6: >tbrt|Tb11.v5.0835 | OG6_105680 | unknown
  v7: >tbrt|Tb927.9.9410 | organism=Trypanosoma brucei | unknown | OG7_0004560.

Usage:
    python 03_fasta_headers_to_csv.py \
        --input  <species_orthogroups.fasta> \
        --version V6|V7 \
        --output <output.csv>
"""

import argparse
import csv


def parse_v6_header(line: str) -> list:
    parts = [p.strip() for p in line.strip().split("|")]
    gene_id = parts[0] + "|" + parts[1]
    orthology_group = parts[2]
    protein = parts[3]
    return [gene_id, protein, orthology_group]


def parse_v7_header(line: str) -> list:
    parts = [p.strip() for p in line.strip().split("|")]
    gene_id = parts[0] + "|" + parts[1]
    protein = parts[3]
    orthology_group = parts[4].rstrip(".")
    return [gene_id, protein, orthology_group]


def convert(input_path: str, version: str, output_path: str) -> int:
    parser_fn = parse_v7_header if version == "V7" else parse_v6_header
    count = 0

    with open(input_path, encoding="utf-8") as infile, \
         open(output_path, "w", newline="", encoding="utf-8") as outfile:
        writer = csv.writer(outfile)
        writer.writerow(["gene_id", "protein", "orthology_group"])
        for line in infile:
            if line.startswith(">"):
                writer.writerow(parser_fn(line))
                count += 1

    return count


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert FASTA headers to structured CSV.")
    parser.add_argument("--input",   required=True)
    parser.add_argument("--version", required=True, choices=["V6", "V7"])
    parser.add_argument("--output",  required=True)
    args = parser.parse_args()

    n = convert(args.input, args.version, args.output)
    print(f"Sequences parsed : {n}")
    print(f"Output written   : {args.output}")


if __name__ == "__main__":
    main()
