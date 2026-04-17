"""
02_remap_ids_to_uniprot.py

Replaces species-specific protein IDs in an InterProScan annotation file
with their corresponding UniProt IDs, using a two-column mapping file
(From -> Entry) downloaded from the UniProt ID mapping tool.

Rows whose original ID has no mapping entry are left unchanged.

Usage:
    python 02_remap_ids_to_uniprot.py \
        --interpro   <species_iprscan_data.txt> \
        --mapping    <idmapping_species.tsv> \
        --output     <species_iprscan_curated.txt>
"""

import argparse
import csv


def load_id_map(mapping_path: str) -> dict:
    """Return a dict {From: Entry} from a UniProt mapping TSV."""
    id_map = {}
    with open(mapping_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            id_map[row["From"].strip()] = row["Entry"].strip()
    return id_map


def remap_file(interpro_path: str, id_map: dict, output_path: str) -> tuple[int, int]:
    remapped = 0
    unchanged = 0

    with open(interpro_path) as in_f, open(output_path, "w") as out_f:
        for line in in_f:
            if not line.strip():
                continue

            parts = line.strip().split()
            especie, old_id = parts[0].split("|")
            new_id = id_map.get(old_id, old_id)

            if new_id != old_id:
                remapped += 1
            else:
                unchanged += 1

            parts[0] = f"{especie}|{new_id}"
            out_f.write("\t".join(parts) + "\n")

    return remapped, unchanged


def main() -> None:
    parser = argparse.ArgumentParser(description="Remap protein IDs to UniProt in InterProScan output.")
    parser.add_argument("--interpro", required=True)
    parser.add_argument("--mapping",  required=True)
    parser.add_argument("--output",   required=True)
    args = parser.parse_args()

    id_map = load_id_map(args.mapping)
    remapped, unchanged = remap_file(args.interpro, id_map, args.output)

    print(f"IDs remapped     : {remapped}")
    print(f"IDs unchanged    : {unchanged}")
    print(f"Output           : {args.output}")


if __name__ == "__main__":
    main()
