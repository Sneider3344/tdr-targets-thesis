"""
01_match_interpro_to_alphafold.py

Joins InterProScan PFAM annotation files with AlphaFold confidence JSON
files for a given species. For each annotated protein, verifies that the
corresponding AlphaFold JSON exists and writes a joined table with columns:
especie, ID, PFAM, AlphaFold_ID, Inicio, Fin.

Proteins without a matching JSON are skipped and logged as warnings.

Usage:
    python 01_match_interpro_to_alphafold.py \
        --interpro_dir  <folder_with_txt_files> \
        --alphafold_dir <folder_with_json_files> \
        --output        <joined_table.txt> \
        --species       atha
"""

import argparse
import os


def build_joined_table(interpro_dir: str, alphafold_dir: str, output_path: str, species: str) -> None:
    missing = 0
    written = 0

    with open(output_path, "w") as out_f:
        out_f.write("especie\tID\tPFAM\tAlphaFold_ID\tInicio\tFin\n")

        for file_name in os.listdir(interpro_dir):
            if not file_name.endswith(".txt"):
                continue

            with open(os.path.join(interpro_dir, file_name)) as f:
                for line in f:
                    if not line.strip():
                        continue

                    parts = line.strip().split()
                    especie_id = parts[0]
                    pfam       = parts[2]
                    inicio     = parts[4]
                    fin        = parts[5]

                    # Expect format 'species|UniprotID'
                    especie, uniprot_id = especie_id.split("|")
                    alphafold_filename  = f"AF-{uniprot_id}-F1-confidence_v4.json"
                    json_path = os.path.join(alphafold_dir, alphafold_filename)

                    if not os.path.exists(json_path):
                        print(f"[WARN] JSON not found: {alphafold_filename}")
                        missing += 1
                        continue

                    out_f.write(f"{especie}\t{uniprot_id}\t{pfam}\t{alphafold_filename}\t{inicio}\t{fin}\n")
                    written += 1

    print(f"Rows written     : {written}")
    print(f"Missing JSONs    : {missing}")
    print(f"Output           : {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Match InterProScan annotations to AlphaFold JSONs.")
    parser.add_argument("--interpro_dir",  required=True)
    parser.add_argument("--alphafold_dir", required=True)
    parser.add_argument("--output",        required=True)
    parser.add_argument("--species",       required=True)
    args = parser.parse_args()

    build_joined_table(args.interpro_dir, args.alphafold_dir, args.output, args.species)


if __name__ == "__main__":
    main()
