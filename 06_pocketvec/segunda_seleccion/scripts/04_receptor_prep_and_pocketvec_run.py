#!/usr/bin/env python3
"""
Receptor preparation and the final PocketVec docking run for the 3
finalist groups (25 species x group combinations, n_runs=25 per molecule).

Receptor prep: strips waters/ions/co-crystallized ligand and alternate
conformations from each PDB, keeping only the protein chain, then
converts to MOL2 with Gasteiger partial charges (OpenBabel). The pocket
centroid is written as a minimal .sd file for PocketVec's cavity-search
step.

Crystal vs. AlphaFold as receptor: for 5 of the 25 combinations, using
the species' own crystal structure as the receptor produced a completely
flat PocketVec descriptor (all 128 library molecules got the same score
-- the docking engine never placed a molecule in a valid pose near the
centroid). All 5 cases involved crystals that were partial/truncated
fragments of the full protein. Substituting the AlphaFold full-length
model as the receptor (keeping the crystal only as the source of the
binding-site position, transferred via US-align) resolved all 5 cases.
Those species are therefore run with the AlphaFold receptor instead of
the crystal, even though a crystal structure exists for them.
"""
import os
import subprocess
import time
from multiprocessing import Pool
from pathlib import Path

import pandas as pd
from Bio.PDB import PDBParser, PDBIO, Select

BASE_DIR   = Path("/home/ssneider/disco-big/ssneider-env/TDR_Targets7.1/PocketVec")
ETAPA3_DIR = BASE_DIR / "Nueva_estrategia_pocketvec" / "etapa3_pocketvec"

POCKETVEC_DIR  = BASE_DIR / "pocketvec"
POCKETVEC_MAIN = POCKETVEC_DIR / "pocketvec_main.py"
RDOCK_COMPILED_FUENTE = BASE_DIR / "rdock" / "rDock-main"

ALPHAFOLD_DIR  = ETAPA3_DIR / "alphafold_pdbs"
MOL2_DIR       = ETAPA3_DIR / "mol2_receptores"
CENTROIDES_DIR = ETAPA3_DIR / "sd_centroides"
POCKETVEC_OUT  = ETAPA3_DIR / "pocketvec_outputs"
GRUPOS_PICKLE  = ETAPA3_DIR / "GRUPOS_checkpoint.pkl"
RESULT_TSV     = ETAPA3_DIR / "resultado_corrida_DEFINITIVA_n25.tsv"

# The 5 species/group combinations where the crystal produced a flat
# descriptor; these use the AlphaFold receptor instead, keeping the
# crystal-derived centroid position
FLAT_DESCRIPTOR_FIX = {
    ("OG7_0001020", "hsap"): "Q13451",
    ("OG7_0006854", "hsap"): "Q02127",
    ("OG7_0006854", "lmaj"): "Q4QEW7",
    ("OG7_0000433", "atha"): "Q9LKB9",
    ("OG7_0000433", "cele"): "P02566",
}


def clean_crystal(pdb_in, chain, pdb_out):
    """Keeps only standard amino-acid residues of the given chain, main
    alternate conformation -- drops waters, ions, and the co-crystallized
    ligand along with everything else non-protein."""
    structure = PDBParser(QUIET=True).get_structure("st", pdb_in)

    class ProteinOnly(Select):
        def accept_residue(self, residue):
            if residue.get_parent().id != chain:
                return 0
            # residue.id[0] == " " means "standard amino acid residue";
            # anything else (H_FMN, H_ORO, W, ...) is a heteroatom
            return residue.id[0] == " "

        def accept_atom(self, atom):
            return atom.get_altloc() in (" ", "A")

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(pdb_out), ProteinOnly())


def convert_to_mol2(pdb_path, mol2_path):
    cmd = ["obabel", str(pdb_path), "-O", str(mol2_path), "--partialcharge", "gasteiger"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0 and Path(mol2_path).exists()


def create_pocket_centroid_sd(centroid, outfile):
    """Writes a minimal single-atom PDB at the centroid position, then
    converts it to .sd with OpenBabel -- PocketVec's expected centroid
    format."""
    x, y, z = (str(round(c, 3)) for c in centroid)
    ctr = " " * (8 - len(x)) + x + " " * (8 - len(y)) + y + " " * (8 - len(z)) + z
    text = "HEADER\nHETATM    1   C  CTR A   1    " + ctr + "  1.00  1.00           C\nEND"
    pdb_tmp = str(outfile).replace(".sd", ".pdb")
    with open(pdb_tmp, "w") as f:
        f.write(text)
    subprocess.run(["obabel", pdb_tmp, "-O", str(outfile)], capture_output=True, text=True)
    os.remove(pdb_tmp)


def prepare_receptors(GRUPOS):
    """Builds the final list of (og7, species, mol2_path, sd_path, outpath)
    tasks, choosing crystal or AlphaFold per the flat-descriptor fix
    above."""
    MOL2_DIR.mkdir(parents=True, exist_ok=True)
    CENTROIDES_DIR.mkdir(parents=True, exist_ok=True)
    POCKETVEC_OUT.mkdir(parents=True, exist_ok=True)

    tasks = []
    for og7, info in GRUPOS.items():
        for species, esp_info in info["especies"].items():
            use_alphafold_fix = (og7, species) in FLAT_DESCRIPTOR_FIX

            if use_alphafold_fix or esp_info["modo"] == "alphafold":
                pdb_source = esp_info.get("alphafold_pdb_path")
                chain_source = "A"
                suffix = "CHECK" if use_alphafold_fix else "FINAL"
            else:
                pdb_source = esp_info.get("pdb_path")
                chain_source = esp_info["chain"]
                suffix = "FINAL"

            if pdb_source is None:
                print(f"  WARNING: no source structure for {og7}|{species}, skipping")
                continue

            pdb_clean = MOL2_DIR / f"{og7}_{species}_clean_{suffix}.pdb"
            clean_crystal(pdb_source, chain_source, pdb_clean)
            mol2_path = MOL2_DIR / f"{og7}_{species}_{suffix}.mol2"
            if not convert_to_mol2(pdb_clean, mol2_path):
                print(f"  WARNING: MOL2 conversion failed for {og7}|{species}")
                continue

            sd_path = CENTROIDES_DIR / f"{og7}_{species}_{suffix}_centroide.sd"
            create_pocket_centroid_sd(esp_info["centroide"], sd_path)

            outpath = POCKETVEC_OUT / f"{og7}_{species}_DEFINITIVO"
            tasks.append((og7, species, str(mol2_path), str(sd_path), str(outpath)))

    return tasks


def run_one_species(args):
    og7, especie, mol2_path, sd_path, outpath = args
    outpath = Path(outpath)
    outpath.mkdir(exist_ok=True, parents=True)
    for f in ["cavity.grd", "cavity_log.log", "st_parameters.as"]:
        (outpath / f).unlink(missing_ok=True)

    cmd = [
        "conda", "run", "-n", "pocketvec_env", "python", str(POCKETVEC_MAIN),
        "-r", str(mol2_path), "-pc", str(sd_path), "-o", str(outpath),
        "-nr", "25", "-radius", "12.0", "-s", "42",
        "-rDock", str(RDOCK_COMPILED_FUENTE) + "/",
    ]
    t0 = time.time()
    result = subprocess.run(cmd, cwd=str(POCKETVEC_DIR), capture_output=True, text=True)
    elapsed = time.time() - t0
    fp_exists = (outpath / "PocketVec_fp.pkl").exists()
    print(f"[{og7}|{especie}] rc={result.returncode} fp={fp_exists} t={elapsed:.0f}s", flush=True)
    return {"og7": og7, "especie": especie, "return_code": result.returncode,
            "fp_generado": fp_exists, "tiempo_seg": round(elapsed, 1)}


def main():
    with open(GRUPOS_PICKLE, "rb") as f:
        GRUPOS = pickle.load(f)

    print("Preparing receptors and centroids...")
    tasks = prepare_receptors(GRUPOS)
    print(f"{len(tasks)} species/group combinations ready to dock, "
          f"{len(FLAT_DESCRIPTOR_FIX)} of them using the AlphaFold-receptor fix")

    print(f"\nRunning PocketVec, 8 in parallel, n_runs=25...")
    with Pool(processes=8) as pool:
        results = pool.map(run_one_species, tasks)

    df = pd.DataFrame(results)
    df.to_csv(RESULT_TSV, sep="\t", index=False)
    print(f"\nDone. Failed runs: {(df['return_code'] != 0).sum()} | "
          f"Missing PocketVec_fp.pkl: {(~df['fp_generado']).sum()}")


if __name__ == "__main__":
    import pickle
    main()
