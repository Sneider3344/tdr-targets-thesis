#!/usr/bin/env python3
"""
Protein-FMN interaction analysis for OG7_0006854, comparing 3
crystallized species (H. sapiens 4oqv/8yhr, L. major 4xq6, P. falciparum
7l01; M. tuberculosis 3tq0 added for the visual comparison only).

Two independent, non-PocketVec ways of comparing "how similar the FMN
binding site looks" across species, used to cross-check the PocketVec
result for this group (see 05_pocketvec_comparison_stage2.py):

1. Raw distance contacts: side-chain atoms (excluding backbone and
   hydrogens) within 4.0 A of any FMN atom.
2. PLIP-classified contacts: the same structures run through PLIP, which
   labels each protein-ligand contact as a hydrogen bond, a hydrophobic
   interaction, or unclassified proximity.

Because residue numbering differs between the crystal structures, contact
positions from each species are translated into a common frame using
H. sapiens as the anchor, via pairwise BLOSUM62 sequence alignment.

ChimeraX renders one visual scene per species (protein + FMN + hydrogen
bonds, via the `hbonds` command) for the figure in the thesis -- that is
a visualization step, not part of the analysis itself, and isn't
reproduced here.
"""
import re
import subprocess
from pathlib import Path

import gzip
import pandas as pd
from Bio import pairwise2
from Bio.Align import substitution_matrices
from Bio.PDB import PDBParser, NeighborSearch

BASE_DIR = Path("/home/ssneider/disco-big/ssneider-env/TDR_Targets7.1/PocketVec")
OUT_DIR = BASE_DIR / "Nueva_estrategia_pocketvec" / "validacion_alineamientos"

STRUCTURES = {
    "hsap": {"pdb_id": "8yhr", "chain": "A"},
    "lmaj": {"pdb_id": "4xq6", "chain": "A"},
    "pfal": {"pdb_id": "7l01", "chain": "A"},
    "mtub": {"pdb_id": "3tq0", "chain": "A"},  # visual comparison only
}
ANCHOR_SPECIES = "hsap"
CONTACT_DISTANCE_A = 4.0

COLS_BIOLIP = [
    "pdb_id", "chain", "resolution", "binding_site",
    "ligand_id", "ligand_chain", "ligand_serial",
    "binding_residues_pdb", "binding_residues_renum",
    "catalytic_pdb", "catalytic_renum",
    "ec_number", "go_terms",
    "affinity_manual", "affinity_moad", "affinity_pdbbind", "affinity_bindingdb",
    "uniprot_id", "pubmed_id", "ligand_seqnum", "receptor_seq",
]

BACKBONE_ATOMS = {"N", "CA", "C", "O"}


def load_biolip_row(pdb_id, chain):
    with gzip.open(BASE_DIR / "analisis_biolip" / "BioLiP.txt.gz", "rt", encoding="utf-8", errors="replace") as f:
        df = pd.read_csv(f, sep="\t", header=None, names=COLS_BIOLIP, low_memory=False)
    df["pdb_id"] = df["pdb_id"].astype(str).str.strip().str.lower()
    row = df[(df["pdb_id"] == pdb_id) & (df["chain"] == chain)]
    return row.iloc[0] if len(row) else None


def download_pdb(pdb_id, outdir):
    import requests
    outpath = Path(outdir) / f"{pdb_id.lower()}.pdb"
    if outpath.exists():
        return outpath
    r = requests.get(f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb", timeout=30)
    if r.status_code == 200:
        outpath.write_text(r.text)
        return outpath
    return None


def raw_contacts_to_fmn(pdb_path, chain, distance=CONTACT_DISTANCE_A):
    """Side-chain (non-backbone, non-hydrogen) atoms within `distance` A
    of any FMN atom."""
    structure = PDBParser(QUIET=True).get_structure("st", pdb_path)
    model = structure[0]

    fmn_atoms = [atom for res in model[chain] if res.get_resname() == "FMN" for atom in res]
    if not fmn_atoms:
        return set()

    protein_atoms = [
        atom for res in model[chain]
        if res.id[0] == " "
        for atom in res
        if atom.get_name() not in BACKBONE_ATOMS and atom.element != "H"
    ]
    ns = NeighborSearch(protein_atoms)

    contact_residues = set()
    for fmn_atom in fmn_atoms:
        for atom in ns.search(fmn_atom.coord, distance):
            contact_residues.add(atom.get_parent().id[1])
    return contact_residues


def run_plip(pdb_path):
    """Runs PLIP on the structure; returns its raw XML report path for
    downstream parsing of bond-type classification."""
    result = subprocess.run(["plip", "-f", str(pdb_path), "-x", "-o", str(pdb_path.parent)],
                             capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  PLIP error on {pdb_path.name}: {result.stderr[:200]}")
    return pdb_path.parent / "report.xml"


def get_sequence(pdb_path, chain):
    structure = PDBParser(QUIET=True).get_structure("st", pdb_path)
    aa3to1 = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
        "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
        "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    }
    residues = [(res.id[1], aa3to1[res.get_resname()])
                for res in structure[0][chain] if res.id[0] == " " and res.get_resname() in aa3to1]
    residues.sort()
    return "".join(r[1] for r in residues), [r[0] for r in residues]


def map_contacts_to_anchor_frame(species, contacts, seq_anchor, seq_species, resnums_species):
    """Pairwise BLOSUM62 alignment of `species` against the anchor
    sequence, to translate its FMN-contact residue numbers into anchor
    sequence positions (0-based)."""
    blosum62 = substitution_matrices.load("BLOSUM62")
    alignments = pairwise2.align.localds(seq_anchor, seq_species, blosum62, -10, -0.5, one_alignment_only=True)
    if not alignments:
        return set()
    aln = alignments[0]

    resnum_to_local_idx = {rn: i for i, rn in enumerate(resnums_species)}
    contact_local_idx = {resnum_to_local_idx[rn] for rn in contacts if rn in resnum_to_local_idx}

    anchor_pos, species_pos, mapped = -1, -1, set()
    for a, b in zip(aln.seqA, aln.seqB):
        if a != "-":
            anchor_pos += 1
        if b != "-":
            species_pos += 1
            if a != "-" and species_pos in contact_local_idx:
                mapped.add(anchor_pos)
    return mapped


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Computing raw distance-based FMN contacts per species...")
    raw_contacts, sequences, resnums = {}, {}, {}
    for species, info in STRUCTURES.items():
        pdb_path = download_pdb(info["pdb_id"], OUT_DIR)
        if pdb_path is None:
            continue
        contacts = raw_contacts_to_fmn(pdb_path, info["chain"])
        raw_contacts[species] = contacts
        seq, resnum_list = get_sequence(pdb_path, info["chain"])
        sequences[species] = seq
        resnums[species] = resnum_list
        print(f"  {species} ({info['pdb_id']}): {len(contacts)} contact residues")

    if ANCHOR_SPECIES not in sequences:
        print(f"ERROR: could not build the anchor sequence for {ANCHOR_SPECIES}")
        return

    print(f"\nMapping contacts onto the {ANCHOR_SPECIES} reference frame (BLOSUM62)...")
    mapped_contacts = {ANCHOR_SPECIES: set(range(len(sequences[ANCHOR_SPECIES])))}  # placeholder, refined below
    anchor_positions_by_species = {}
    for species in STRUCTURES:
        if species == ANCHOR_SPECIES:
            # Anchor's own contacts, expressed directly as 0-based local positions
            resnum_to_idx = {rn: i for i, rn in enumerate(resnums[ANCHOR_SPECIES])}
            anchor_positions_by_species[species] = {
                resnum_to_idx[rn] for rn in raw_contacts[ANCHOR_SPECIES] if rn in resnum_to_idx
            }
            continue
        if species not in sequences:
            continue
        anchor_positions_by_species[species] = map_contacts_to_anchor_frame(
            species, raw_contacts[species], sequences[ANCHOR_SPECIES], sequences[species], resnums[species]
        )

    binding_site_union = set()
    for positions in anchor_positions_by_species.values():
        binding_site_union |= positions
    print(f"\nUnion of binding-site positions across all species "
          f"(anchor frame): {len(binding_site_union)}")

    print("\n=== Table 5-style overlap, relative to L. major ===")
    if "lmaj" in anchor_positions_by_species:
        lmaj_positions = anchor_positions_by_species["lmaj"]
        for other in ("hsap", "pfal"):
            if other not in anchor_positions_by_species:
                continue
            overlap = lmaj_positions & anchor_positions_by_species[other]
            pct = 100 * len(overlap) / len(lmaj_positions) if lmaj_positions else 0
            print(f"  lmaj contacts also present in {other}: {len(overlap)}/{len(lmaj_positions)} ({pct:.1f}%)")

    print("\nRunning PLIP on each structure (bond-type classification)...")
    for species, info in STRUCTURES.items():
        pdb_path = OUT_DIR / f"{info['pdb_id'].lower()}.pdb"
        if pdb_path.exists():
            run_plip(pdb_path)


if __name__ == "__main__":
    main()
