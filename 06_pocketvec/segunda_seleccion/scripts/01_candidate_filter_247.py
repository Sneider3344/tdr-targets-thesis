#!/usr/bin/env python3
"""
Second-stage candidate filter: 992 -> 247 OG7 groups.

Two parts:
1. Ligand categorization: classifies every BioLiP ligand as a universal
   cofactor, a confirmed KEGG metabolite, a confirmed DrugCentral drug, or
   unconfirmed -- using a small hardcoded cofactor list, KEGG's compound
   database (InChIKey matching), and DrugCentral's structure table.
2. Candidate shortlist: from the 992 corrected candidates, keeps groups
   with at least 2 crystallized species from at least 2 distinct taxonomic
   lineages, and at least 8 species present in the network overall.
"""
from pathlib import Path

import pandas as pd
from rdkit.Chem.inchi import InchiToInchiKey

IN_DIR  = "biolip_exploracion2"
OUT_DIR = "Nueva_estrategia_pocketvec"
KEGG_INCHI_FILE    = "/big/lab/ssneider/ssneider-env/kegg_compound/compound/compound.inchi"
KEGG_COMPOUND_FILE = "/big/lab/ssneider/ssneider-env/kegg_compound/compound/compound"
DRUGCENTRAL_FILE   = "/big/lab/ssneider/ssneider-env/drugcentral_structures.tsv"
LIGAND_FILE        = "/home/ssneider/disco-big/ssneider-env/TDR_Targets7.1/PocketVec/analisis_biolip/ligand.tsv.gz"

COFACTORS_UNIVERSAL = {
    "ATP", "ADP", "AMP", "GTP", "GDP", "GMP", "CTP", "CDP", "CMP", "UTP", "UDP", "UMP",
    "NAD", "NAP", "NDP", "FAD", "FMN", "SAM", "SAH", "COA", "ACO", "HEM", "HEC",
    "PLP", "PMP", "TPP", "THF", "FE2", "MG", "ZN", "CA", "MN", "CU", "NI",
    "PYR", "OXA", "CIT", "MAL", "FUM", "SUC", "AKG", "OXL", "PGA", "PEP", "FBP",
    "G3P", "2PG", "IMP", "GSH", "DUT", "DUP", "RBF", "FES", "SF4", "F3S",
}

# Lineage assignment per species code, used to require BioLiP coverage
# across at least 2 distinct lineages (not just 2 species of the same kind)
LINEAGE = {
    "hsap": "vertebrate", "mmus": "vertebrate",
    "dmel": "model_invertebrate", "cele": "model_invertebrate",
    "scer": "fungus", "calb": "fungus",
    "atha": "plant", "osat": "plant",
    "ddis": "model_amoeba",
    "ecol": "bacteria", "mtub": "bacteria", "kpm": "bacteria", "sao": "bacteria",
    "ldon": "kinetoplastid", "lmaj": "kinetoplastid",
    "tbrt": "kinetoplastid", "tcru": "kinetoplastid",
    "pfal": "apicomplexan", "pvip": "apicomplexan", "tgon": "apicomplexan",
    "bmaa": "nematode", "loa": "nematode", "ovo": "nematode",
    "egr": "cestode",
    "gass": "other_parasite", "ehia": "other_parasite", "tvag": "other_parasite",
    "gmur": "other",
}
IS_PARASITE = {"ldon", "lmaj", "tbrt", "tcru", "pfal", "pvip", "tgon",
                "bmaa", "loa", "ovo", "egr", "gass", "ehia", "tvag"}
IS_HOST = {"hsap", "mmus"}


def classify_ligands():
    """Parses KEGG compound + DrugCentral to categorize every BioLiP ligand."""
    print("Parsing KEGG compound...")
    kegg_is_metabolite = {}
    entry, has_pathway, has_enzyme, has_br08001, section = None, False, False, False, None
    with open(KEGG_COMPOUND_FILE, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("///"):
                if entry:
                    kegg_is_metabolite[entry] = has_pathway or has_enzyme or has_br08001
                entry, has_pathway, has_enzyme, has_br08001, section = None, False, False, False, None
                continue
            if not line.startswith(" ") and line.strip():
                parts = line.split(None, 1)
                section = parts[0]
                rest = parts[1].strip() if len(parts) > 1 else ""
            else:
                rest = line.strip()
            if section == "ENTRY":
                entry = rest.split()[0]
            elif section == "PATHWAY" and rest:
                has_pathway = True
            elif section == "ENZYME" and rest:
                has_enzyme = True
            elif section == "BRITE" and "br08001" in rest:
                has_br08001 = True

    print("Converting KEGG InChI to InChIKey...")
    inchikey_is_metabolite, prefix_is_metabolite = {}, {}
    with open(KEGG_INCHI_FILE, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or "\t" not in line:
                continue
            kid, inchi = line.split("\t", 1)
            kid = kid.strip()
            if kid not in kegg_is_metabolite:
                continue
            try:
                ikey = InchiToInchiKey(inchi.strip())
            except Exception:
                continue
            if not ikey:
                continue
            is_met = kegg_is_metabolite[kid]
            inchikey_is_metabolite[ikey] = is_met
            prefix = ikey[:14]
            prefix_is_metabolite[prefix] = prefix_is_metabolite.get(prefix, False) or is_met

    print("Loading DrugCentral...")
    df_dc = pd.read_csv(DRUGCENTRAL_FILE, sep="\t")
    df_dc["InChIKey"] = df_dc["InChIKey"].astype(str).str.strip()
    drugcentral_inchikeys = set(df_dc["InChIKey"].dropna())
    drugcentral_prefixes = {k[:14] for k in drugcentral_inchikeys if len(k) >= 14}

    print("Loading ligand.tsv.gz...")
    lig_cols = ["lig_id", "formula", "inchi", "inchikey", "smiles", "name", "chebi", "drugbank", "zinc"]
    import gzip
    with gzip.open(LIGAND_FILE, "rt", encoding="utf-8", errors="replace") as f:
        df_lig = pd.read_csv(f, sep="\t", header=None, names=lig_cols, on_bad_lines="warn", engine="python")
    df_lig["lig_id"] = df_lig["lig_id"].astype(str).str.strip().str.upper()
    df_lig["inchikey"] = df_lig["inchikey"].astype(str).str.strip()
    df_lig = df_lig[df_lig["lig_id"] != "#CCD"].copy()

    def classify(row):
        lig, ikey = row["lig_id"], row["inchikey"]
        prefix = ikey[:14] if len(ikey) >= 14 else ""
        if lig in COFACTORS_UNIVERSAL:
            return "cofactor_universal"
        if (ikey in inchikey_is_metabolite and inchikey_is_metabolite[ikey]) or \
           (prefix in prefix_is_metabolite and prefix_is_metabolite[prefix]):
            return "metabolito_confirmado_kegg"
        if (ikey in drugcentral_inchikeys) or (prefix in drugcentral_prefixes):
            return "farmaco_confirmado_drugcentral"
        return "no_confirmado"

    df_lig["categoria"] = df_lig.apply(classify, axis=1)
    print(df_lig["categoria"].value_counts())
    return dict(zip(df_lig["lig_id"], df_lig["categoria"]))


def parse_species(s):
    if pd.isna(s):
        return []
    return [x.strip() for x in str(s).split(",") if x.strip()]


def filter_247_candidates(lig_categoria):
    df = pd.read_csv(f"{IN_DIR}/candidatos_final.tsv", sep="\t")

    # Each candidate's dominant ligand (from stage 2's enrichment) tagged
    # with its confirmed category
    df["top_lig"] = df["top_ligandos"].astype(str).str.split(",").str[0].str.strip()
    df["top_lig_categoria"] = df["top_lig"].map(lig_categoria).fillna("no_confirmado")

    df["lista_biolip"] = df["especies_biolip"].apply(parse_species)
    df["lista_red"] = df["especies_red"].apply(parse_species)
    df["n_esp_biolip"] = df["lista_biolip"].apply(len)
    df["n_esp_red"] = df["lista_red"].apply(len)
    df["linajes_biolip"] = df["lista_biolip"].apply(lambda L: sorted(set(LINEAGE.get(e, "other") for e in L)))
    df["n_linajes_biolip"] = df["linajes_biolip"].apply(len)
    df["tiene_hosp_biolip"] = df["lista_biolip"].apply(lambda L: any(e in IS_HOST for e in L))
    df["tiene_para_biolip"] = df["lista_biolip"].apply(lambda L: any(e in IS_PARASITE for e in L))

    print("Total OG7:", len(df))
    print("  >=1 species in BioLiP:", (df["n_esp_biolip"] >= 1).sum())
    print("  >=2 species in BioLiP:", (df["n_esp_biolip"] >= 2).sum())
    print("  >=2 species AND >=2 lineages:", ((df["n_esp_biolip"] >= 2) & (df["n_linajes_biolip"] >= 2)).sum())

    cand = df[
        (df["n_esp_biolip"] >= 2) &
        (df["n_linajes_biolip"] >= 2) &
        (df["n_esp_red"] >= 8)
    ].copy()

    cand["ligando_defendible"] = cand["top_lig_categoria"].isin(
        ["cofactor_universal", "metabolito_confirmado_kegg", "farmaco_confirmado_drugcentral"]
    )

    # Heuristic score used only to order the shortlist for manual review,
    # not a final selection metric
    cand["score_candidato"] = (
        cand["n_esp_biolip"] * 2 +
        cand["n_linajes_biolip"] * 3 +
        cand["tiene_hosp_biolip"].astype(int) * 4 +
        cand["tiene_para_biolip"].astype(int) * 4 +
        cand["ligando_defendible"].astype(int) * 3
    )
    cand = cand.sort_values(
        ["tiene_hosp_biolip", "tiene_para_biolip", "score_candidato"], ascending=False
    ).reset_index(drop=True)

    cols = ["og7", "n_esp_biolip", "n_linajes_biolip", "tiene_hosp_biolip", "tiene_para_biolip",
            "especies_biolip", "n_esp_red", "top_lig", "top_lig_categoria",
            "ligando_defendible", "score_candidato"]

    Path(OUT_DIR).mkdir(exist_ok=True)
    cand[cols].to_csv(f"{OUT_DIR}/candidatos_divergencia_stage1.tsv", sep="\t", index=False)
    print(f"\nCandidates passing the 247 filter: {len(cand)}")
    print(cand[cols].head(30).to_string(index=False))


def main():
    lig_categoria = classify_ligands()
    filter_247_candidates(lig_categoria)


if __name__ == "__main__":
    main()
