# TDR Targets Therapeutic Target Identification — Thesis Pipeline

Bioinformatics pipeline developed for the identification and curation of therapeutic targets in neglected tropical disease (NTD) pathogens, as part of an undergraduate Biotechnology thesis. The project updates and extends the TDR Targets chemogenomic platform (v7) across five specific objectives: biological/chemical layer curation, identifier interoperability, AlphaFold structural integration, PFAM-to-InterPro annotation comparison, and PocketVec/US-align structural validation of binding pockets.

---

## Biological context

Neglected tropical diseases affect over one billion people worldwide, yet drug development remains underfunded. Computational target identification — combining orthology data, structural predictions, and multi-layer biological networks — allows prioritising proteins that are essential to the pathogen, absent or dissimilar in the human host, and structurally well-characterised enough to be druggable.

This pipeline processes data from 29 parasite and model-organism genomes in the network, including *Trypanosoma cruzi* (CL Brener and Dm28c), *Brugia malayi*, *Onchocerca volvulus*, *Leishmania* spp., *Plasmodium falciparum*, *Homo sapiens*, and *Mus musculus*, among others. Five species (*Loa loa*, *Onchocerca volvulus*, *Echinococcus granulosus*, *Klebsiella pneumoniae*, *Staphylococcus aureus*) are not part of OrthoMCL's native catalog; their orthology groups are assigned via Diamond BLASTp against the existing network.

---

## Repository structure

```
tdr-targets-thesis/
├── 02_orthomcl_version_comparison/  # Model genome benchmarking (scer, tbrt): OrthoMCL v6 vs v7
├── 03_alphafold_pfam/               # AlphaFold confidence × PFAM domain integration
│   └── scripts/species_specific/    # Non-standard ID mapping pipelines
├── 04_gene_mapper/                  # OrthoMCL <-> UniProt identifier resolution
├── 05_interpro_family_domain/       # InterPro Family vs. Domain comparison, BioLiP validation
├── 06_pocketvec/                    # PocketVec / US-align structural pocket validation
│   ├── primera_seleccion/           # Stage 1: 2 pilot OG7 groups
│   └── segunda_seleccion/           # Stage 2: 3 finalist OG7 groups (sequence-divergence selection)
├── requirements.txt
├── .gitignore
└── README.md
```

Numbering starts at `02`: an earlier module (`01_orthomcl_comparison`) compared native OrthoMCL v6/v7 group files for the five non-core species, an approach later superseded by the Diamond BLASTp assignment described above, and was removed.

Each module has its own README with a full usage guide. Modules are numbered in the order they were developed; 04-06 depend on outputs from Gene Mapper (04) and, for 06, on the candidate lists produced in 05.

---

## Requirements

- Python ≥ 3.10
- See `requirements.txt` for Python dependencies
- External tools: `gsutil` (AlphaFold download), `InterProScan` (PFAM annotation), `MAFFT` (sequence alignment), `US-align` (structural alignment), `OpenBabel`, `PLIP`, PocketVec + rDock (docking-based pocket descriptors)

```bash
pip install -r requirements.txt
```

---

## Quick start

```bash
git clone https://github.com/Sneider3344/tdr-targets-thesis.git
cd tdr-targets-thesis
pip install -r requirements.txt
```

Then follow the README in each module directory, in numerical order.

---

## Development notes

- Raw data files (BioLiP, ChEMBL, AlphaFold JSON/PDB, InterProScan outputs, OrthoMCL TSVs) are excluded from version control (see `.gitignore`).
- Some species required non-standard ID mapping before AlphaFold integration (see `03_alphafold_pfam/scripts/species_specific/`).
- Modules 04-06 depend on `GeneMapper` and `seqpandas`, in-house/external tools cited in the thesis and not included in this repository.
- This is an active thesis project — additional modules and refinements will be added as the analysis progresses.
