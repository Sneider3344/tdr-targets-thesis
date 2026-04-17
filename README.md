# TDR Targets Therapeutic Target Identification — Thesis Pipeline

Bioinformatics pipeline developed for the identification and curation of therapeutic targets in neglected tropical disease (NTD) pathogens, as part of an undergraduate Biotechnology thesis. The project updates and extends the TDR Targets platform (v7.1) by incorporating a new ortholog database version, structural confidence data from AlphaFold, and multi-layer network analysis across five or more parasite genomes.

---

## Biological context

Neglected tropical diseases affect over one billion people worldwide, yet drug development remains underfunded. Computational target identification — combining orthology data, structural predictions, and multi-layer biological networks — allows prioritising proteins that are essential to the pathogen, absent or dissimilar in the human host, and structurally well-characterised enough to be druggable.

This pipeline processes data from the following genomes (core OrthoMCL + non-core species assigned via DiamondBLAST):

| Code | Species |
|------|---------|
| kpm  | *Leishmania panamensis* |
| loa  | *Loa loa* |
| ovo  | *Onchocerca volvulus* |
| sao  | *Strongyloides* sp. |
| tcru | *Trypanosoma cruzi* dm28c |
| scer | *Saccharomyces cerevisiae* (reference) |
| tbrt | *Trypanosoma brucei* (reference) |

---

## Repository structure

```
tdr-targets-thesis/
├── 01_orthomcl_comparison/         # OrthoMCL v6 vs v7 group stability analysis
│   ├── scripts/
│   └── README.md
├── 02_orthomcl_version_comparison/ # Model genome benchmarking (scer, tbrt)
│   ├── scripts/
│   └── README.md
├── 03_alphafold_pfam/              # AlphaFold confidence × PFAM domain integration
│   ├── scripts/
│   │   └── species_specific/       # Non-standard ID mapping pipelines
│   └── README.md
├── requirements.txt
├── .gitignore
└── README.md
```

Each module has its own README with a full usage guide. Run modules in the order listed above, as later modules depend on outputs from earlier ones.

---

## Requirements

- Python ≥ 3.10
- See `requirements.txt` for Python dependencies
- External tools: `gsutil` (AlphaFold download), `InterProScan` (PFAM annotation)

```bash
pip install -r requirements.txt
```

---

## Quick start

```bash
git clone https://github.com/<your-username>/tdr-targets-thesis.git
cd tdr-targets-thesis
pip install -r requirements.txt
```

Then follow the README in each module directory in order.

---

## Development notes

- Raw data files (`data/raw/`) are excluded from version control (see `.gitignore`). These include OrthoMCL TSVs, AlphaFold JSON files, and InterProScan outputs.
- Some species required non-standard ID mapping before AlphaFold integration (see `03_alphafold_pfam/scripts/species_specific/`).
- This is an active thesis project — additional modules will be added as the analysis progresses.
