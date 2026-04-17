# Module 03 — AlphaFold Confidence Scores × PFAM Domains

Integrates per-residue AlphaFold structural confidence scores with PFAM domain annotations from InterProScan. For each annotated domain in each protein, the mean confidence score over the domain's residue range is computed and stored. Results are aggregated across all species to enable PFAM-level quality assessment.

AlphaFold proteomes are downloaded from Google Cloud:
```bash
gsutil -m cp gs://public-datasets-deepmind-alphafold-v4/proteomes/proteome-tax_id-[TAX_ID]-*_v4.tar .
tar xvf proteome-tax_id-[TAX_ID]-0_v4.tar
```
Tax IDs are obtained from the [NCBI Taxonomy Browser](https://www.ncbi.nlm.nih.gov/taxonomy).

## Pipeline overview

```
InterProScan output (.txt)
        │
        ▼
02_remap_ids_to_uniprot.py       ← species-specific ID mapping (see below)
        │
        ▼
01_match_interpro_to_alphafold.py
        │
        ▼
03_calculate_confidence_scores.py
        │
        ▼
04_merge_all_species.py
        │
        ├─▶ 05_pfam_aggregate_stats.py
        ├─▶ 06_plot_confidence_distributions.py
        └─▶ 07_filter_by_pfam_list.py
```

## Scripts

| # | Script | Description |
|---|--------|-------------|
| 1 | `01_match_interpro_to_alphafold.py` | Join InterProScan annotations with AlphaFold JSON paths |
| 2 | `02_remap_ids_to_uniprot.py` | Replace species IDs with UniProt IDs using mapping file |
| 3 | `03_calculate_confidence_scores.py` | Compute mean confidence score per PFAM domain |
| 4 | `04_merge_all_species.py` | Concatenate per-species tables into one file |
| 5 | `05_pfam_aggregate_stats.py` | Compute per-PFAM mean score, count, and species coverage |
| 6 | `06_plot_confidence_distributions.py` | Generate all distribution figures |
| 7 | `07_filter_by_pfam_list.py` | Filter results by a list of PFAMs of interest |

## Species-specific ID mapping

Some species required non-standard mapping strategies before running the main pipeline:

| Species | Script | Reason |
|---------|--------|--------|
| *Onchocerca volvulus* (ovo) | `species_specific/idmapping_ovo.py` | Not in UniProt mapper; requires GFF3 → WBGene → UniProt two-step mapping |
| *Loa loa* (loa) | `species_specific/idmapping_loa.py` | Uses Ensembl Metazoa BioMart to obtain Gene → Transcript → UniProt mapping |
| *T. cruzi* dm28c (tcru) | `species_specific/idmapping_tcru.py` | Multiple UniProt IDs per gene; applies curated > longest JSON > first ID priority rule |

## Usage

```bash
# Per species (after ID remapping if needed):
python scripts/01_match_interpro_to_alphafold.py \
  --interpro_dir  data/atha/interpro/ \
  --alphafold_dir data/atha/alphafold/ \
  --output        data/atha/atha_table_joined.txt \
  --species       atha

python scripts/03_calculate_confidence_scores.py \
  --table         data/atha/atha_table_joined.txt \
  --alphafold_dir data/atha/alphafold/ \
  --output        data/atha/atha_confidenceScore.txt

# Cross-species:
python scripts/04_merge_all_species.py \
  --pattern "data/*/*_confidenceScore.txt" \
  --output  data/all_species_confidence.tsv

python scripts/05_pfam_aggregate_stats.py \
  --input  data/all_species_confidence.tsv \
  --output data/PFAM_confidence_avg.tsv

python scripts/06_plot_confidence_distributions.py \
  --all_species data/all_species_confidence.tsv \
  --pfam_stats  data/PFAM_confidence_avg.tsv \
  --outdir      figures/
```
