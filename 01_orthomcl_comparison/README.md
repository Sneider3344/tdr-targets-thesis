# Module 01 — OrthoMCL v6/v7 Group Stability

Evaluates how stable ortholog group assignments are across OrthoMCL database versions for the four non-core genomes (kpm, loa, ovo, sao). For each species, sequences are merged across versions, collapsed into per-orthogroup rows, and classified as conserved or reassigned based on whether the numeric group identifier is consistent between v6 and v7.

## Scripts

| # | Script | Input | Output |
|---|--------|-------|--------|
| 1 | `01_merge_orthomcl_versions.py` | v6 TSV, v7 TSV | Combined TSV |
| 2 | `02_group_by_orthogroup.py` | Raw OrthoMCL TSV | Grouped TSV + ungrouped TSV |
| 3 | `03_compare_orthogroup_versions.py` | v6 grouped, v7 grouped | Conserved TSV + reassigned TSV |

## Usage

```bash
# Step 1 — merge
python scripts/01_merge_orthomcl_versions.py \
  --old data/kpm_orthomcl6.tsv --new data/kpm_orthomcl7.tsv \
  --out data/kpm_combined.tsv

# Step 2 — group (repeat for v6 and v7)
python scripts/02_group_by_orthogroup.py \
  --input data/kpm_orthomcl7.tsv --out data/kpm_v7_grouped.tsv \
  --ungrouped data/kpm_v7_nongrouped.tsv --prefix kpm

# Step 3 — compare
python scripts/03_compare_orthogroup_versions.py \
  --old data/kpm_v6_grouped.tsv --new data/kpm_v7_grouped.tsv \
  --same data/kpm_conserved.tsv --diff data/kpm_reassigned.tsv
```
