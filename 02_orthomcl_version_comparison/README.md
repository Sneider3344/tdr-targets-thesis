# Module 02 — OrthoMCL v6 vs v7: Model Genome Comparison

Determines whether OrthoMCL v6 or v7 produces more comprehensive ortholog groupings, using *Saccharomyces cerevisiae* (scer) and *Trypanosoma brucei* (tbrt) as reference genomes. Both have known gene counts (9,660 for T. brucei; 5,907 for S. cerevisiae), making them reliable benchmarks.

The expected result is that v7 produces fewer, larger orthogroups — meaning it identifies more associations between proteins than v6.

## Scripts

| # | Script | Description |
|---|--------|-------------|
| 1 | `01_filter_fasta_by_species.py` | Extract scer/tbrt sequences from full OrthoMCL FASTA |
| 2 | `02_split_fasta_by_species.py` | Split filtered FASTA into one file per species |
| 3 | `03_fasta_headers_to_csv.py` | Parse headers into structured CSV (gene_id, protein, orthology_group) |
| 4 | `04_plot_genes_per_orthogroup.py` | Bar chart of orthogroup size distribution |

## Usage

```bash
# Step 1 — filter from compressed database
python scripts/01_filter_fasta_by_species.py \
  --input aa_seqs_OrthoMCL-6.21.fasta.gz \
  --output filtered_scer_tbrt.fasta \
  --species scer tbrt

# Step 2 — split by species
python scripts/02_split_fasta_by_species.py \
  --input filtered_scer_tbrt.fasta \
  --species scer tbrt \
  --outdir data/split/

# Step 3 — parse headers to CSV (run for each species/version combination)
python scripts/03_fasta_headers_to_csv.py \
  --input data/split/scer_orthogroups.fasta \
  --version V7 \
  --output data/scer_OrtgroupsV7.csv

# Step 4 — plot
python scripts/04_plot_genes_per_orthogroup.py \
  --input data/scer_OrtgroupsV7.csv \
  --title "S. cerevisiae — OrthoMCL v7" \
  --output figures/scer_v7_distribution.png
```

## Notes

- This analysis informed the decision to adopt OrthoMCL v7 as the reference database for all downstream analyses in this project.
- The comparison with OrthoFinder is documented separately.
