# 05 — InterPro Family vs. Domain comparison

Compares the InterPro `Family` and `Domain` entry categories across the
29 species of the network, to inform a future migration of the network's
functional annotation from PFAM to InterPro (Objective 4 of the thesis).

## Scripts

### `family_domain_coverage.py`

Main script. Parses the per-species InterProScan TSVs and computes:

- **Pillar 1 — Coverage**: of the proteins annotated with Family, how
  many also have a Domain.
- **Pillar 2a — Heterogeneity**: within a given Family, do its members
  share the same Domain, different Domains, or no Domain at all.
- **Pillar 2b — Positional overlap**: for proteins that have both
  annotations, do the Family and Domain windows share exact coordinates,
  is the Domain nested inside the Family, or do they not overlap.

It also writes `all_family_domain_windows.tsv` and `protein_lengths.tsv`,
which are the inputs for `biolip_enrichment.py` below.

### `domain_architecture_diagram.py`

Generates the schematic domain-architecture figures used to illustrate
individual Family/Domain cases (e.g. Figure 18 in the thesis), from the
Pillar 2b output.

### `biolip_enrichment.py`

Validates the Family/Domain comparison against ligand-binding data from
BioLiP: are binding residues over- or under-represented inside Domain vs.
Family windows, relative to what their size alone would predict.

This script excludes residues whose `uniprot_id` maps ambiguously to more
than one raw protein (the same ID-collision issue handled in the Gene
Mapper step) — an earlier version of this check, without that correction,
gave different enrichment values than the ones reported in the thesis.

## Results reported in the thesis

- 283,070 proteins with at least one Pfam-derived InterPro annotation;
  82,395 with a Family entry, 205,167 with a Domain entry, 4,492 with
  both (5.5% of those with Family).
- Of the Families surveyed: 89.7% with no associated Domain, 9.3%
  heterogeneous, 1.0% homogeneous.
- Of 6,416 Family–Domain pairs (proteins with both annotations): 0% exact
  coordinate match, 5.3% nested, 94.6% no overlap.
- General BioLiP enrichment: Domain 1.17, Family 1.07, unannotated 0.63
  (n = 970,826 binding residues).
- Restricted to proteins with both annotations, clean uniprot_id only:
  Domain 1.09, Family 1.22 (n = 38,362 binding residues).
