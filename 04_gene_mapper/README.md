# 04 — Gene Mapper: identifier interoperability

Resolves the correspondence between OrthoMCL and UniProt identifiers for
each species in the project, using GeneMapper (an in-house tool developed
in the lab, not included in this repository; see citation in the thesis).

## Mapping strategy

For each OrthoMCL gene, we look for its match in UniProt, in this order
of priority:

1. Exact ID match
2. Sequence hash (MD5) match
3. BLAST match (minimum coverage 0.99)

Results are split into three files: `mapped`, `unmapped`, and
`problematic`.

## Scripts

- `map_orthomcl_to_uniprot.py`: runs the full mapping for one species.
- `clean_mapped_output.py`: keeps only the final columns
  (`OrthoMCL_ID`, `Uniprot_ID`, `Length`) for downstream analysis.

## Usage

This was run individually for each species in the project by changing
`SPECIES_CODE` (and `UNIPROT_TAXON_ID` in the mapping script) at the top
of each script.

```
SPECIES_CODE = "bmaa"
```

## Dependencies not included

`GeneMapper` and `seqpandas` are lab tools already cited in the thesis;
they are not included in this repository.
