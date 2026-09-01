#%%
"""
Mapeo de identificadores OrthoMCL -> UniProt para una especie.

Usa GeneMapper (herramienta desarrollada en el labo, citada en la tesis)
para matchear por id exacto, luego por hash MD5 de secuencia, y por último
por BLAST (umbral de cobertura 0.99). Se corre una vez por especie:
alcanza con cambiar SPECIES_CODE y UNIPROT_TAXON_ID abajo.
"""

from modules.gene_mapper import GeneMapper
import seqpandas as spd
import pandas as pd
import numpy as np
from pathlib import Path

#%%
# --- Config: lo único que cambia entre especies ---
SPECIES_CODE = "bmaa"        # código de especie (ver tabla en el README del repo)
UNIPROT_TAXON_ID = "6279"    # id taxonómico de UniProt para esta especie

BASE_DIR = Path("/big/lab/ssneider/ssneider-env/TDR_Targets7.1")
GENE_MAPPER_DIR = BASE_DIR / "gene_mapper" / "Ortho_vs_Uniprot" / SPECIES_CODE

fasta_new = BASE_DIR / "actualizados_short" / "sequences" / f"{SPECIES_CODE}_aa_seqs_OrthoMCL-7.fasta"
fasta_old = GENE_MAPPER_DIR / f"uniprotkb_organism_id_{UNIPROT_TAXON_ID}_{SPECIES_CODE}.fasta"
tsv_uniprot = GENE_MAPPER_DIR / f"uniprotkb_organism_id_{UNIPROT_TAXON_ID}_{SPECIES_CODE}.tsv"

#%%
mapper = GeneMapper()

#%%
# OrthoMCL es el query (más corto, "new") y UniProt es la referencia (grande, "old").
# Por cada gen de OrthoMCL buscamos su match en UniProt.

# seqpandas parsea el fasta a un DataFrame usando SeqRecord de Biopython como base.
df_old_crudo = spd.read_seq(str(fasta_old), format='fasta')
df_new_crudo = spd.read_seq(str(fasta_new), format='fasta')

#%%
# Nos quedamos solo con secuencia, id y descripción
columnas = ['_seq', 'id', 'description']

df_old = df_old_crudo[columnas].copy()
df_new = df_new_crudo[columnas].copy()

#%%
# Limpieza de IDs: en ambos archivos el accession real está en split("|")[1]
#   UniProt:  "tr|A0A024QYR6|A0A024QYR6_HUMAN" -> "A0A024QYR6"
#   OrthoMCL: "hsap|ENSG00000136830"           -> "ENSG00000136830"
df_old["id"] = df_old["id"].apply(lambda x: x.split("|")[1])
df_new["id"] = df_new["id"].apply(lambda x: x.split("|")[1])

#%%
df_old.to_csv(f'uniprot_{SPECIES_CODE}_1.csv', index=False)
df_new.to_csv(f'orthomcl_{SPECIES_CODE}_1.csv', index=False)

#%%
# Mapeo: query = df_new (OrthoMCL), referencia = df_old (UniProt)
genemap = mapper.genes_mapper(df_new, df_old, blast_coverage_threshold=0.99)

#%%
mapped      = genemap['mapped']
unmapped    = genemap['unmapped']
problematic = genemap['problematic']

#%%
# Reseteamos índices para que coincidan con los idx que devuelve GeneMapper
df_new_idx = df_new.reset_index(drop=True)
df_old_idx = df_old.reset_index(drop=True)

#%%
def resolver_id(idx, df):
    """Devuelve el id de df en la posición idx, o NaN si idx es <NA>/NaN."""
    if pd.isna(idx):
        return np.nan
    return df.iloc[int(idx)]['id']

# ID del query (OrthoMCL), viene de idx_gtm
mapped['orthomcl_id'] = mapped['idx_gtm'].apply(lambda i: resolver_id(i, df_new_idx))

# IDs de UniProt según cada estrategia de match
mapped['uniprot_id_by_id']    = mapped['idx_db_id'].apply(lambda i: resolver_id(i, df_old_idx))
mapped['uniprot_id_by_md5']   = mapped['idx_db_hash_md5'].apply(lambda i: resolver_id(i, df_old_idx))
mapped['uniprot_id_by_blast'] = mapped['idx_db_blast'].apply(lambda i: resolver_id(i, df_old_idx))

# ID consolidado: prioridad id exacto > md5 > blast
mapped['uniprot_id'] = (
    mapped['uniprot_id_by_id']
      .fillna(mapped['uniprot_id_by_md5'])
      .fillna(mapped['uniprot_id_by_blast'])
)

# De qué método vino cada match (para control de calidad)
def metodo_match(row):
    if pd.notna(row['uniprot_id_by_id']):    return 'id'
    if pd.notna(row['uniprot_id_by_md5']):   return 'md5'
    if pd.notna(row['uniprot_id_by_blast']): return 'blast'
    return 'none'
mapped['match_method'] = mapped.apply(metodo_match, axis=1)

# Cargamos el TSV de UniProt para sumar nombres legibles al resultado
uniprot_meta = pd.read_csv(tsv_uniprot, sep='\t')
uniprot_meta = uniprot_meta.rename(columns={
    'Entry': 'uniprot_id',
    'Entry Name': 'uniprot_entry_name',
    'Protein names': 'uniprot_protein_names',
    'Reviewed': 'uniprot_reviewed',
})

#%%
mapped_enriched = mapped.merge(uniprot_meta, on='uniprot_id', how='left')

# Para unmapped/problematic alcanza con resolver el id de OrthoMCL,
# así los CSV quedan legibles.
unmapped['orthomcl_id']    = unmapped['idx_gtm'].apply(lambda i: resolver_id(i, df_new_idx))
problematic['orthomcl_id'] = problematic['idx_gtm'].apply(lambda i: resolver_id(i, df_new_idx))

#%%
mapped_enriched.to_csv(GENE_MAPPER_DIR / 'mapped_orthomcl_to_uniprot.csv', index=False)
unmapped.to_csv(GENE_MAPPER_DIR / 'unmapped_orthomcl_to_uniprot.csv', index=False)
problematic.to_csv(GENE_MAPPER_DIR / 'problematic_orthomcl_to_uniprot.csv', index=False)

# Chequeo rápido de resultados
print(mapped_enriched['match_method'].value_counts())
print(f"Total mapped:      {len(mapped_enriched)}")
print(f"Total unmapped:    {len(unmapped)}")
print(f"Total problematic: {len(problematic)}")
