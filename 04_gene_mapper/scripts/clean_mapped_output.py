"""
Limpieza del resultado de map_orthomcl_to_uniprot.py.

Nos quedamos solo con las columnas que se usan más adelante en el análisis:
ID de OrthoMCL, ID de UniProt y longitud de la secuencia. Se corre después
del script de mapeo, para la misma especie.
"""
import pandas as pd
from pathlib import Path

SPECIES_CODE = "bmaa"  # mismo código usado en map_orthomcl_to_uniprot.py

BASE_DIR = Path("/big/lab/ssneider/ssneider-env/TDR_Targets7.1/gene_mapper/Ortho_vs_Uniprot") / SPECIES_CODE
in_path  = BASE_DIR / "mapped_orthomcl_to_uniprot.csv"
out_path = BASE_DIR / "mapped_clean.csv"

df = pd.read_csv(in_path)

clean = df[['orthomcl_id', 'uniprot_id', 'Length']].rename(columns={
    'orthomcl_id': 'OrthoMCL_ID',
    'uniprot_id':  'Uniprot_ID',
    'Length':      'Length',
})

clean.to_csv(out_path, index=False)
print(clean.head())
print(f"\nGuardado en: {out_path}")
