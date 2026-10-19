# Datos (no se suben a GitHub)

Esta carpeta se llena localmente. El `.gitignore` evita que los datos se publiquen.

## Nivel A: melanoma acral mexicano
- Fuente: Basurto-Lozada et al., *Nature* 2026;651:221-230.
- Los autores indican, en su repositorio (https://github.com/CGBio-Lab/Mex-acral-exomes-transcriptomes, archivo «Input files download»), una carpeta pública de Dropbox con los insumos de sus análisis:
  https://www.dropbox.com/scl/fo/okg557dsdosxw11hu7f25/AFrwyU62kR5k_kStb8iMXUU?rlkey=u3f58zc7hf51ca8gxto40hsro&st=w4m160bb&dl=0
- Archivos que el script de clustering de los autores carga (confirmar que están en la carpeta):
  - `acralMX_expression_data_update44s.rds`
  - `acralMX_metadata_cc_update.csv`
- Guardarlos en `data/raw/`.
- Counts crudos: EGA, dataset EGAD00001015756 (acceso controlado, requiere solicitud al comité de acceso).

## Nivel B: TCGA-SKCM
- Se descarga automáticamente con el script `02_nivelB_tcga_skcm.R` (TCGAbiolinks, GDC, flujo STAR - Counts). Requiere ~20 GB libres.

## Nivel C: GTEx
- Descargar de https://gtexportal.org (sección Downloads) los archivos de la versión 8 y guardarlos en `data/raw/`:
  - TPM de genes (`..._gene_tpm.gct.gz`)
  - Anotación de muestras (`..._SampleAttributesDS.txt`)
- Verificar los nombres exactos en el portal y, si difieren, ajustarlos en `scripts/config.R`.
