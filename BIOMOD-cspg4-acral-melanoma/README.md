# BIOMOD 2026 «Caja de Pandora»: validación bioinformática de CSPG4 en melanoma acral

Evaluación de la expresión de ARNm de **CSPG4** en piel sana, melanoma cutáneo y melanoma acral de pacientes mexicanos, con datos públicos de RNA-seq (GTEx, TCGA-SKCM) y la cohorte publicada por Basurto-Lozada et al. (*Nature*, 2026). Responsable: Rubí.

> **Hipótesis (no conclusión):** la expresión de CSPG4 en melanoma acral es mayor que en piel sana y comparable a la de melanoma cutáneo. Se reporta el resultado sea cual sea, con tamaños de efecto.

## Cómo está organizado el análisis

| Nivel | Qué | Script | Validez |
|---|---|---|---|
| A | CSPG4 en 44 melanomas acrales mexicanos (ulceración, estadio, cluster) | `01_nivelA_acral_mexicano.R` | Análisis interno, sin batch entre fuentes |
| B | CSPG4 en TCGA-SKCM, primario vs metástasis (DESeq2) | `02_nivelB_tcga_skcm.R` | Formal dentro de TCGA |
| C | CSPG4 en piel sana (GTEx) y otros tejidos | `03_nivelC_gtex.R` | Referencia descriptiva |
| — | Integración de A + B + C y figura principal | `04_integracion_figuras.R` | **Exploratoria** (inter-estudio, sin corrección de batch) |

**Por qué no se corrige batch entre fuentes:** cada grupo proviene de una sola fuente, así que fuente y condición son la misma variable y la corrección borraría la diferencia biológica. Wang et al. (*Sci Data* 2018) pudieron corregir entre GTEx y TCGA porque comparten tejidos normales; melanoma quedó fuera de su conjunto.

## Cómo correrlo
Desde la raíz del repositorio:
```
Rscript scripts/01_nivelA_acral_mexicano.R
Rscript scripts/02_nivelB_tcga_skcm.R
Rscript scripts/03_nivelC_gtex.R
Rscript scripts/04_integracion_figuras.R
```
Antes: descargar los datos (ver `data/README.md`) y revisar `scripts/config.R` (nombres de archivos y columnas).

## Estado de las pruebas (leer antes de usar)
| Script | Estado |
|---|---|
| 01, 03, 04, `helpers.R` | Probados con **datos simulados** (R 4.3.3, ggplot2). Funcionan, pero **no se han corrido con los datos reales**. |
| 02 | **Solo revisada la sintaxis.** Necesita descargar TCGA-SKCM y paquetes de Bioconductor. Revisar cada mensaje al correrlo. |
| Prueba de Dunn | Implementada a mano; con 2 grupos reproduce exactamente el p de `wilcox.test` (verificado). |

Paquetes: R base, `ggplot2`; para el script 02: `TCGAbiolinks`, `SummarizedExperiment`, `DESeq2`, `apeglm` (Bioconductor).

## Documentación
- `docs/VERIFICACION.md`: lista de cosas por confirmar antes de citar o interpretar.
- `docs/PLAN_B_DATOS.md`: fuentes alternativas si no llega el acceso a los datos.
- `docs/BITACORA.md`: registro de cada sesión de análisis.
- Protocolo completo: «Protocolo 1: Validación bioinformática de CSPG4 en melanoma acral (BIOMOD 2026)» (exportar a Markdown y guardar en `docs/`).

## Fuentes principales
- Basurto-Lozada P, et al. *Nature* 2026;651:221-230. Código de los autores: https://github.com/CGBio-Lab/Mex-acral-exomes-transcriptomes
- Zubair H. protocols.io 2026, doi:10.17504/protocols.io.kqdg3mn6el25/v1 (flujo DESeq2 con TCGA; validado por su autor en TCGA-LUAD).
- Love MI, Huber W, Anders S. *Genome Biol* 2014;15:550. Colaprico A, et al. *Nucleic Acids Res* 2016;44:e71.
- NCI GDC, mRNA Analysis Pipeline: https://docs.gdc.cancer.gov
- Wang Q, et al. *Sci Data* 2018, doi:10.1038/sdata.2018.61.

## Aviso sobre datos
Este repositorio es público. **No se suben datos crudos ni procesados** (ver `.gitignore`). La tabla por muestra de los pacientes mexicanos tampoco se publica hasta contar con el visto bueno de las autoras y autores del estudio original.
