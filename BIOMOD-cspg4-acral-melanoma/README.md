# BIOMOD 2026 «Caja de Pandora» · Protocolo 1: expresión de CSPG4 en melanoma acral mexicano

Responsable: Rubí

> **Estado: borrador de trabajo.** Los resultados son exploratorios y hay pendientes de terceros (ver «Pendientes»). No son una conclusión definitiva.

## Qué se evalúa

Si el patrón de expresión del ARNm de **CSPG4** (ENSG00000173546) respalda o limita su uso como blanco del direccionamiento activo del nanovehículo del equipo. Se compara la expresión, en escala log2(TPM + 1), entre tres condiciones y se contextualiza con otros tejidos normales:

| Grupo | Condición | Fuente |
|---|---|---|
| G1a / G1b | Piel sana no acral (no expuesta / expuesta al sol) | GTEx v8 |
| G2 | Melanoma cutáneo primario | TCGA-SKCM (GDC) |
| G3 | Melanoma acral primario (44 muestras; 65 en sensibilidad) | Cohorte mexicana, Basurto-Lozada et al. (referencia completa por verificar; el repositorio de los autores indica 2024) |
| Contexto | 54 tejidos normales (medianas) | GTEx v8 |

Código de los autores de la cohorte: https://github.com/CGBio-Lab/Mex-acral-exomes-transcriptomes

## Resultados principales (exploratorios)

Cada condición viene de un solo estudio, así que **fuente y condición están confundidas**: no se aplicó corrección de lote y las diferencias entre fuentes no se pueden atribuir solo a la biología.

- Medianas de log2(TPM + 1): piel no expuesta 3.576, piel expuesta 3.830, melanoma cutáneo primario 7.673, melanoma acral primario 5.511.
- Acral frente a piel sana: mayor, con efecto mediano a grande (delta de Cliff −0.66 y −0.61; el signo negativo indica que la piel tiene menor expresión).
- Acral frente a melanoma cutáneo primario: menor, con efecto en el límite entre mediano y grande según los umbrales del protocolo (delta de Cliff +0.48, IC 95 % 0.31 a 0.63). Esta diferencia **depende de la prueba usada**: con Brunner-Munzel p ajustado = 1.5e-7, y con Kruskal-Wallis + Dunn p ajustado = 0.077.
- CSPG4 también se expresa en piel sana y en otros tejidos normales; 10 de 54 tejidos tienen una mediana mayor que la del grupo acral (medianas de tejido completo, sin prueba estadística).
- Sin asociación detectable con estadio, ulceración, sexo, clúster, pureza tumoral ni fracciones celulares estimadas en las 44 muestras acrales (poder bajo).

## Qué no se puede concluir

La expresión de ARN no demuestra proteína en la superficie celular, accesibilidad del epítopo, unión ni internalización del nanovehículo, ni especificidad tumoral. GTEx no incluye piel acral ni melanocitos aislados. La diferencia entre piel expuesta y no expuesta mezcla sitio anatómico y exposición solar.

## Cómo reproducir

Los datos **no están en el repositorio**. Hay que obtenerlos y colocarlos según `data/README.md`. Después, desde la raíz del proyecto:

```
Rscript scripts/00_preparar_tablas.R
Rscript scripts/05_objetivo1_final.R
Rscript scripts/06_apendice_clinico_acral.R
```

Versión de R y paquetes: ver `entrega_objetivo1/3_reporte_y_figuras/sessionInfo.txt`. El script 05 incluye 17 verificaciones automáticas de auditoría.

## Estructura

| Carpeta | Contenido |
|---|---|
| `scripts/` | Los tres scripts del análisis |
| `docs/` | `DECISIONES.md` (decisiones y su evidencia), `BITACORA.md`, `CAMBIOS_PROTOCOLO_v2.md` |
| `entrega_objetivo1/` | Protocolo y resultados, reporte, figuras y tablas resumen agregadas |

## Datos y privacidad

Este repositorio es público. No se suben datos crudos ni tablas por muestra de pacientes mexicanos hasta contar con el visto bueno de las autoras y autores del estudio original. Los archivos que contienen datos por muestra están excluidos en `.gitignore`.

## Pendientes

- Respuesta de la colaboración (Stef / Dra. Vázquez): matriz de 44 muestras, cuáles de las 65 primarias están en la publicación y pipeline de cuantificación.
- Verificar la referencia completa de la cohorte (Basurto-Lozada et al.): revista, año y volumen.
- Versión de GENCODE de GTEx v8, y fecha de consulta de la identidad del gen en NCBI Gene.
- Aprobación del equipo de los cambios al protocolo.
- Visto bueno para publicar figuras y tablas que muestran valores individuales de la cohorte mexicana.
- Licencia del repositorio: sin decidir (no se agrega archivo LICENSE hasta entonces).

## Fuentes principales
- Basurto-Lozada P, et al. *Nature* 2026;651:221-230. Código de los autores: https://github.com/CGBio-Lab/Mex-acral-exomes-transcriptomes
- Zubair H. protocols.io 2026, doi:10.17504/protocols.io.kqdg3mn6el25/v1 (flujo DESeq2 con TCGA; validado por su autor en TCGA-LUAD).
- Love MI, Huber W, Anders S. *Genome Biol* 2014;15:550. Colaprico A, et al. *Nucleic Acids Res* 2016;44:e71.
- NCI GDC, mRNA Analysis Pipeline: https://docs.gdc.cancer.gov
- Wang Q, et al. *Sci Data* 2018, doi:10.1038/sdata.2018.61.

## Aviso sobre datos
Este repositorio es público. **No se suben datos crudos ni procesados** (ver `.gitignore`). La tabla por muestra de los pacientes mexicanos tampoco se publica hasta contar con el visto bueno de las autoras y autores del estudio original.
