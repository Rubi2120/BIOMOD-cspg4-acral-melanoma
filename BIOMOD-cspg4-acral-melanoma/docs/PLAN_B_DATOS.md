# Plan B: si no llega el acceso a los datos o la matriz no trae CSPG4

Todo lo de abajo viene de descripciones de catálogos y artículos revisados en una búsqueda rápida.
**Cada fuente debe revisarse (plataforma, número de muestras, disponibilidad) antes de usarse.**

| Fuente | Qué es (según su descripción) | Uso posible | Advertencia |
|---|---|---|---|
| GSE5591 (GEO, plataforma GPL4213) | Piel humana normal por sitio anatómico (acral vs no acral), citada en un conjunto de datos de Choi et al. 2023 | **Piel sana acral**: comparar CSPG4 entre piel acral y no acral dentro del mismo estudio | Es microarreglo: se analiza por separado, sin mezclar con RNA-seq. Verificar que CSPG4 tenga sonda |
| GSE162682 (GEO) | RNA-seq de melanoma acral vs no acral en población asiática | Cohorte acral independiente para replicar (no mexicana) | Verificar n, si hay counts y si se puede bajar sin restricción |
| EGAD00001005500 (EGA) | Newell et al.: 87 casos acrales con DNA; 63 con RNA-seq | Cohorte acral grande | Acceso controlado (solicitud al comité) |
| phs001036 (dbGaP) | 34 pacientes con melanoma acral, RNA-seq tumoral | Cohorte acral adicional | Acceso controlado (solicitud) |
| GSE215120 y GSE228421 (GEO) | Célula única de melanoma acral y de piel sana, usadas juntas en un estudio publicado | Ver CSPG4 por tipo celular (melanocitos vs fibroblastos, etc.), lo que corrige la limitación de bulk | Análisis distinto (célula única); evaluar si hay tiempo |
| GSE28909 / GSE28910 | Según su descripción, un primario y una metástasis de melanoma acral | **No sirven** para este análisis (n = 1 por grupo) | Retirar del protocolo |

## Orden sugerido
1. Pedir la matriz completa a Stef o a Ferro si CSPG4 no está en el archivo del Dropbox.
2. Mientras tanto, avanzar con los Niveles B y C (no dependen del acceso).
3. Agregar GSE5591 para «piel sana acral» (análisis aparte).
4. Si se necesita otra cohorte acral, preferir GSE162682 (sin trámite, por confirmar) antes que las controladas.
